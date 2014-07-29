#!/usr/bin/perl
########################################################################
# Script name  :    nph-blast.pl
#
# Date created :    August 2003
#
# Author       :    Shuai Weng <shuai@genome.stanford.edu>
# 
# This script is a simple client interface for displaying the blast 
# search form, running the blast search, and generating the search 
# result. 
#
# 
########################################################################
use strict;
use warnings;
#use diagnostics;

########################################################################
select(STDOUT); 
$| = 1;  # to prevent buffering problems
########################################################################

use CGI;
use CGI qw/:html :form -nph/;
$CGI::POST_MAX=1024 * 100; #max 100k posts
use CGI::Carp qw(fatalsToBrowser);
use Bio::SearchIO;
use MyHTMLResultWriter;
use DBI;
use Bio::GMOD::Blast::Graph;
use Bio::GMOD::Blast::Util;
use Bio::DB::Das::Chado;
use File::Temp qw/tempfile/;
use URI::Escape;
## Browser must have constant ip during session
use CGI::Session '-ip_match';

CGI::Session->name("BLAST-WEB-SESSION"); # set blast web session cookie name
use CGI::Cookie;
use Net::SAML;


## Change this variable to point to your own location and the name 
## for the configuration file
my $CONF_FILE = '../conf/Blast.conf';
### global variables



my $idp = "";
my $embedded = 0;
my $cssurl = "";

my $chadoDb = ""; 
my $chadoDbHost = "";

my $chadoDbUser = "";

my $chadoDbPass = "";



my $title = 'BLAST Search';


########################################################################
my ($datasetDir, $seqtmp, $sequence, $program, $dataset, $options, 
    $filtering);

my (@program, %programLabel, %programType, @db, %dbLabel, %dbType, @matrix);

my ($blastBinDir, $blatBinDir, $blastOutputFile, %port, %host);

my ($imageDir, $imageUrl);


my $needAuthSearchForm = 1;
my $needAuthChado = 1;
my ($sid, $username);

&setVariables; ## read configuration


my $needAuthAll = defined $idp;


### make a new CGI session and generate the cookie
my $q = new CGI();
my $url = $q->url();

### discard old session if it's there
#my $s = CGI::Session->load("driver:db_file", $q) or die CGI::Session->errstr();
#if ( $s->is_expired ) {
#    print $s->header(),
#    $q->start_html(),
#    $q->p("Your session timed out! Refresh the screen to start new session!"),
#    $q->end_html();
#    exit(0);
#}

my $session = new CGI::Session("driver:db_file", $q) or die CGI::Session->errstr();

my $q = $session->query();
## requests a fresh session:
if ($q->param('o') =~ /C|Q|R/) {
    warn "requesting fresh session";
    $session->delete();
    $session->flush();
    $session = new CGI::Session("driver:db_file", $q);    
}

$session->expire('30m');

my $new_session = $session->is_new;
my $cookie = CGI::Cookie->new(-name=>$session->name, -value=>$session->id);
#### store all control paramters in the session

my @names = $q->param();
## store all cgi parameters in session parameters, they wouldn't survive the
## SSO process:
$session->save_param($q);
if ($new_session) {
## make sure the user gets the cookie:
    warn "new session, reloading\n";
    print $q->redirect(-uri=>$q->url(), -cookie=>$cookie); exit;
} else {
##   warn "resuming session ", $session->id(),". featureid: ". $session->param('feature_id');
}

($sid, $username) = &doSSO(); 
### SSO killed our parameters, so load them again
$session->load_param($q);
#### ICICIC check for valid session, keep SAML and CGI sesssions synched

$session->param('samlid', $sid);
$session->param('username', $username);

$embedded = $session->param('embedded');

if (_get_param('feature_id') && _get_param('start')  
    && _get_param('end') && _get_param('locus') ) {
    
    &writeCoordsToDb;

    }
    elsif (! $q->param('sequence') && (! $q->param('file'))) {

    &printSearchForm;
   
}
else {

    &checkArgsAndDoSearch;
  

}
$session->flush(); 
exit(0);





####################################################################
sub printSearchForm {
####################################################################

    &setVariables;

    &printStartPage;

    

    print &blastForm('#CCCCFF');
   
    print end_html;

}


sub writeCoordsToDb {

### do not do anything, unless we have a valid session with a session id:
    unless ($session && ! $session->is_expired && $session->param('samlid') == $sid) {
	warn "attempt to write coords with invalid session\n";
	exit(0);
    }
### we need a valid DB connection first:
  
    my $dbh;
    if ($chadoDb) {
	$dbh = DBI->connect( "DBI:Pg:database=$chadoDb;host=$chadoDbHost", 
			     $chadoDbUser, $chadoDbPass, { 
				 RaiseError => 1,
			     }  
	    ) or die "Unable to connect to the chado DB: $chadoDb !\n $! \n $@\n$DBI::errstr";
	
    } else {
	warn "chadoDb is required for db updates, please define it in the config file!\n";
	exit (0);
    };
    
### yes, do it, user issued a confirm!
    if (_get_param('confirm') eq 'true'){
	my $testsql="SELECT COUNT(*) FROM featureloc where feature_id="._get_param('feature_id')."";
	my $count = $dbh->selectrow_array($testsql);
	if ($count == 0) {
	    my $sql=
		"INSERT into chado.featureloc(feature_id, srcfeature_id, 
                 fmin,fmax,strand) VALUES ("._get_param('feature_id').", 
                 (SELECT feature_id from feature where name='"._get_param('locus')."'),"
                 ._get_param('start')."-1,"._get_param('end')."-1, "._get_param('strand')."  )";
	    my $prep = $dbh->prepare($sql) or die $dbh->errstr;
	    	      
	    $prep->execute(); 
	    $prep->finish();

	    $session->clear(['locus', 'start', 'end', 'confirm', 'strand']);
	    print $session->header();
	    print start_html();
	    print qq| <script type="text/javascript">
            alert("Location saved");
	         window.location.replace("$url");
                </script> |;
	    print end_html();
	}
	else {
	    print start_html();
	    print ' <script type="text/javascript">
		 alert("This feature has already a location. You have to modify the location manually ");
		 window.location.replace("'.$url.'");
		 </script>';
	    $session->clear(['locus', 'start', 'end', 'confirm', 'strand']);
	    print end_html();
	}
    }
    ### ask user for confirmation in a dialogue
    elsif (! _get_param('confirm')) {
	
	my $rul = $q->url();
	print start_html();
	print '<script type="text/javascript">
<!--
var confirmurl = "'.$rul.'?confirm=true";
var returnurl = "'.$rul.'?confirm=nope";
;
var answer = confirm ("Are you sure you want to save these coordinates permanently in the database? This can only be done once. Source feature: ' 
._get_param('locus') . ' start= '._get_param('start').', end = '._get_param('end').', strand = '._get_param('strand').
' ");

if (!answer){
 window.location.replace(returnurl);
}else{        
 
 window.location.replace(confirmurl);
	
}
// -->
</script> ';
	print end_html();
 
    } else {
	
	### user wanted to break out, so delete the parameters
	$session->clear(['locus', 'start', 'end', 'confirm']);
	
	print ($q->redirect(-uri=>$q->url(), -cookie=>$cookie));
       
    }
    exit (0);
}









####################################################################
sub checkArgsAndDoSearch {
####################################################################
    &setVariables;

    &printStartPage;

    #####################################################
    

    &checkParameters;

    #####################################################
    &createTmpSeqFile;

    &setOptions;

    #####################################################
    &runBlast;
    
    #####################################################
    &showGraph;

    #####################################################
    &showResult;
    
    #####################################################
 
    unlink($seqtmp);

    unlink($blastOutputFile);

    print end_html;

}

####################################################################
sub runBlast {
####################################################################

     
    my  ($fh, $tmpfile) = tempfile(UNLINK => 1);
    unless ($program =~ /^blat/i) {

	if ($filtering) {

	    open(OUT, ">$blastOutputFile") ||
		die "Can't open '$blastOutputFile' for writing:$!";

	    print OUT "Filtering On\n";
	
	    close(OUT);

	}

    }

    my $cmd;

    if ($program =~ /^(blat|tblat)/i) {

	my $port = $port{$program}{$dataset};

	my $host = $host{$program}{$dataset};

	if ($program =~ /^blat/i) {

	    $cmd = "$blatBinDir/gfClient $host $port / $seqtmp -out=blast $blastOutputFile >> /dev/null 2>&1";
	
	}
	else {

	    $cmd = "$blatBinDir/gfClient $host $port / $seqtmp -out=blast -q=prot -t=dnax $blastOutputFile >> /dev/null 2>&1"

	}


    }
    else {
       
       
	$program = $blastBinDir.$program;
        $cmd = "$program -db $dataset -query $seqtmp  -out $blastOutputFile  $options > $tmpfile  2>&1";
 #	die $cmd;
	
   #	$cmd =" $program $dataset $seqtmp $options -cpus=2 -progress=20 >> $blastOutputFile 2>&1";
   #    die "Erreur blast :$cmd";
    }

    my $err = system($cmd);

    if ($err) {

	print "Error occurred when running BLAST/BLAT program. See the following message:", p;

	print '<pre>';

	while (<$fh>) {
	    print;
	}
	close ($fh);
	print '</pre>';

	#$session->delete();
	#$session->flush();
	
	exit;

    }
   

}

####################################################################
sub showGraph {
####################################################################

   my $graph;
 
   eval {
       $graph = Bio::GMOD::Blast::Graph->new(-outputfile=>$blastOutputFile,
					     -dstDir=>$imageDir,
					     -dstURL=>$imageUrl);
   ;
       if (ref $graph) {
	   $graph->showGraph;
       }
   };
   if ($@) {
       print p("error parsing blast output");
   }
   
}

####################################################################
sub showResult {
####################################################################

    print p, hr;
    my $q = new CGI;

    my $in = new Bio::SearchIO(-format=>'blast',
			       -file=>$blastOutputFile,
			       -verbose=>0,
			       -signif=>1);

    my $result = $in->next_result;
    return unless $result->num_hits > 0;
    ### check if the feature already has a feature_loc
    ### then we do not display the buttons:
    my $fid = undef;
    my $dbh;
    if ($chadoDb && _get_param('feature_id')) {
	$dbh = DBI->connect( "DBI:Pg:database=$chadoDb;host=$chadoDbHost", 
			     $chadoDbUser, $chadoDbPass, { 
				 RaiseError => 1,
			     }  
	    ) or die "Unable to connect to the chado DB: $chadoDb !\n $! \n $@\n$DBI::errstr";
    
	if ($dbh) {	
	    my $testsql="SELECT COUNT(*) FROM featureloc where feature_id="._get_param('feature_id')."";
	    my $count = $dbh->selectrow_array($testsql);
	    if ($count == 0) {
		$fid = _get_param('feature_id');
	    } else {
		$session->clear('feature_id');
	    }
	 }   
    }
    ### buttons will only be displayed if feature_id is provided
    my $writer =  MyHTMLResultWriter->new($fid); 
   
    my $out = new Bio::SearchIO(-writer => $writer);

    
    eval { $out->write_result($result); };
    #$session->delete();
    #$session->flush();


    if ($@) {

	#die $@  # "<pre> error=$@</pre> ",p;

    }

}

#####################################################################
sub blastForm {
#####################################################################
# This method simply calls 'blastSearchBox' method to display the 
# popup menu for databases and search programs, and calls 
# 'blastSearchOptions' to display the options for the blast program. 

    my ($optionBg) = @_;

    $optionBg ||= 'white';

    return 

	start_multipart_form(-action=>$url).
	table({-width=>600,
	       -border=>0,
	       -rules=>'none',
	       -cellpadding=>0,
	       -cellspacing=>0},
	      Tr(td({-colspan=>2},
		 )).
	      Tr(td({-colspan=>2},
		    &blastSearchBox)).
	      Tr(td({-colspan=>2},
		    &blastSearchOptions($optionBg)))).
		    end_multipart_form();


}

#####################################################################
sub blastSearchBox {
#####################################################################
# This method displays the popup menu for the blast search databases
# and search programs. You should update this method to include your
# database names and search programs.
#DAMIENGENESTE
    my $name;
    my $value;
    my $buffer;
    my @pairs;
    my $pair;
    my %tab=();
    my $residue="";
    my $db;
    my $dbh;
    if ($chadoDb) {
	$dbh = DBI->connect( "DBI:Pg:database=$chadoDb;host=$chadoDbHost", 
			     $chadoDbUser, $chadoDbPass, { 
				 RaiseError => 1,
			     }  
	    ) or die "Unable to connect to the chado DB: $chadoDb !\n $! \n $@\n$DBI::errstr";
	
	$db  = Bio::DB::Das::Chado->new(
	    -dsn  => "DBI:Pg:database=$chadoDb;host=$chadoDbHost",
	    -user => $chadoDbUser,
	    -pass => $chadoDbPass,
	    -srcfeatureslice => 1,
	    -recursivMapping => 1,
	    -tripal => 1
	    );
    };


    my $displayName;
    my $description;
    my $primaryTag;

#Get residues an add automatic annotation to the query comment if available
    if ($dbh && $db && (my $id = _get_param('feature_id'))){

	    my $feature = $db->get_feature_by_id($id);
	    #die ("feature $id not found") unless $feature;
	    if (ref $feature && eval { ref $feature->seq() } ) {
		$residue = $feature->seq()->seq 
		    || die "missing sequence for feature ".$feature->display_name;
		$displayName = $feature->display_name();
		($description) =  $feature->get_tag_values('note');
		$primaryTag = $feature->primary_tag();
	    } else {
#Requete to get residues
		my $prep = $dbh->prepare("SELECT residues FROM chado.feature where feature_id="._get_param('feature_id')."") or die $dbh->errstr;
		$prep->execute() or die "Query failed\n"; 
		
		my @ar = $prep->fetchrow_array();
		$residue = $ar[0];
		
		
		$prep->finish();

		my $prep = $dbh->prepare("SELECT name FROM chado.feature where feature_id="._get_param('feature_id')."") or die $dbh->errstr;
		$prep->execute() or die "Query failed\n"; 		
		my @ar = $prep->fetchrow_array();
		$displayName = $ar[0];		
	    }  
    }



    return b("Query Comment (optional, will be added to output for your use):").br.
	textfield(-name=>'seqname', -value=>join(" ", ($displayName, $primaryTag, $description)),
		  -size=>60
	).p."\n".
	b(font({-color=>'red'},
	       "NOTE: If the input sequence is less than 30 letters you should change the default Cutoff Score value to something less than 100 or you can miss matches.")).p.
	       b("Upload Local TEXT File: FASTA, GCG, and RAW sequence formats are okay").br.
	       "WORD Documents do not work unless saved as TEXT.".
	       filefield(-name=>'filename',
			 -size=>60,
			 -accept=>"text/*").p."\n".
			 b("Type or Paste a Query Sequence : (No Comments, Numbers Okay)").br."\n".
			 textarea(-name=>'sequence',
				  -columns=>'80',
				  -rows=>'5',
				  -value=>"\n$residue").p."\n".
				  b("Choose the Appropriate Search Program:").br.
				  popup_menu(-name=>'program',
					     -values=>\@program,
					     -labels=>\%programLabel).p."\n".
					     b("Choose a Sequence Database:").br.
					     popup_menu(-name=>'database',
							-values=>\@db,
							-labels=>\%dbLabel).br."\n".
							submit(-value=>'Run BLAST').' or '.reset();


}


######################################################################
sub blastSearchOptions {
######################################################################
# This method is used to display the options for the blast search.

    my ($optionBg) = @_; 

    my $ncbiBlastHelp = 'http://www.ncbi.nlm.nih.gov/books/NBK1763/';

    ######## output format 
    my @format = ('gapped', 'nongapped');

    my %formatLabel = ('gapped'=>'gapped alignments',
		       'nongapped'=>'nongapped alignments');

    my $formatDefault = 'gapped';


     ######## taskt 
    my @task = ('','blastn-short', 'megablast','dc-megablast','blastp-short');

    my %taskLabel = (''=>'',
	             'blastn-short'=>'BLASTN program optimized for sequences shorter than 50 bases',
		       'megablast'=>'Traditional megablast used to find very similar (e.g., intraspecies or closely related species) sequences',
	                 'dc-megablast'=>'Discontiguous megablast used to find more distant (e.g., interspecies) sequences',
	                  'blastp-short' =>'BLASTP optimized for queries shorter than 30 residues');

    my $taskDefault = "";   

    

    ######## cutoff score 
    my @cutoff = ('default', '30', '50', '70', '90', '110');

    my $cutoffDefault = 'default';

    ######## word length 
    my @wordLength = ('default');

    for (my $i = 15; $i >= 2; $i--) {
	push(@wordLength, $i);
    }

    my $wordLengthDefault = 'default';
    ######## Perc_ identity
    my @percident = ('default','50','60','70','80','90','95');

    my $P_ident_Default = 'default';

    ######## Expect threshold
    my @eValue = ('default', '0.0001', '0.01', '1', '10', '100', '1000');
    
    my $eValueDefault = 'default';

    ######## Number of best alignments to show
    my @alignNum = ('0', '25', '50', '100', '200', '400', '800', '1000');

    my $alignNumDefault = 100;
    
    ######## Sort output by 
#    my @sortBy = ('pvalue', 'count', 'highscore', 'totalscore');

#    my $sortByDefault = 'pvalue';

    return b('Options: ').
	   'For descriptions of BLAST options and parameters, refer to the '. 
	   a({-href=>$ncbiBlastHelp}, 
	     'BLAST documentation at NCBI.').p.
	   table({-bgcolor=>$optionBg,
		  -cellspacing=>0},
		  Tr(th({-align=>'left'},
		       'Output format :').
		    td(popup_menu(-name=>'output',
				  -values=>\@format,
				  -default=>$formatDefault,
				  -labels=>\%taskLabel)).
		    td(br)).
                     Tr(th({-align=>'left'},
		       'Task :').
		    td(popup_menu(-name=>'task',
				  -values=>\@task,
				  -default=>$taskDefault,
				  -labels=>\%formatLabel)).
		    td(br)).
	#	    Tr(th({-align=>'left'},
	##	       'Comparison Matrix :').
	#	    td(popup_menu(-name=>'matrix',
       #    				  -values=>\@matrix)).
	#	    td(br)).
		 Tr(th({-align=>'left'},
		       'Cutoff Score (S value) :').
		    td(textfield(-name=>'sthr',
				  
				  -default=>0)).
		    td(br)).
		 Tr(th({-align=>'left'},
		       'Word Length (W value) :').
		    td(popup_menu(-name=>'wordlength',
				  -values=>\@wordLength,
				  -default=>\$wordLengthDefault)).
		    td('Default = 11 for BLASTN, 3 for all others')).
		 Tr(th({-align=>'left'},
		       '% Identity  :').
		    td(textfield(-name=>'perc_ident',
				  
				  -default=>0)).
		    td(br)).
		 Tr(th({-align=>'left'},
		       'Expect threshold (E threshold) :').
		    td(textfield(-name=>'ethr',				  
				  -default=>0)).
		    td(br)).
	#	 Tr(th({-align=>'left'},
	#	       'Number of best alignments to show :').
	#	    td(popup_menu(-name=>'showal',
	#			  -values=>\@alignNum,
	#			  -default=>$alignNumDefault)).
	#	    td(br)).
		 Tr(th({-align=>'left'},
		       'Others options :').
		    td(textfield(-name=>'options_bis',
			         -values=>"")).
		    td('-options=value ')));

#		 Tr(th({-align=>'left'},
#		       'Sort output by :').
#		    td(popup_menu(-name=>'sortop',
#				  -values=>\@sortBy,
#				  -default=>\$sortByDefault))));

}

####################################################################
sub setVariables {
####################################################################

    open(CONF, "$CONF_FILE") ||
	die "Can't open '$CONF_FILE' for reading:$!";

    while(<CONF>) {
	
	if (/^\#/ || /^$/) { next; }

	chomp;

	my ($name, $value);
	
	if (/^([^\=]+) *= *(.+)$/) {

	    $name = $1;

	    $value = $2;
	    $value =~ s/^ *//;
	    $value =~ s/ *$//;

	}

	if ($name =~ /^tmpDir/i) {

	    $seqtmp = $value."blastseq.$$.tmp";
    
	    $blastOutputFile = $value."blast.$$.output";

	}
	elsif ($name =~ /^imageDir/i) {

	    $imageDir = $value;

	}
	elsif ($name =~ /^imageUrl/i) {

	    $imageUrl = $value;

	}
	elsif ($name =~ /^databaseDir/i) {

	    $datasetDir = $value;

	    $ENV{'BLASTDB'} = $value;

	}
	elsif ($name =~ /^database/i) {

	    my ($db, $type, $desc) = split(/=>/, $value);

	    push(@db, $db);

	    $dbType{$db} = $type;

	    $dbLabel{$db} = $desc;

        }
	elsif ($name =~ /^blastBinDir/i) {

	    $blastBinDir = $value;

	}
	elsif ($name =~ /^blatBinDir/i) {

	    $blatBinDir = $value;

	}
	elsif ($name =~ /^port/i) {

	    my ($port, $host, $program, $dataset) = split(/=>/, $value);

	    $port{$program}{$dataset} = $port;

	    $host{$program}{$dataset} = $host;

	}
	elsif ($name =~ /^program/i) {

	    my ($program, $type, $desc) = split(/=>/, $value);

	    push(@program, $program);

	    $programType{$program} = $type;

	    $programLabel{$program} = $desc;

	}
	elsif ($name =~ /^blastMAT/i) {

	    $ENV{'BLASTMAT'} = $value;

	}
	elsif ($name =~ /^matrix/i) {

	    push(@matrix, $value);

	}
	elsif ($name =~ /^blastFILTER/i) {
 
	    $ENV{'BLASTFILTER'} = $value;

	}elsif ($name =~ /^idp/i) {
	    $idp =  uri_escape($value);

       	}elsif ($name =~ /^chadoDb(\s)*$/i) {
	
	    $chadoDb = $value;

	}
	elsif ($name =~ /^chadodbhost/) {
	    $chadoDbHost = $value;

	}elsif ($name =~ /^chadoDbUser/i) {
	    
	    $chadoDbUser = $value;
	}elsif ($name =~ /chadoDbPass/i) {
	    $chadoDbPass = $value;
	}elsif ($name =~ /cssurl/i) {
	    $cssurl = $value;
	} else {
	   # warn " unknown configuration option $name ";
	}
    }
    close(CONF) or die "error closing conf";

    ############ 
    $program = _get_param('program');
    $dataset = _get_param('database');    

}

####################################################################
sub checkParameters {
####################################################################

    &checkDatabase;
    
    &checkSequence;

    &checkSeqLengthAndSvalue;

    &checkSeqLengthAndWordLength;
    
    &checkDatasetAndProgram;

    &checkProgramAndSeqlength;

}

####################################################################
sub createTmpSeqFile {
####################################################################

    my $seqname = _get_param('seqname');

    $seqname .= "  (Length: ".length($sequence).")";

    Bio::GMOD::Blast::Util->createTmpSeqFile($seqtmp, $seqname, $sequence);

}

####################################################################
sub checkSeqLengthAndSvalue {
####################################################################
    
    if (_get_param('sthr') !=0  && _get_param('sthr') >length($sequence) ) {

	print "The Cutpff Score is higher than the  sequence length ", p,
	       "Return to the form to adjust either the cutoff score value or ",
               " sequence.",p;
    
	print end_html;

	exit;

    }

}

####################################################################
sub checkTextfield {
####################################################################
    
    my $text  =$_[0];
    if ($text=~ /;/  || $text=~/&/ ) {

	print  "Error in some textfield  ", p,
	"Return to the form to adjust it";
    
	print end_html;

	exit;

    }

}

####################################################################
sub checkSeqLengthAndWordLength {
####################################################################

    if ($program eq "blastn" && _get_param('wordlength') ne "default" && 
	_get_param('wordlength') < 11 && length($sequence) > 10000) {

	print "The maximum sequence length for a word length of less than 11 is ", b("10000"), ".", p,
	      "Return to the form to adjust either the word length ",
	      "or sequence.",p;
	
	print end_html;

	exit;
    }

}


####################################################################
sub checkDatasetAndProgram {
####################################################################

    if ($dbType{$dataset} eq $programType{$program}) {

	return;

    }

    #### add your own checking code here to make sure 
    #### the selected blast search program matches the database.

    print "Your choice of Database (".b($dataset).") does not match the ",
          "choice of BLAST search program (".b($program).").",p,
	  "BLASTP and BLASTX require a protein sequence database and ",
	  "other BLAST programs require a nucleotide  sequence database. ",p,
	  "Return to the form and adjust either the program ",
          "or database selection.",p;
    
    print end_html;

    exit;

}

####################################################################
sub checkProgramAndSeqlength {
####################################################################
    
#    if (!param('email') && $program =~ /(tblastx|tblastn)/ &&
#	length($sequence) > 5001) {

#	print "The maximum sequence length for TBLASTN and TBLASTX ",
#              "is 5,000 bp unless the Email option is used.",p,
#	      "Return to the form and reduce the sequence length, ",
#	      "select the Email option or choose another BLAST program.";

#        print end_html;
    
#        exit;

#    }

}

####################################################################
sub checkEmail {
####################################################################
    
    if (!_get_param('email')) { return; }

    if (!Bio::GMOD::Blast::Util->validateEmail(_get_param('email'))) {

	print "You requested that the results be sent to your e-mail ",
	      "account. However, your email address is missing, ",
	      "appears incomplete, or does not contain a valid hostname.",p,
	      "You entered this email address: ".b(_get_param('email')),p,
	      "Please return to the form and check that your email address ",
	      "is correct.",p;

	print end_html;
    
	exit;

    }

}

####################################################################
sub setOptions {
####################################################################
    
    return if ($program =~ /blat/i);

   # $options = &blastOptions($program, length($sequence));
    $options="";
    my @params=(_get_param('sortop'),_get_param('ethr'),_get_param('sthr'),_get_param('output'),_get_param('wordlength'),_get_param('perc_ident'),_get_param('options_bis'));
    foreach my $text( @params){
    
	&checkTextfield($text); 
	
    }
    if (_get_param('sortop') && _get_param('sortop') ne "pvalue") { 

	$options .= " -sort_by_"._get_param('sortop'); 

    }
    if (_get_param('ethr') && _get_param('ethr') ne "default") {

	$options .= " -evalue="._get_param('ethr');

    }
    if (_get_param('sthr') && _get_param('sthr') ne "default") {

	$options .= " -min_raw_gapped_score="._get_param('sthr')." ";

    }
 #   $options .= " B="._get_param('showal')." V="._get_param('showal');

    if (_get_param('output') ne "gapped") { $options .= " -ungapped"; }
    
    if (_get_param('task') ne "" ){ $options.="-task=\""._get_param('task')."\"";}

 #   if ($program ne "blastn" && _get_param('matrix') ne "BLOSUM62") { 

#	$options .= " -matrix="._get_param('matrix'); 

  #  }

    if (_get_param('wordlength') ne "default") { 

	$options .= " -word_size="._get_param('wordlength');

    }
    if (_get_param('perc_ident') ne "default") { 

	$options .= " -perc_identity="._get_param('perc_ident');

    }

    if (_get_param('options_bis') ne "") {
    
	$options.=" "._get_param('options_bis')." ";
    }


}

#######################################################################
sub blastOptions {
#######################################################################
    my ($program, $seqlen) = @_;

    return if ($program =~ /^blat/i);

    my $hspmax;
    my $gapmax;

    if ( $seqlen < 10000 ) {

	if ($program eq "blastn") {

	    $hspmax = 6000;

	    $gapmax = 3000;

	} 
	else {

	    $hspmax = 2000;

	    $gapmax = 1000;

	}
    } 
    else {

	$hspmax = 10000;

	if ($program eq "blastn") {

	    $gapmax = 3000;

	} 

	else {

	    $gapmax = 1000;
	}
    }

    return " -hspsepsmax=" . $hspmax . " -hspsepqmax=" . $hspmax . " -gapsepsmax=" . $gapmax . " -gapsepqmax=" . $gapmax . " ";

}
		
####################################################################
sub checkDatabase {
####################################################################
    
    if (_get_param('database') eq '-') {

	print b("No Database Selected."),p;

	print "Please return to the form and select a database.",p;

	print end_html;

	exit;

    }

    my $datasetLockFile = $datasetDir._get_param('database').".update";

    if (-e "$datasetLockFile") {

	print b("SORRY the "._get_param('database')." dataset is currently being UPDATED. Please try again in a few minutes or select another dataset."),p;
	
	print end_html;

	exit;

    }

}

####################################################################
sub checkSequence {
####################################################################

    my $filehandle = $q->upload('filename');

    if ($filehandle) {

	while (<$filehandle>) {

	    $sequence .= $_;
	}
       
    }
    else { $sequence = $q->param('sequence'); }

    Bio::GMOD::Blast::Util->deleteUnwantedCharFromSequence(\$sequence);

    if (!$sequence) {

	print b("No Sequence Provided."),p;
	  
	print "Please return to the form and enter a sequence.",p;

	print end_html;

	exit;

    }

}

####################################################################
sub printStartPage {
####################################################################
   
    print $session->header;

 
    print start_html(-title=>$title, -style=>{'src'=>$cssurl});
    print '<a name="pagetop" />';

    unless ($embedded) {
	print center(h2({-style=>"margin:20px"}, $title)), hr;
       
	if ($sid) {
	    ## we are logged in!
	    my $conf = "PATH=/var/zxid/&URL=$url";
	    my $cf = Net::SAML::new_conf_to_cf($conf);
	    
	    my $ses = Net::SAML::fetch_ses($cf, $sid);
	    my $cgi = Net::SAML::new_cgi($cf, $q->query_string());
	    my $redir =  Net::SAML::sp_slo_redir($cf, $cgi ,$ses); 
	    $redir =~ s/^Location: //;
	    $redir =~ s/[\r|\n]//g;
	    print p("Welcome: $username, you are logged in.",
	    a({-href=>$redir}, "Logout") );
	    #print p(Net::SAML::fed_mgmt_cf($cf, undef, -1, $sid, 0x1900));
	}
	
    }

    


    if (_get_param('program') =~ /^(blat|tblatn)$/i) {

 	print center('If there are no hits found using BLAT/TBLATN'.br.'the remainder of this page will be blank.'), hr({-width=>'20%'});

    }


    
}

####################################################################
sub doSSO {
####################################################################
    require Net::SAML;




    my $conf = "PATH=/var/zxid/&URL=$url";
    my $cf = Net::SAML::new_conf_to_cf($conf);
    my $qs = $ENV{'QUERY_STRING'};
    my $q = new CGI;

    if ($qs =~ /o=P/) {
	$qs = $q->query_string();
    } else {     		
	$qs .= "&e=$idp&l0=TRUE"
	    unless ($qs =~ /s=|(SAMLart=)|o=B/);
    };
    my $res = Net::SAML::simple_cf($cf, -1, $qs, undef, 0x1828); # keep the flags 0x1828 !!! 
    my ($redirecturl) = $res =~ /^Location:\s+(.+)/;
    $redirecturl =~ s/\r|\n|\s//g;
    print STDERR "RESULT: $res\n,$redirecturl\n";
    my $op = substr($res, 0, 1);
    if ($op eq 'L') { die "$res" unless $redirecturl;
		      print ($q->redirect($redirecturl)); 
		      exit } # LOCATION (Redir) or CONTENT
    if ($op eq 'C') { print ($res); exit }
    if ($op eq 'n') { exit; } # already handled
    if ($op eq 'e') { #die "cannot choose IDP"; exit;
 } # not logged in
    elsif ($op ne 'd') { die "Unknown Net::SAML::simple() res($res)"; }
    my ($sid) = $res =~ /^sesid: (.*)$/m;  # Extract a useful attribute from SSO output
    die "invalid session id" unless $sid;
    my ($username) =   $res =~ /^displayName: (.*)$/m ;

    return ($sid, $username);
}

####################################################################
sub renderLogOutButton {

}
####################################################################



####################################################################
# sub writeLog {
####################################################################
#    my ($Cuser, $Csystem) = @_;

#    if (_get_param('email')) { return; }
 
#    Bio::GMOD::Blast::Util->writeLog($program, $dataset, $options, 
#                   $Cuser, $Csystem, length($sequence), $remoteLink);
		      
# }

####################################################################

sub _get_param {
    my $p = shift;
    return unless $q || $session;
    return ($q->param($p) || ($session) ?  $session->param($p) : undef);

}























