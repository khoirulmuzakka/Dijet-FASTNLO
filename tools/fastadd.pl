#!/usr/bin/env perl 
#
# fastNLO table addition script:
#    This perl script creates either a total sum table (scen.tab) or
#    sum tables for statistical uncertainty evaluation: 
#    - For each single LO table with *one* 'dummy' NLO table
#    - For each single NLO table with *all* LO tables or
#    a summary warmup table from multiple warmup runs
#
# Version:
# 
# created by K. Rabbertz: 05.03.2006
# last modified:
# 02.05.2011 KR: Added warmup table part
#
#-----------------------------------------------------------------------
# Todo:
#

use Cwd;
use English;
#use FindBin qw($Bin);
#use lib "$Bin/modules";
use Getopt::Std;
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n##################################################\n";
print "# fastadd.pl: Starting table addition for fastNLO: FASTADD_$date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h, $opt_l, $opt_n, $opt_s, $opt_v, $opt_w ) = ( "", "", "", "", "", "1", "" );
getopts('dhl:n:sv:w') or die "fastadd.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "fastadd.pl\n";
    print "Usage: fastadd.pl [switches/options] scenario\n";
    print "  -d              Verbose output\n";
    print "  -h              Print this text\n";
    print "  -l dir          Directory for LO tables, (def.=scenario)\n";
    print "  -n dir          Directory for NLO/NNLO tables, (def.=scenario)\n";
    print "  -s              Produce tables for statistical evaluation,\n".
	"                  i.e. combinations of each LO with 1 NLO table and\n".
	"                  all LO with each NLO table\n";
    print "  -v #            Choose between fastNLO version 1, 2 or 2.1 (def.=1)\n";
    print "  -w              Write summary table from multiple warmup runs (def.=no)\n\n";
    exit;
}

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fastadd.pl: Error! Need one scenario name!\n";
}
my $scen   = shift;
my $debug  = $opt_d;
my $vers   = $opt_v;

#
# Initialization
#
my $lodir    = "${scen}"; 
if ( $opt_l ) {$lodir = $opt_l;}
my $nlodir   = "${scen}"; 
if ( $opt_n ) {$nlodir = $opt_n;}
my $nnlodir  = $nlodir; 
my $wrmdir   = "${scen}wrm"; 
if ( $opt_n ) {$wrmdir = $opt_n;}
my $tabext = "raw";
if ($vers == 2) {$tabext = "tab";} 
my $loglob   = "${scen}*born*.${tabext}*";
my $nloglob  = "${scen}*nlo*.${tabext}*";
my $nnloglob = "${scen}*thrcor*.${tabext}*";
my $wrmglob = "${scen}*wrm*.dat";
my $wrmglobdef = "fastNLO-warmup*.dat";

# Directory
my $sdir = getcwd();

#
# Check on nlofast-add resp. fnlo-merge
#
my $merger = "nlofast-add";
my $cmd;
unless ( $opt_w ) {
    if ( $vers == 2 ) {
#	$merger = "fnlo-merge";
	$merger = "fnlo-merge2";
    }
    $cmd = `which ${merger}`;
    chomp $cmd;
    unless ( $cmd ) {
	if ( $ENV{NLOJET} ) {
	    if ( -f "$ENV{NLOJET}/bin/${merger}" ) {
		$cmd = "$ENV{NLOJET}/bin/${merger}";
	    } else {
		die "fastadd.pl: ERROR! ${merger} command not found, ".
		    "neither via \`which\` nor in $ENV{NLOJET}/bin, aborted!\n";
	    }
	} else {
	    die "fastadd.pl: ERROR! ${merger} command not found ".
		"via \`which\` and \$NLOJET is not set, aborted!\n";
	}
    }	
}

#
# Analyze warmup tables
#
if ( $opt_w ) {
    print "\nfastadd.pl: Warmup mode\n";
    if ( -d "$wrmdir" ) {
	chdir $wrmdir;
	my @files = glob $wrmglob;  
	chomp @files;
	unless ( @files ) {
	    print "fastadd.pl: Warning! No warm-up files found for scenario $scen,\n";
	    print "            looking for generic filename fastNLO-warmup\* instead.\n";
	    @files = glob "fastNLO-warmup*";  
	    chomp @files;
	    unless ( @files ) {
		die "fastadd.pl: Warning! No warm-up files found, stopped\n";
	    } 
	} 
	my $ifil = 0;
	my $nent = 0;
	my $stat = 0;
	my @xmin;
	my @blow;
	my @bhig;
	my @clow;
	my @chig;
	foreach my $file ( @files ) {
	    open(INFILE,"< $file") or die "fastadd.pl: Error! Could not open $file!\n";
	    my $ient = 0; 
	    if ( $debug ) {print "fastadd.pl: Analyzing file no.: $ifil\n";}
	    while ( my $in = <INFILE> ) {	
		if ( $in =~ "//" ) {
		    my $tmp = $in;
		    chomp $tmp;
		    my @tmps = split(" ",$tmp);
		    $stat = $stat + $tmps[1];
		    if ( $debug ) {print "fastadd.pl: Accumulated statistics: $stat\n";}
		} else {
		    if ( $vers == 2 ) {
			my $tmp = $in;
			chomp $tmp;
			my @tmps = split(" ",$tmp);
			if ( !$ifil || $tmps[4] < $xmin[$ient] ) {$xmin[$ient] = $tmps[4];} 
			if ( !$ifil || $tmps[10] < $blow[$ient] ) {$blow[$ient] = $tmps[10];} 
			if ( !$ifil || $tmps[16] > $bhig[$ient] ) {$bhig[$ient] = $tmps[16];} 
			$ient++;
		    } else {
			my $tmp = $in;
			chomp $tmp;
			my @tmps = split(" ",$tmp);
			if ( !$ifil || $tmps[4] < $xmin[$ient] ) {$xmin[$ient] = $tmps[4];} 
			if ( !$ifil || $tmps[10] < $blow[$ient] ) {$blow[$ient] = $tmps[10];} 
			if ( !$ifil || $tmps[16] > $bhig[$ient] ) {$bhig[$ient] = $tmps[16];} 
			if ( !$ifil || $tmps[22] < $clow[$ient] ) {$clow[$ient] = $tmps[22];} 
			if ( !$ifil || $tmps[28] > $chig[$ient] ) {$chig[$ient] = $tmps[28];} 
			$ient++;
		    }
		}
	    }
	    if ( !$ifil ) {$nent = $ient};
	    if ( $debug ) {print "fastadd.pl: Number of entries: $nent\n";}
	    if ( $ifil && $ient != $nent ) {
		print "fastadd.pl: Error! Inconsistent line numbers in warm-up files found \n";
		die "            (ient = $ient, nent = $nent), stopped\n";
	    }
	    close INFILE;
	    $ifil++;
	}
	my $outfile = "${scen}wrm.dat";
	open(OUTFILE,"> $outfile") or die "fastadd.pl: Error! Could not open $outfile!\n";
	my $line = "      // $stat contributions (!= events) in warm-up run\n";
	print OUTFILE $line;
	for (my $ient = 0; $ient < $nent; $ient++) {
	    my $line;
	    if ( $vers == 2 ) {
		$line = "      xlim[ $ient ] = $xmin[$ient] , mulo[ $ient ] = $blow[$ient] , muup[ $ient ] = $bhig[$ient] ;\n";
	    } else {
		$line = "      xlim[ $ient ] = $xmin[$ient] , scale1lo[ $ient ] = $blow[$ient] , scale1hi[ $ient ] = $bhig[$ient] , scale2lo[ $ient ] = $clow[$ient] , scale2hi[ $ient ] = $chig[$ient];\n";
	    }
	    print OUTFILE $line;
	}
	close OUTFILE;
	print "fastadd.pl: Finished table addition, result file is: $outfile\n";
	chdir $sdir;
	exit(0);
    }
}

#
# Find LO tables
#
my @lotabs;
if ( -d "$lodir" ) {
    chdir $lodir;
    @lotabs = glob $loglob;
    chdir $sdir;
}
unless ( @lotabs ) {
    print "fastadd.pl: No LO table found in $lodir, now looking in $sdir ...\n";
    $lodir  = ".";
    @lotabs = glob $loglob;
} 
unless ( @lotabs ) {
    die "fastadd.pl: ERROR! No LO table found, aborted!\n";
}
my $ntab = scalar @lotabs;
print "fastadd.pl: $ntab LO tables found.\n";
if ( $debug ) {print "fastadd.pl: DEBUG! lotabs @lotabs\n";}

#
# Find NLO tables
#
my @nlotabs;
if ( -d "$nlodir" ) {
    chdir $nlodir;
    @nlotabs = glob $nloglob;
    chdir $sdir;
}
unless ( @nlotabs ) {
    print "fastadd.pl: No NLO table found in $nlodir, now looking in $sdir ...\n";
    $nlodir  = ".";
    @nlotabs = glob $nloglob;
} 
unless ( @nlotabs ) {
    die "fastadd.pl: ERROR! No NLO table found, aborted!\n";
}
$ntab = scalar @nlotabs;
print "fastadd.pl: $ntab NLO tables found.\n";
if ( $debug ) {print "fastadd.pl: DEBUG! nlotabs @nlotabs\n";}

#
# Find NNLO tables
#
my @nnlotabs;
if ( -d "$nnlodir" ) {
    chdir $nnlodir;
    @nnlotabs = glob $nnloglob;
    chdir $sdir;
}
unless ( @nnlotabs ) {
    print "fastadd.pl: No NNLO table found in $nnlodir, now looking in $sdir ...\n";
    $nnlodir  = ".";
    @nnlotabs = glob $nnloglob;
} 
unless ( @nnlotabs ) {
    print "fastadd.pl: WARNING! No NNLO table found!\n";
}
$ntab = scalar @nnlotabs;
print "fastadd.pl: $ntab NNLO tables found.\n";
if ( $debug ) {print "fastadd.pl: DEBUG! nnlotabs @nnlotabs\n";}

#
# nlofast-add resp. fnlo-merge
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: ${merger}: $date\n";
# Statistics mode: Each LO with first NLO table & all LO with each NLO table
if ( $opt_s ) {
    print "\nfastadd.pl: Statistics mode\n";
    my $scmd2 = $cmd;
    foreach my $lotab ( @lotabs ) {
	my $scmd1 = $cmd;
	$scmd1 .= " $lodir/$lotab $nlodir/$nlotabs[0]";
	$scmd2 .= " $lodir/$lotab";
	$scmd1 .= " $lotab";
	$scmd1 =~ s/\.${tabext}$/\.tab/;
	print "fastadd.pl: Creating sum table for $lotab ...\n";
	if ( $debug ) {print "fastadd.pl: Running command $scmd1\n";}
	system("$scmd1 >> ${scen}_addst.log");
    }
    foreach my $nlotab ( @nlotabs ) {
	my $scmd1 = $scmd2;
	$scmd1 .= " $nlodir/$nlotab";
	$scmd1 .= " $nlotab";
	$scmd1 =~ s/\.${tabext}$/\.tab/;
	print "fastadd.pl: Creating sum table for $nlotab ...\n";
	if ( $debug ) {print "fastadd.pl: Running command $scmd1\n";}
	system("$scmd1 >> ${scen}_addst.log");
    }
# Normal mode: All LO with all NLO/NNLO tables
} else {
    my $scmd1 = $cmd;
    foreach my $lotab ( @lotabs ) {
	$scmd1 .= " $lodir/$lotab";
    }
    foreach my $nlotab ( @nlotabs ) {
	$scmd1 .= " $nlodir/$nlotab";
    }
    foreach my $nnlotab ( @nnlotabs ) {
	$scmd1 .= " $nnlodir/$nnlotab";
    }
    $scmd1 .= " $scen.tab";
    print "fastadd.pl: Creating total sum table for $scen ...\n";
    if ( $debug ) {print "fastadd.pl: Running command $scmd1\n";}
    system("$scmd1 > ${scen}_add.log");
}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: fastNLO finished: $date\n";
exit 0;
