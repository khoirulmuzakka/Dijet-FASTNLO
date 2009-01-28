#!/usr/bin/env perl 
#
# fastNLO table addition script:
#    This perl script creates either a total sum table (scen.tab) or
#    sum tables for statistical uncertainty evaluation: 
#    - For each single LO table with *one* 'dummy' NLO table
#    - For each single NLO table with *all* LO tables
#
# Version:
# 
# created by K. Rabbertz: 05.03.2006
# last modified:
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
our ( $opt_h, $opt_l, $opt_n, $opt_s, $opt_v ) = ( "", "", "", "", "" );
getopts('hl:n:sv') or die "fastadd.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastadd.pl\n";
    print "Usage: fastadd.pl [switches/options] scenario\n";
    print "  -h              Print this text\n";
    print "  -l dir          Directory for LO tables, (def.=scenario)\n";
    print "  -n dir          Directory for NLO/NNLO tables, (def.=scenario)\n";
    print "  -s              Produce tables for statistical evaluation,\n".
	"                  i.e. combinations of each LO with 1 NLO table and\n".
	"                  all LO with each NLO table\n";
    print "  -v              Verbose output\n\n";
    exit;
}

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fastadd.pl: Error! Need one scenario name!\n";
}
my $scen   = shift;

#
# Initialization
#
my $lodir    = "${scen}"; 
if ( $opt_l ) {$lodir = $opt_l;}
my $nlodir   = "${scen}"; 
if ( $opt_n ) {$nlodir = $opt_n;}
my $nnlodir  = $nlodir; 
my $loglob   = "${scen}*born*.raw*";
my $nloglob  = "${scen}*nlo*.raw*";
my $nnloglob = "${scen}*thrcor*.raw*";
# Directory
my $sdir = getcwd();

#
# Check on nlofast-add
#
my $cmd = `which nlofast-add`;
chomp $cmd;
unless ( $cmd ) {
    if ( $ENV{NLOJET} ) {
	if ( -f "$ENV{NLOJET}/bin/nlofast-add" ) {
	    $cmd = "$ENV{NLOJET}/bin/nlofast-add";
	} else {
	    die "fastadd.pl: ERROR! nlofast-add command not found, ".
		"neither via \`which\` nor in $ENV{NLOJET}/bin, aborted!\n";
	}
    } else {
	die "fastadd.pl: ERROR! nlofast-add command not found ".
	    "via \`which\` and \$NLOJET is not set, aborted!\n";
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
if ( $opt_v ) {print "fastadd.pl: DEBUG! lotabs @lotabs\n";}

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
if ( $opt_v ) {print "fastadd.pl: DEBUG! nlotabs @nlotabs\n";}

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
if ( $opt_v ) {print "fastadd.pl: DEBUG! nnlotabs @nnlotabs\n";}

#
# nlofast-add
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: nlofast-add: $date\n";
# Statistics mode: Each LO with first NLO table & all LO with each NLO table
if ( $opt_s ) {
    print "\nfastadd.pl: Statistics mode\n";
    my $scmd2 = $cmd;
    foreach my $lotab ( @lotabs ) {
	my $scmd1 = $cmd;
	$scmd1 .= " $lodir/$lotab $nlodir/$nlotabs[0]";
	$scmd2 .= " $lodir/$lotab";
	$scmd1 .= " $lotab";
	$scmd1 =~ s/\.raw$/\.tab/;
	print "fastadd.pl: Creating sum table for $lotab ...\n";
	if ( $opt_v ) {print "fastadd.pl: Running command $scmd1\n";}
	system("$scmd1 >> ${scen}_addst.log");
    }
    foreach my $nlotab ( @nlotabs ) {
	my $scmd1 = $scmd2;
	$scmd1 .= " $nlodir/$nlotab";
	$scmd1 .= " $nlotab";
	$scmd1 =~ s/\.raw$/\.tab/;
	print "fastadd.pl: Creating sum table for $nlotab ...\n";
	if ( $opt_v ) {print "fastadd.pl: Running command $scmd1\n";}
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
    if ( $opt_v ) {print "fastadd.pl: Running command $scmd1\n";}
    system("$scmd1 > ${scen}_add.log");
}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: fastNLO finished: $date\n";
exit 0;
