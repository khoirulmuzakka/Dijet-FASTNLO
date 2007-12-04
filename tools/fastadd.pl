#!/usr/bin/env perl 
#
# fastNLO table addition script
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
print "# fastadd.pl: Starting table addition for fastNLO: $date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_l, $opt_n ) = ( "", "", "" );
getopts('hl:n:') or die "fastadd.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastadd.pl\n";
    print "Usage: fastadd.pl [switches/options] scenario\n";
    print "  -h              Print this text\n";
    print "  -l dir          Directory for LO tables, (def.=scenario_LO_tables)\n";
    print "  -n dir          Directory for NLO tables, (def.=scenario_NLO_tables)\n\n";
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
my $lodir   = "${scen}_LO_tables"; 
if ( $opt_l ) {$lodir = $opt_l;}
my $nlodir  = "${scen}_NLO_tables"; 
if ( $opt_n ) {$nlodir = $opt_n;}
my $loglob  = "${scen}*born*.raw*";
my $nloglob = "${scen}*nlo*.raw*";
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
print "fastadd.pl: DEBUG! lotabs @lotabs\n";

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
print "fastadd.pl: DEBUG! nlotabs @nlotabs\n";

#
# nlofast-add
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: nlofast-add: $date\n";
foreach my $lotab ( @lotabs ) {
    $cmd .= " $lodir/$lotab";
}
foreach my $nlotab ( @nlotabs ) {
    $cmd .= " $nlodir/$nlotab";
}
$cmd .= " $scen.tab";
print "fastrun.pl: Running command $cmd\n";
system("$cmd >& $scen.log");
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: fastNLO finished: $date\n";
exit 0;
