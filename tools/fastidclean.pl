#!/usr/bin/env perl 
#
# Move identical fastNLO tables out of the way into subdirectory
# Version:
# 
# K. Rabbertz: 10.10.2006
# last modified:
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
#use FindBin qw($Bin);
#use lib "$Bin/modules";
use Getopt::Std;
use strict;
use warnings;

my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n######################################################\n";
print "# fastidclean.pl: Starting table clean-up for fastNLO: FASTIDCLEAN_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h, $opt_v ) = ( ".", "", "1" );
getopts('d:hv:') or die "fastidclean.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastidclean.pl\n";
    print "Usage: fastidclean.pl scenario\n";
    print "  -d dir          Directory for scenario and scenarioref subdirs, (def.=.)\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO version 1 or 2 (def.=1)\n\n";
    exit;
}

#
# Initialization
#
my $scenario = shift;
if ( ! $scenario ) {
    die "fastidclean.pl: ERROR! No scenario specified\n";
}
chomp $scenario;

my @subdirs = ($scenario, $scenario."ref");
my @globs = ("\\\*born\\\*.raw", "\\\*nlo\\\*.raw");
my $vers  = $opt_v;
if ($vers == 2) {
    push(@subdirs,$scenario."wrm");
    @globs = ("\\\*born\\\*.tab", "\\\*nlo\\\*.tab", "\\\*warmup\\\*.dat");
}

#
# Clean
#
chdir $opt_d or die "fastidclean.pl: ERROR! Could not cd to directory $opt_d, aborted.\n";
my $cwd = `pwd`;
chomp $cwd;
foreach my $subdir (@subdirs) {
    chdir $subdir or die "fastidclean.pl: ERROR! Could not cd to subdirectory $subdir, aborted.\n";;
    foreach my $glob (@globs) {
	print "\nfastidclean.pl: Calling fastidcheck.pl for filename glob $glob in directory $subdir ...\n";
	system("fastidcheck.pl -v $vers $glob");
    }
    chdir $cwd;
}

exit 0;
