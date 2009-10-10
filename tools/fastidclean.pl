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
our ( $opt_d, $opt_h ) = ( ".", "" );
getopts('d:h') or die "fastidclean.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastidclean.pl\n";
    print "Usage: fastidclean.pl scenario\n";
    print "  -d dir          Directory for scenario and scenarioref subdirs, (def.=.)\n";
    print "  -h              Print this text\n\n";
    exit;
}

#
# Initialization
#
my $scenario = shift;
chomp $scenario;
my @subdirs = ($scenario, $scenario."ref");
my @globs = ("\\\*born\\\*.raw", "\\\*nlo\\\*.raw");

#
#
#
chdir $opt_d or die "fastidclean.pl: ERROR! Could not cd to directory $opt_d, aborted.\n";
my $cwd = `pwd`;
chomp $cwd;
foreach my $subdir (@subdirs) {
    chdir $subdir or die "fastidclean.pl: ERROR! Could not cd to subdirectory $subdir, aborted.\n";;
    foreach my $glob (@globs) {
	print "\nfastidclean.pl: Calling fastidcheck.pl for filename glob $glob in directory $subdir ...\n";
	system("fastidcheck.pl $glob");
    }
    chdir $cwd;
}

exit 0;
