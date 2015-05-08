#!/usr/bin/env perl
#
# Clean up fastNLO table files with respect to duplicates or inf or NaN entries
# Version:
#
# created by K. Rabbertz: 10.10.2006
# adapted by K. Rabbertz from fastidclean.pl: 23.04.2015
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-clean-tables.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################################\n";
print "# fnlo-clean-tables.pl: Starting table clean-up for fastNLO: FASTIDCLEAN_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h, $opt_i, $opt_n, $opt_v ) = ( ".", "", "", "", "2.3" );
getopts('d:hinv:') or die "fnlo-clean-tables.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-clean-tables.pl\n";
    print "Usage: fnlo-clean-tables.pl scenario\n";
    print "  -d dir          Directory for scenario and scenarioref subdirs, (def.=.)\n";
    print "  -h              Print this text\n";
    print "  -i              Run check on duplicates (id) (def.=no)\n";
    print "  -n              Run check on inf or NaN content (def.=no)\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
my $dir = "";
if ( $opt_d ne "." ) {$dir = $opt_d;}
my $id  = $opt_i;
my $nan = $opt_n;
unless ( $id || $nan ) {
    print "fnlo-clean-tables.pl: Nothing to be done, exiting.\n";
    print "fnlo-clean-tables.pl: Specify at least one check option -i or -n\n";
    exit(0);
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-clean-tables.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-clean-tables.pl: Error! No scenario specified!\n";
}
my $scenario = shift;
chomp $scenario;
my @subdirs = ($scenario, $scenario."ref");
my @globs = ("\\\*-born\\\*.tab", "\\\*-nlo\\\*.tab", "\\\*-nnlo\\\*.tab");

#
# Clean
#
if ($dir) {chdir $dir or die "fnlo-clean-tables.pl: ERROR! Could not cd to directory $dir, aborted.\n";}
my $cwd = `pwd`;
chomp $cwd;
foreach my $subdir (@subdirs) {
    unless (chdir $subdir) {
        print "fnlo-clean-tables.pl: Warning! Could not cd to subdirectory $subdir, skipped.\n";
        next;
    }
# Clean "nan" and "inf"
    foreach my $glob (@globs) {
        if ( $nan ) {
            print "\nfnlo-clean-tables.pl: Calling fnlo-check-infnan.pl for filename glob $glob in directory $subdir ...\n";
            my $ret = system("fnlo-check-infnan.pl -v $vers $glob");
            if ( $ret ) {die "fnlo-clean-tables.pl: fnlo-check-infnan.pl failed: $ret, aborted!\n";}
        }
    }
# Clean duplicates
    foreach my $glob (@globs) {
        if ( $id ) {
            print "\nfnlo-clean-tables.pl: Calling fnlo-check-duplicates.pl for filename glob $glob in directory $subdir ...\n";
            my $ret = system("fnlo-check-duplicates.pl -v $vers $glob");
            if ( $ret ) {die "fnlo-clean-tables.pl: fnlo-check-duplicates.pl failed: $ret, aborted!\n";}
        }
    }
    chdir $cwd;
}

exit 0;
