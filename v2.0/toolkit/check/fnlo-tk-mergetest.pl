#!/usr/bin/env perl
#
# Unit test to catenate two fastNLO tables
# Version:
#
# created by K. Rabbertz: 20.10.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

# Remove potentially left-over temporary files
my @tabs = ("InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_0000", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_0001", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_0002", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_0100", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_0101", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_0102", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix" );
my @tabfls;
my @tabgzs;
foreach my $tab ( @tabs ) {
    push @tabfls, ${tab}.".tab";
    push @tabgzs, ${tab}.".tab".".gz";
}
foreach my $file ( "mergetest.tab", "mergediff.log", @tabfls, @tabgzs ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
foreach my $tabgz ( @tabgzs ) {
    my $cmd = "cp -f ../data/check/${tabgz} .";
    print "Executing command: $cmd\n";
    my $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-mergetest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
    $cmd = "gunzip ${tabgz}";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-mergetest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}
}

# Merge
my $cmd = "fnlo-tk-merge $tabfls[0] $tabfls[1] $tabfls[2] $tabfls[3] $tabfls[4] $tabfls[5] mergetest.tab";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-mergetest: Merging test tables $tabfls[0] to $tabfls[5] failed: $ret, aborted!\n";}

# Determine difference to default catenated table
$cmd = "diff $tabfls[6] mergetest.tab > mergediff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-mergetest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "mergediff.log" ) {
    print "fnlo-tk-mergetest: Merged test table differs from default:\n";
    $cmd = "cat mergediff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-mergetest: Catenation unit test failed, please fix!\n";
}

print "fnlo-tk-mergetest: Catenation unit test passed.\n";

exit 0;
