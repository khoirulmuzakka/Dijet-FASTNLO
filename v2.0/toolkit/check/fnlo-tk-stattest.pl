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
my $tabl = "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat";
my $tabn = "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat";
my @tabs = ("InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0000", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0001", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0002", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0100", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0101", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0102" );
my @tabfls;
my @tabgzs;
foreach my $tab ( @tabs ) {
    push @tabfls, ${tab}.".tab";
    push @tabgzs, ${tab}.".tab".".gz";
}
foreach my $file ( "statlotest.log", "statnlotest.log", "statlodiff.log", "statnlodiff.log", @tabfls, @tabgzs ) {
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
    if ( $ret ) {die "fnlo-tk-stattest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
    $cmd = "gunzip ${tabgz}";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-stattest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}
}

# Statistical evaluation
my $cmd = "fnlo-tk-statunc $tabl | tail -42 > statlotest.log";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Statistical evaluation of LO test tables ${tabl} failed: $ret, aborted!\n";}
my $cmd = "fnlo-tk-statunc $tabn | tail -42 > statnlotest.log";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Statistical evaluation of NLO test tables ${tabn} failed: $ret, aborted!\n";}

# Determine difference to default statistical uncertainties
$cmd = "diff ../data/check/${tabl}.log statlotest.log > statlodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Result comparison with LO 'diff' failed: $ret, aborted!\n";}
$cmd = "diff ../data/check/${tabn}.log statnlotest.log > statnlodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Result comparison with NLO 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "statlodiff.log" ) {
    print "fnlo-tk-stattest: Statistical evaluation of LO test tables differs from default:\n";
    $cmd = "cat statlodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-stattest: Statistical evaluation unit test failed, please fix!\n";
}
if ( ! -z "statnlodiff.log" ) {
    print "fnlo-tk-stattest: Statistical evaluation of NLO test tables differs from default:\n";
    $cmd = "cat statnlodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-stattest: Statistical evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-stattest: Statistical evaluation unit test passed.\n";

exit 0;
