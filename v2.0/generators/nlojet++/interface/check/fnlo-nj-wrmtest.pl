#!/usr/bin/env perl
#
# Unit test to produce a warmup file for a NLOjet++ scenario
# Version:
#
# created by K. Rabbertz: 23.10.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

# Remove potentially left-over temporary files
my $tabb  = "fnr0001midpHT_I723509_v23_fix";
my $steer = ${tabb}.".str";
my $tabw  = ${tabb}."_InclusiveNJets_warmup.txt";

foreach my $file ( $steer, $tabw, "InclusiveNJets.str", "wrmdiff.log" ) {
    if ( -e $file || -l $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy defaults from data/check storage
my $cmd = "cp -f ../data/check/${steer} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-wrmtest: Copying test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "ln -s ${steer} InclusiveNJets.str";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-wrmtest: Linking test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "fnlo-tk-config --libdir";
print "Executing command: $cmd\n";
my $libdir = `fnlo-tk-config --libdir`;
if ( $ret ) {die "fnlo-nj-wrmtest: Determining lib dir failed: $ret, aborted!\n";}
chomp $libdir;

# Produce warmup file
$cmd = "nlojet++ --calculate -cnlo -n fnr0001midpHT_I723509_v23_fix -u ${libdir}/fastnlo_interface_nlojet/libInclusiveNJets.la --max-event=10000 --save-after=1000 -s 12345";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-wrmtest: Producing warmup table for test scenario failed: $ret, aborted!\n";}

# Determine difference to default warmup table
$cmd = "diff ../data/check/${tabw} ${tabw} > wrmdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-wrmtest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "wrmdiff.log" ) {
    print "fnlo-nj-wrmtest: Warmup test table differs from default:\n";
    $cmd = "cat wrmdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-nj-wrmtest: Warmup table production unit test failed, please fix!\n";
}

print "fnlo-nj-wrmtest: Warmup table production unit test passed.\n";

exit 0;
