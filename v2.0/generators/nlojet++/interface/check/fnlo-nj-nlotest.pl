#!/usr/bin/env perl
#
# Unit test to produce a NLO fastNLO table for a NLOjet++ scenario
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
my $tabn  = ${tabb}."-hhc-nlo-2jet.tab";
my $tabgz = ${tabn}.".gz";
my $tabw  = ${tabb}."_InclusiveNJets_warmup.txt";

foreach my $file ( $steer, $tabw, $tabn, $tabgz, "InclusiveNJets.str", "nlodiff.log" ) {
    if ( -e $file || -l $file ) {
        system("rm -f $file");
    }
}
if ( -e "output/${tabn}" ) {system("rm -f output/${tabn}");}

# Prepare test setup
# Copy defaults from data/check storage
my $cmd = "cp -f ../data/check/${steer} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-nlotest: Copying test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "ln -s ${steer} InclusiveNJets.str";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-nlotest: Linking test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "fnlo-tk-config --libdir";
print "Executing command: $cmd\n";
my $libdir = `fnlo-tk-config --libdir`;
if ( $ret ) {die "fnlo-nj-nlotest: Determining lib dir failed: $ret, aborted!\n";}
chomp $libdir;
my $cmd = "cp -f ../data/check/${tabw} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-nlotest: Copying test warmup ${tabw} failed: $ret, aborted!\n";}
my $cmd = "cp -f ../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-nlotest: Copying NLO table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "gunzip ${tabgz}";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-nlotest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}

# Produce NLO table
$cmd = "nlojet++ --calculate -cnlo -n fnr0001midpHT_I723509_v23_fix -u ${libdir}/fastnlo_interface_nlojet/libInclusiveNJets.la --max-event=10000 --save-after=1000 -s 12345";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-nlotest: Producing NLO table for test scenario failed: $ret, aborted!\n";}

# Determine difference to default NLO table
$cmd = "diff ./${tabn} output/${tabn} > nlodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-nlotest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "nlodiff.log" ) {
    print "fnlo-nj-nlotest: NLO test table differs from default:\n";
    $cmd = "cat nlodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-nj-nlotest: NLO table production unit test failed, please fix!\n";
}

print "fnlo-nj-nlotest: NLO table production unit test passed.\n";

exit 0;
