#!/usr/bin/env perl
#
# Unit test to produce a LO fastNLO table for a NLOjet++ scenario
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
my $tabl  = ${tabb}."-hhc-born-2jet.tab";
my $tabgz = ${tabl}.".gz";
my $tabw  = ${tabb}."_InclusiveNJets_warmup.txt";

foreach my $file ( $steer, $tabw, $tabl, $tabgz, "InclusiveNJets.str", "lodiff.log" ) {
    if ( -e $file || -l $file ) {
        system("rm -f $file");
    }
}
if ( -e "output/${tabl}" ) {system("rm -f output/${tabl}");}

# Prepare test setup
# Copy defaults from data/check storage
my $cmd = "cp -f ../data/check/${steer} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-lotest: Copying test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "ln -s ${steer} InclusiveNJets.str";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-lotest: Linking test steering ${steer} failed: $ret, aborted!\n";}
$cmd = "fnlo-tk-config --libdir";
print "Executing command: $cmd\n";
my $libdir = `fnlo-tk-config --libdir`;
if ( $ret ) {die "fnlo-nj-lotest: Determining lib dir failed: $ret, aborted!\n";}
chomp $libdir;
$cmd = "cp -f ../data/check/${tabw} .";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-lotest: Copying test warmup ${tabw} failed: $ret, aborted!\n";}
$cmd = "cp -f ../data/check/${tabgz} .";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-lotest: Copying LO table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "gunzip ${tabgz}";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-lotest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}

# Produce LO table
$cmd = "nlojet++ --calculate -cborn -n fnr0001midpHT_I723509_v23_fix -u ${libdir}/fastnlo_interface_nlojet/libInclusiveNJets.la --max-event=10000 --save-after=1000 -s 12345";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-lotest: Producing LO table for test scenario failed: $ret, aborted!\n";}

# Determine difference to default LO table
$cmd = "diff ./${tabl} output/${tabl} > lodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-nj-lotest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "lodiff.log" ) {
    print "fnlo-nj-lotest: LO test table differs from default:\n";
    $cmd = "cat lodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-nj-lotest: LO table production unit test failed, please fix!\n";
}

print "fnlo-nj-lotest: LO table production unit test passed.\n";

exit 0;
