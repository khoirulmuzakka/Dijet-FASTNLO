#!/usr/bin/env perl
#
# Unit test to read a normalisable fastNLO table
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
my $tab   = "InclusiveNJetEvents_fnl5662i_v23_fix_norm_1-12";
my $tabfl = ${tab}.".tab";
my $tabgz = ${tabfl}.".gz";
foreach my $file ( "cppnormtest.log", "cppnormdiff.log", "${tabfl}", "${tabgz}" ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
my $cmd = "cp -f ../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppnormtest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "gunzip ${tabgz}";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppnormtest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}

# Evaluate, keep only last 21 lines with numerical results
$cmd = "fnlo-tk-cppread ${tabfl} _ _ _ norm | tail -21 > cppnormtest.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppnormtest: Evaluating test table ${tabfl} failed: $ret, aborted!\n";}

# Determine difference to default evaluation output
$cmd = "diff ../data/check/${tab}_cppnorm.log cppnormtest.log > cppnormdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppnormtest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "cppnormdiff.log" ) {
    print "fnlo-tk-cppnormtest: Evaluation of test table differs from default:\n";
    $cmd = "cat cppnormdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    print "fnlo-tk-cppnormtest: Do you use LHAPDF version 6? Do you have the CT10nlo PDF set installed? Is it found?\n";
    die "fnlo-tk-cppnormtest: Table evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-cppnormtest: Table evaluation unit test passed.\n";

exit 0;
