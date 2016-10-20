#!/usr/bin/env perl
#
# Unit test to read a fastNLO table
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

my $tab   = "fnl1014_I902309_2-2";
my $tabfl = ${tab}.".tab";
my $tabgz = ${tabfl}.".gz";
if ( -e "${tabfl}" ) {system("rm -f ${tabfl}");}
if ( -e "${tabgz}" ) {system("rm -f ${tabgz}");}
my $cmd = "cp -f ../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppunittest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "gunzip ${tabgz}";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppunittest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "fnlo-tk-cppread ${tabfl} > ${tab}_cpptest.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppunittest: Cppreading test table ${tabfl} failed: $ret, aborted!\n";}
$cmd = "diff ../data/check/${tab}_cpp.log ${tab}_cpptest.log > ${tab}_cppdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppunittest: Diffing test table result ${tab}_cpptest.log failed: $ret, aborted!\n";}
if ( ! -z "${tab}_cppdiff.log" ) {
    print "fnlo-tk-cppunittest: Evaluation of test table ${tab} differs from default:\n";
    $cmd = "cat ${tab}_cppdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-cppunittest: cpp unit test failed, please fix!\n";
}
system("rm -f ${tabfl}");
system("rm -f ${tab}_cpptest.log");
system("rm -f ${tab}_cppdiff.log");
print "fnlo-tk-cppunittest: cpp unit test passed.\n";

exit 0;
