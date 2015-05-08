#!/usr/bin/env perl
#
# CIDIJET run script
# Version:
#
# created by K. Rabbertz: 08.03.2011
# last modified:
#
#-----------------------------------------------------------------------
# Todo:
#
use Cwd;
use English;
#use FindBin qw($Bin);
#use lib "$Bin/modules";
use Getopt::Std;
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n##################################################################\n";
print "# cidijetrun.pl: Start running CIDIJET: CIDIJETRUN_$date\n";
print "##################################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "cidijetprep.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\ncidijetrun.pl\n";
    print "Usage: cidijetprep.pl [switches/options]\n";
    print "  -h              Print this text\n\n";
    exit;
}
my $arcname = "cidijet-bin";

#
# Unpack archive in cwd
#
print "cidijetrun.pl: Unpacking CIDIJET archive in current directory $ENV{PWD}/..\n";
my $cmd = "tar xzvf ${arcname}.tgz";
my $ret = system("$cmd");
if ( $ret ) {die "cidijetrun.pl: ERROR! Could not unpack archive ${arcname}.tgz: $ret\n";}
#
# Run cidijet
#
print "cidijetrun.pl: Running CIDIJET in current directory $ENV{PWD}/..\n";
$cmd = "./cidijet";
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "fastrun.pl: Starting calculation: CIDIJET0_$date\n";
print "fastrun.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "cidijetrun.pl: ERROR! Error $ret in CIDIJET run step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\ncidijetrun.pl: Calculation finished: CIDIJET1_$date\n";

$cmd = "mv output.dat @CONF@_Lambda-@LAMBDA@_MassBin-@MASSLOW@-@MASSUPP@_Order-@ORDER@_xmu-@XMU@_@MY_JOBCNT@.dat";
$ret = system("$cmd");
if ( $ret ) {die "cidijetrun.pl: ERROR! Error $ret in CIDIJET output move step, aborted!\n";}

exit 0;
