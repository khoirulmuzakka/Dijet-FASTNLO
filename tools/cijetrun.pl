#!/usr/bin/env perl
#
# CIJET run script
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
print "# cijetrun.pl: Start running CIJET: CIJETRUN_$date\n";
print "##################################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "cijetprep.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\ncijetrun.pl\n";
    print "Usage: cijetprep.pl [switches/options]\n";
    print "  -h              Print this text\n\n";
    exit;
}
my $arcname = "cijet-bin";

#
# Unpack archive in cwd
#
print "cijetrun.pl: Unpacking CIJET archive in current directory $ENV{PWD}/..\n";
my $cmd = "tar xzvf ${arcname}.tgz";
my $ret = system("$cmd");
if ( $ret ) {die "cijetrun.pl: ERROR! Could not unpack archive ${arcname}.tgz: $ret\n";}
#
# Run cijet
#
$ENV{LHAPATH} = $ENV{PWD};
print "cijetrun.pl: Running CIJET in current directory $ENV{PWD} with LHAPATH set to $ENV{LHAPATH}.\n";
$cmd = "./dijets4ci_sin CIJET_MassBin-@MASSLOW@-@MASSUPP@_ChiBin-@CHILOW@-@CHIUPP@";
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "fastrun.pl: Starting calculation: CIJET0_$date\n";
print "fastrun.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "cijetrun.pl: ERROR! Error $ret in CIJET run step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\ncijetrun.pl: Calculation finished: CIJET1_$date\n";

#$cmd = "mv output.dat @CONF@_Lambda-@LAMBDA@_MassBin-@MASSLOW@-@MASSUPP@_@MY_JOBCNT@.dat";
$ret = system("$cmd");
if ( $ret ) {die "cijetrun.pl: ERROR! Error $ret in CIJET output move step, aborted!\n";}

exit 0;
