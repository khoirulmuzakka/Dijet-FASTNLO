#!/usr/bin/env perl
#
# CIJET xsection script
# Version:
#
# created by K. Rabbertz: 11.07.2016
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
print "# cijet-xsec.pl: Start running CIJETXSEC: CIJETXSEC_$date\n";
print "##################################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "cijet-xsec.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\ncijet-xsec.pl\n";
    print "Usage: cijet-xsec.pl [switches/options] convfile label\n";
    print "  -h              Print this text\n\n";
    exit;
}

#
# Parse arguments
#
unless ( @ARGV == 2 ) {
    die "cijet-xsec.pl: Error! Need two arguments!\n";
}
my $convfile = shift;
my $label    = shift;

#
# Print some info at beginning of job
#
print "cijet-xsec.pl: Print initial listing of cwd:\n";
my $cmd = "ls -la";
my $ret = system("$cmd");
print "\n";

#
# Unpack
#
my $outfile = $convfile;
$outfile    =~ s/.conv$//;
$outfile   .= "_${label}.xsc";

#
# Run cijet xsection
#
print "cijet-xsec.pl: Running CIJET xsection in current directory $ENV{PWD}.\n";
$cmd = "./cixsec $convfile $outfile";
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "cijet-xsec.pl: Starting xsection: CIJETXSEC0_$date\n";
print "cijet-xsec.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "cijet-xsec.pl: ERROR! Error $ret in CIJET xsection step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\ncijet-xsec.pl: Xsection finished: CIJETXSEC1_$date\n";

#
# Print some info at end of job
#
print "\ncijet-xsec.pl: Print final listing of cwd:\n";
$cmd = "ls -la";
$ret = system("$cmd");
print "\n";

exit 0;
