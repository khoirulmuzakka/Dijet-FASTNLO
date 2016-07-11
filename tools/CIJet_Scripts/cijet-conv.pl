#!/usr/bin/env perl
#
# CIJET conversion script
# Version:
#
# created by K. Rabbertz: 09.07.2016
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
print "# cijet-conv.pl: Start running CIJETCONV: CIJETCONV_$date\n";
print "##################################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "cijet-conv.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\ncijet-conv.pl\n";
    print "Usage: cijet-conv.pl [switches/options] PDFset PDFmember gridfile pdflabel\n";
    print "  -h              Print this text\n\n";
    exit;
}

#
# Parse arguments
#
unless ( @ARGV == 4 ) {
    die "cijet-conv.pl: Error! Need four arguments!\n";
}
my $pdfset   = shift;
my $pdfmem   = shift;
my $gridfile = shift;
my $gridgz   = 0;
if ( $gridfile =~ m/.gz$/ ) {$gridgz = 1};
my $pdflab   = shift;

#
# Print some info at beginning of job
#
print "cijet-conv.pl: Print initial listing of cwd:\n";
my $cmd = "ls -la";
my $ret = system("$cmd");
print "\n";

#
# Unpack
#
if ( $gridgz == 1 ) {
    $cmd = "gunzip ${gridfile}";
    $ret = system("$cmd");
    $gridfile =~ s/.gz$//;
}
$cmd = "mkdir fgrid";
$ret = system("$cmd");
$cmd = "mv $gridfile fgrid";
$ret = system("$cmd");
my $outfile = $gridfile;
$outfile =~ s/.sum$//;
$outfile .= "_${pdflab}_${pdfmem}.conv";

#
# Run cijet conversion
#
#
# LHAPDF needs to be found together with requested PDF sets!
# $ENV{LHAPATH} = $ENV{PWD};
# print "cijet-conv.pl: Running CIJET conversion in current directory $ENV{PWD} with LHAPATH set to $ENV{LHAPATH}.\n";
#
print "cijet-conv.pl: Running CIJET conversion in current directory $ENV{PWD}.\n";
$cmd = "./ciconv $pdfset $pdfmem $gridfile $outfile";
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "cijet-conv.pl: Starting conversion: CIJETCONV0_$date\n";
print "cijet-conv.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "cijet-conv.pl: ERROR! Error $ret in CIJET conversion step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\ncijet-conv.pl: Conversion finished: CIJETCONV1_$date\n";

#
# Print some info at end of job
#
print "\ncijet-conv.pl: Print final listing of cwd:\n";
$cmd = "ls -la";
$ret = system("$cmd");
print "\n";

exit 0;
