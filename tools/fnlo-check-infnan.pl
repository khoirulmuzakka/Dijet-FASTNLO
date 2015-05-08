#!/usr/bin/env perl
#
# Check on fastNLO tables contaminated by 'inf' or 'NaN' entries
# Version:
#
# created by K. Rabbertz: 10.10.2006
# adapted by K. Rabbertz from fastidclean.pl: 23.04.2015
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-check-infnan.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################################\n";
print "# fnlo-check-infnan.pl: Starting table contamination check for fastNLO: FASTINFNANCHECK_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_v ) = ( "", "2.3" );
getopts('hv:') or die "fnlo-check-infnan.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-check-infnan.pl\n";
    print "Usage: fnlo-check-infnan.pl glob (selects all files matching glob)\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-check-infnan.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-check-infnan.pl: Error! Need one file name glob to select files to check!\n";
}
my $glstr = shift;
chomp $glstr;
print "fnlo-check-infnan.pl: Checking for file glob $glstr ...\n";
my @files = glob "*${glstr}*";
chomp @files;
if ( ! -d "Problems" ) {
    print "fnlo-check-infnan.pl: Creating subdirectory Problems ...\n";
    system("mkdir Problems");
}

#
# Grep
#
foreach my $file (@files) {
    print "fnlo-check-infnan.pl: Checking file $file\n";
    my $hasinfnan = 0;
    my $cmd = "grep -i -l \"nan\" $file";
    system($cmd);
    my $ret = $? >> 8;
    if ( ! $ret ) {
        print "fnlo-check-infnan.pl: WARNING! Found nan in file $file\n";
        $hasinfnan = 1;
    } else {
        my $cmd = "grep -i -l \"inf\" $file";
        system($cmd);
        my $ret = $? >> 8;
        if ( ! $ret ) {
            print "fnlo-check-infnan.pl: WARNING! Found inf in file $file\n";
            $hasinfnan = 1;
        }
    }
    if ( $hasinfnan ) {
        my $logfil = $file;
        $logfil =~ s/\.tab/\.log/;
        my $errfil = $file;
        $errfil =~ s/\.tab/\.err/;
        my $ret = system("mv $file Problems");
        if ( $ret ) {die "fnlo-check-infnan.pl: Couldn't move file $file into ".
                         "Problems: $ret, aborted!\n";}
        if (-f $logfil) {
            my $ret = system("mv $logfil Problems");
            if ( $ret ) {die "fnlo-check-infnan.pl: Couldn't move file $logfil into ".
                             "Problems: $ret, aborted!\n";}
        }
        if (-f $errfil) {
            my $ret = system("mv $errfil Problems");
            if ( $ret ) {die "fnlo-check-infnan.pl: Couldn't move file $errfil into ".
                             "Problems: $ret, aborted!\n";}
        }
    }
}

exit 0;
