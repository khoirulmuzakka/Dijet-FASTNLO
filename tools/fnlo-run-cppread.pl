#!/usr/bin/env perl
#
# Run fnlo-tk-cppread on all matching fastNLO tables in cwd and store log file
# Version:
#
# created by K. Rabbertz: 10.10.2006
# adapted by K. Rabbertz from fnlo-check-infnan.pl: 02.07.2017
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
if ( $? ) {die "fnlo-run-cppread.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################################\n";
print "# fnlo-run-cppread.pl: Starting cppread evaluation of all matching tables: FASTCPPREADEVAL_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_v ) = ( "", "2.3" );
getopts('hv:') or die "fnlo-run-cppread.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-run-cppread.pl\n";
    print "Usage: fnlo-run-cppread.pl glob (selects all files matching glob)\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-run-cppread.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-run-cppread.pl: Error! Need one file name glob to select files to check!\n";
}
my $glstr = shift;
chomp $glstr;
print "fnlo-run-cppread.pl: Checking for file glob $glstr ...\n";
my @files = glob "*${glstr}*.tab*";
chomp @files;

# tmp
my $pdf    = "CT14nnlo";
my $nvars  = "_";
my $ascode = "LHAPDF";
my $norm   = "_";
my $flex   = "scale2";

#
# fnlo-tk-cppread
#
foreach my $file (@files) {
    print "fnlo-run-cppread.pl: Evaluating table $file\n";
    my $logfil = $file;
    $logfil =~ s/tab\.gz$/log/;
    $logfil =~ s/tab$/log/;
    my $cmd = "fnlo-tk-cppread $file $pdf $nvars $ascode $norm $flex > $logfil";
    system($cmd);
    my $ret = $?;
    if ( ! $ret ) {
        print "fnlo-run-cppread.pl: WARNING! Problem in evaluation of table $file\n";
    }
}

exit 0;
