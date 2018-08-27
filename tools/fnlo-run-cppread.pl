#!/usr/bin/env perl
#
# Run fnlo-tk-cppread on all matching fastNLO tables in cwd and
# store log file for numerical crosschecks
#
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
our ( $opt_a, $opt_h, $opt_n, $opt_p, $opt_s ) = ( "LHAPDF", "", 1, "CT10nlo", "scale1" );
getopts('a:hn:p:s:') or die "fnlo-run-cppread.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-run-cppread.pl\n";
    print "Usage: fnlo-run-cppread.pl glob (selects all files matching glob)\n";
    print "  -a ascode       Choose alpha_s evolution code (def.=LHAPDF)\n";
    print "                  For alternatives check output of tnlo-tk-cppread -h\n";
    print "  -h              Print this text\n";
    print "  -n nvars        Choose scale variations (def.=1)\n";
    print "  -p pdfset       Choose PDF set (def.=CT10nlo)\n";
    print "  -s scale        Choose mur & muf scales (def.=scale1)\n";
    print "                  For comparison of flex-scale tables to NNLOJET use:\n";
    print "                  - DIS: scale12, DIS inclusive jets with pTjet: scale1\n";
    print "                  - pp : scale21, pp  inclusive jets with pTjet: scale2\n\n";
    exit;
}
my $ascode = $opt_a;
my $nvars  = $opt_n;
my $pdf    = $opt_p;
my $norm   = "_";
my $scale  = $opt_s;



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


#
# fnlo-tk-cppread
#
foreach my $file (@files) {
    print "fnlo-run-cppread.pl: Evaluating table $file\n";
    my $logfil = $file;
    $logfil =~ s/\.tab\.gz$//;
    $logfil =~ s/\.tab$//;
    my $log = $logfil.".log";
    if ( ! -f $logfil ) {
        my $cmd = "fnlo-tk-cppread $file $pdf $nvars $ascode $norm $scale &> $log";
        system($cmd);
        my $ret = $?;
        if ( ! $ret ) {
            #        print "fnlo-run-cppread.pl: WARNING! Problem in evaluation of table $file\n";
        }
    } else {
        print "fnlo-run-cppread.pl: WARNING! log file $logfil exists already, skipped!\n";
    }
}

exit 0;
