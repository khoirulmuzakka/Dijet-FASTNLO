#!/usr/bin/env perl 
#
# fastNLO grid submission archiver
# Version:
# 
# created by K. Rabbertz: 10.12.2007
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
print "\n#####################################################\n";
print "# fastprep.pl: Starting archive creation for fastNLO: $date\n";
print "#####################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_p ) = ( "", "" );
getopts('hp:') or die "fastprep.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastprep.pl\n";
#    print "Usage: fastprep.pl [switches/options] scenario\n";
    print "Usage: fastprep.pl [switches/options]\n";
    print "  -h              Print this text\n";
    print "  -p pdf          Add CTEQ parton densities or LHAPDF\n\n";
    exit;
}
unless ( $opt_p eq "" || $opt_p eq "CTEQ" || $opt_p eq "LHAPDF" ) {
    die "fastrun.pl: Error! Illegal option -p $opt_p, aborted.\n";
}
my $pdf   = $opt_p;
my $arcname = "fastNLO-bin.tgz";
if ( $pdf ) {
    $arcname = "fastNLO-bin-${pdf}.tgz";
}

#
# Starting archive creation
#
if ( ! $ENV{FASTNLO} ) {die "fastprep.pl: ERROR! Environment variable ".
			    "FASTNLO not set, aborted!\n";}
print "fastprep.pl: Preparing fastNLO archive for submission in directory $ENV{FASTNLO}/..\n";
print "fastprep.pl: Only hh collisions, no LHAPDF for now!\n";
chdir "$ENV{FASTNLO}/.." or die "fastprep.pl: ERROR! Could not cd to $ENV{FASTNLO}/..!\n";

my $cmd = "tar cfz $arcname lib/lib*.so* nlojet/bin nlojet/lib fastjet/lib fastjet/plugins/SISCone/.libs fastjet/plugins/SISCone/siscone/siscone/.libs fastNLO/trunk/common/* fastNLO/trunk/tools/* fastNLO/trunk/v1.4/author1c/hadron/*.la fastNLO/trunk/v1.4/author1c/hadron/.libs";
if ( -d "lib64" ) {
    $cmd .= " lib64/lib*.so*";
}
if ( $pdf eq "CTEQ" ) {
    $cmd .= " fastNLO/trunk/v1.4/author1c/common/ctq61.00.tbl".
	" fastNLO/trunk/v1.4/author1c/hadron/common".
	" fastNLO/trunk/v1.4/author1c/hadron/ctq61.00.tbl";
} elsif ( $pdf eq "LHAPDF" ) {
    $cmd .= " lhapdf/lib";
}
my $ret = system("$cmd");
if ( $ret ) {die "fastprep.pl: ERROR! Could not create archive fastNLO-bin.tgz: $ret\n";}
exit 0;
