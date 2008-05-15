#!/usr/bin/env perl 
#
# Check on identical fastNLO tables (aka identical random seeds ...)
# Version:
# 
# modified by K. Rabbertz from first version of T. Kluge: 13.03.2006
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

my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n######################################################\n";
print "# fastidcheck.pl: Starting table id check for fastNLO: FASTIDCHECK_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "fastidcheck.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastidcheck.pl\n";
    print "Usage: fastidcheck.pl glob (selects all files matching glob)\n";
    print "  -h              Print this text\n\n";
    exit;
}

#
# Initialization
#
my $glstr = shift;
chomp $glstr;
my @files = glob "*${glstr}*";  
chomp @files;

#
# Diff
#
for ( my $i=0; $i < @files; $i++) {
    print "fastidcheck.pl: Checking $files[$i] ...\n";
    for( my $j=$i+1; $j < @files; $j++){
	if (system("diff -q $files[$i] $files[$j] > /dev/null") == 0 ) {
	    print "fastidcheck.pl: Identical: $files[$i] and $files[$j]\n";
	}
    }
}

exit 0;
