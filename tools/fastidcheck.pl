#!/usr/bin/env perl 
#
# Check on identical fastNLO tables (aka identical random seeds ...)
# Version:
# 
# modified by K. Rabbertz from first version of T. Kluge: 13.03.2006
# last modified:
#
#-----------------------------------------------------------------------
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
our ( $opt_h, $opt_v ) = ( "", "1" );
getopts('hv:') or die "fastidcheck.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastidcheck.pl\n";
    print "Usage: fastidcheck.pl glob (selects all files matching glob)\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO version 1 or 2 (def.=1)\n\n";
    exit;
}

#
# Initialization
#
my $glstr = shift;
chomp $glstr;
my @files = glob "*${glstr}*";  
chomp @files;
if ( ! -d "Dubletten" ) {
    print "fastidcheck.pl: Creating subdirectory Dubletten ...\n";
    system("mkdir Dubletten");
}
my $vers  = $opt_v;
my $tabext = "raw";
if ($vers == 2) {$tabext = "tab";} 

#
# Diff
#
for ( my $i=0; $i < @files; $i++) {
# File $i has not been moved already into Dubletten ...
    if ( -e $files[$i] ) {
	print "fastidcheck.pl: Checking $files[$i] ...\n";
	for( my $j=$i+1; $j < @files; $j++){
# File $j has not been moved already into Dubletten ...
	    if ( -e $files[$j] ) {
		if (system("diff -q $files[$i] $files[$j] > /dev/null") == 0 ) {
		    print "fastidcheck.pl: Identical: $files[$i] and $files[$j]\n";
		    print "fastidcheck.pl: Moving $files[$j] into Dubletten\n";
		    my $logfil = $files[$j];
		    $logfil =~ s/${tabext}/log/;
#		    my $errfil = $files[$j];
#		    $errfil =~ s/${tabext}/err/;
		    system("mv $files[$j] Dubletten");
		    system("mv $logfil Dubletten");
#		    system("mv $errfil Dubletten");
		}
	    }
	}
    }
}

exit 0;
