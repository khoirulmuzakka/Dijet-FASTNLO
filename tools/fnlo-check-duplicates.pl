#!/usr/bin/env perl
#
# Check on identical fastNLO tables (aka identical random seeds ...)
# Version:
#
# created by T. Kluge: 13.03.2006
# modified by K. Rabbertz: 14.03.2006
# adapted by K. Rabbertz from fastidcheck.pl: 23.04.2015
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
if ( $? ) {die "fnlo-check-duplicates.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################################\n";
print "# fnlo-check-duplicates.pl: Starting table id check for fastNLO: FASTIDCHECK_$date\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_v ) = ( "", "2.3" );
getopts('hv:') or die "fnlo-check-duplicates.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-check-duplicates.pl\n";
    print "Usage: fnlo-check-duplicates.pl glob (selects all files matching glob)\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-check-duplicates.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-check-duplicates.pl: Error! Need one file name glob to select files to check!\n";
}
my $glstr = shift;
chomp $glstr;
print "fnlo-check-duplicates.pl: Checking for file glob $glstr ...\n";
my @files = glob "*${glstr}*tab*";
chomp @files;
if ( ! -d "Duplicates" ) {
    print "fnlo-check-duplicates.pl: Creating subdirectory Duplicates ...\n";
    system("mkdir Duplicates");
}

#
# Diff
#
for ( my $i=0; $i < @files; $i++) {
# File $i has not been moved already into Duplicates ...
    if ( -e $files[$i] ) {
        print "fnlo-check-duplicates.pl: Checking $files[$i] ...\n";
        for( my $j=$i+1; $j < @files; $j++){
# File $j has not been moved already into Duplicates ...
            if ( -e $files[$j] ) {
                if (system("zdiff -q $files[$i] $files[$j] > /dev/null") == 0 ) {
                    print "fnlo-check-duplicates.pl: Identical: $files[$i] and $files[$j]\n";
                    print "fnlo-check-duplicates.pl: Moving $files[$j] into Duplicates\n";
		    my $base = $files[$j];
		    $base =~ s/\.gz$//;
		    $base =~ s/\.tab$//;
		    my @prbfiles = glob "${base}*";
		    chomp @prbfiles;
		    foreach my $prbfile ( @prbfiles ) {
			my $ret = system("mv $prbfile Duplicates");
			if ( $ret ) {die "fnlo-check-duplicates.pl: Couldn't move file $prbfile into ".
					 "Problems: $ret, aborted!\n";}
		    }
                }
            }
        }
    }
}

exit 0;
