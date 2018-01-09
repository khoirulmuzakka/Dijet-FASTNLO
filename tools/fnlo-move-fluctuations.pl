#!/usr/bin/env perl
#
# Re-/Move fastNLO fluctuation tables:
#    Move tables with largest statistical deviation in any observable bin
#    from default location into 'Fluctuations' subfolder. Also delete
#    corresponding temporary table for statistical evaluation.
#    The table locations are written into a 'kill' file by the
#    evaluation job. CHECK that your directory setup matches!
#
# Version:
#
# created by K. Rabbertz: 09.04.2013
# last modified:
# adapted by K. Rabbertz from fastmv.pl: 02.06.2015
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
chomp $date;
print "\n##################################################\n";
print "# fnlo-move-fluctuations.pl: Starting table removal: TABDEL_$date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "fnlo-move-fluctuations.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-move-fluctuations.pl: Re-/Move fluctuating fastNLO tables\n";
    print "Usage: fnlo-move-fluctuations.pl killfilename\n";
    print "  -h              Print this text\n";
    print "\n";
    print "Absolutely check the kill file first for the table locations and\n";
    print "compatibility with your directory setup!\n";
    exit;
}

#
# Initialization
#
my $infile = shift;
open(KILLIN,"< $infile") or die "fnlo-move-fluctuations.pl: Could not open $infile!\n";

#
# Read all lines
#
while ( my $in = <KILLIN> ) {
    my $file = $in;
    chomp $file;
    if ( ! -f $file ) {
	$file .= ".gz";
    }
    if ( ! -f $file ) {
	die "File $file not found, aborted!";
    }
    print "Found file $file.\n";

    my @parts = split("/",$file);
    my $scenname;
    foreach my $part (@parts) {
        if ( $part eq "stat" ) {
            last;
        } else {
            $scenname = $part;
        }
    }
    print "Derived scenario to be: $scenname\n";
    my $cmd1 = "rm -f $file";
    print "Executing rm command: $cmd1\n";
    my $ret = 0;
    my $ret = system("$cmd1");
    if ( $ret ) {die "fnlo-move-fluctuations.pl: Re-/moving command failed: $ret, aborted!\n";}
    $file =~ s/results/tables/;
    my $sfil = $file;
    my $tfil = $file;
    $sfil =~ s/stat/$scenname/;
    my $flcdir = $scenname."/Fluctuations";
    $tfil =~ s/stat/$flcdir/;
    my $tdir = `dirname $tfil`;
    chomp $tdir;
    print "tdir $tdir\n";
    if (! -d $tdir) {
        my $ret = system("mkdir $tdir");
        if ( $ret ) {die "fnlo-move-fluctuations.pl: Creating target directory $tdir failed: $ret, aborted!\n";}
    }
    my $cmd2 = "mv $sfil $tfil";
    print "Executing table mv command: $cmd2\n";
    $ret = system("$cmd2");
    if ( $ret ) {die "fnlo-move-fluctuations.pl: Re-/moving command failed: $ret, aborted!\n";}
    $sfil =~ s/\.tab\./\.log\./;
    $tfil =~ s/\.tab\./\.log\./;
    $cmd2 = "mv $sfil $tfil";
    if (-f $sfil) {
        print "Executing log file mv command: $cmd2\n";
        $ret = system("$cmd2");
        if ( $ret ) {die "fnlo-move-fluctuations.pl: Re-/moving command failed: $ret, aborted!\n";}
    }
}

print "fnlo-move-fluctuations.pl: Finished.\n";

exit 0;
