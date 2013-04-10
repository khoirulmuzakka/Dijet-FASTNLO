#!/usr/bin/env perl
#
# Re-/Move fastNLO fluctuation tables
#
# created by K. Rabbertz: 09.04.2013
# last modified:
#
# modification log:
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

#
# Parse options
#
our ( $opt_h ) = ( "" );
getopts('h') or die "fastmv.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastmv.pl: Re-/Move fastNLO fluctuation tables\n";
    print "Usage: fastmv.pl killfilename\n";
    print "  -h              Print this text\n\n";
    exit;
}

#
# Initialization
#
my $infile = shift;
open(KILLIN,"< $infile") or die "fastmv.pl: Could not open $infile!\n";

#
# Read all lines
#
while ( my $in = <KILLIN> ) {
    my $file = $in;
    chomp $file;
    my @parts = split("/",$file);
    my $scenname;
    foreach my $part (@parts) {
        if ( $part eq "stat" ) {
            last;
        } else {
            $scenname = $part;
        }
    }
    my $cmd1 = "rm -f $file";
    print "Executing rm command: $cmd1\n";
    my $ret = system("$cmd1");
    if ( $ret ) {die "fastmv.pl: Re-/moving command failed: $ret, aborted!\n";}
    $file =~ s/results/tables/;
    my $sfil = $file;
    my $tfil = $file;
    $sfil =~ s/stat/$scenname/;
    my $flcdir = $scenname."/Fluctuations";
    $tfil =~ s/stat/$flcdir/;
    my $tdir = `dirname $tfil`;
    chomp $tdir;
    if (! -d $tdir) {
        my $ret = system("mkdir $tdir");
        if ( $ret ) {die "fastmv.pl: Creating target directory $tdir failed: $ret, aborted!\n";}
    }
    my $cmd2 = "mv $sfil $tfil";
    print "Executing mv command: $cmd2\n";
    $ret = system("$cmd2");
    if ( $ret ) {die "fastmv.pl: Re-/moving command failed: $ret, aborted!\n";}
}

print "fastmv.pl: Finished.\n";

exit 0;
