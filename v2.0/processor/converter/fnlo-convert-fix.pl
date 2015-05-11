#!/usr/bin/env perl
#
# Eliminate some residual differences between converted v1.4 tables and
# original v2.0 tables like leading/trailing blanks, trailing decimal points
# and inline space
#
# created by K. Rabbertz: 05.09.2009
# adapted by K. Rabbertz from tablefix.pl: 08.05.2015
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
if ( $? ) {die "fnlo-convert-fix.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################################\n";
print "# fnlo-convert-fix.pl: Fixing residual differences of tables converted from v1.4\n";
print "######################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_n, $opt_r ) = ( "", "", "" );
getopts('hn:r') or die "fnlo-convert-fix.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-convert-fix.pl\n";
    print "Usage: fnlo-convert-fix.pl table file\n";
    print "  -h              Print this text\n";
    print "  -n              Name of cleaned table (def.=\"oldname\".new)\n";
    print "  -r              Replace old table by new one\n";
    print "                  No further options available\n\n";
    exit;
}
if ( $opt_n && $opt_r ) {
    die "fnlo-convert-fix.pl: Error! Specifying both options, -n $opt_n and -r does not make sense!\n";
}
my $newname;
if ( $opt_n ) {$newname = $opt_n;}
my $lrep = $opt_r;


#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-convert-fix.pl: Error! No table file specified!\n";
}
my $infile = shift;
chomp $infile;
my $outfile = ${infile}.".new";
if ( $newname ) {
    $outfile = $newname;
} elsif ( $lrep ) {
    $outfile = $infile;
}
my $tmpout = ${outfile}.".tmp";
if ( -e $tmpout ) {die "fnlo-convert-fix.pl: Error! File $tmpout for temporary write out exists already, please remove first!\n";}
open(INFILE ,"< $infile") or die "fnlo-convert-fix.pl: Error! Could not open infile $infile!\n";
open(OUTFILE,"> $tmpout") or die "fnlo-convert-fix.pl: Error! Could not open temporary outfile $tmpout!\n";

while ( my $line = <INFILE> ) {
# Remove all leading blanks
    $line =~ s/^\ +//;
# Remove all trailing blanks
    $line =~ s/\ +$//;
# Remove a trailing decimal point
    if ( ! $line =~ m/Nagy/ ) {
        $line =~ s/\.$//;
    }
# Replace inline blanks by underscores
    if ( ! $line =~ m/Nagy/ ) {
        $line =~ s/\ /_/;
    }
    print OUTFILE "$line";
}
close INFILE;
close OUTFILE;

if ( -e $outfile ) {
    if ( $lrep ) {
        my $ret = system("rm -f $outfile");
        if ( $ret ) {die "fnlo-convert-fix.pl: Error! Removal of $outfile failed: $ret, aborted!\n";}
    } else {
        die "fnlo-convert-fix.pl: Error! $outfile exists already, please remove first!\n";
    }
}
my $ret = system("mv $tmpout $outfile");
if ( $ret ) {die "fnlo-convert-fix.pl: Error! Renaming $tmpout to $outfile failed: $ret, aborted!\n";}

print "fnlo-convert-fix.pl: Finished.\n";

exit;
