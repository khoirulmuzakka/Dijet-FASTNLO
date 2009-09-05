#!/usr/bin/env perl 
#
# Eliminate remaining differences between converted v1.4 tables and 
# original v2.0 tables like leading/trailing blanks, trailing decimal points
# and inline space
#
# created by K. Rabbertz: 05.09.2009
# last modified:
# 
# modification log:
#
#-----------------------------------------------------------------------
# Todo:
#

use Cwd;
use English;
use strict;
use warnings;

my $infile = shift;
open(INFILE,"< $infile") or die "tablefix.pl: Could not open infile $infile!\n";
open(OUTFILE,"> ${infile}.new") or die "tablefix.pl: Could not open outfile $infile!\n";

while ( my $line = <INFILE> ) {	
# Remove all leading blanks
    $line =~ s/^\ +//;
# Remove all trailing blanks
    $line =~ s/\ +$//;
# Remove a trailing decimal point
    $line =~ s/\.$//;
# Replace inline blanks by underscores
    $line =~ s/\ /_/;
    print OUTFILE "$line";
}
close INFILE;
close OUTFILE;

print "Done\n";

exit;
