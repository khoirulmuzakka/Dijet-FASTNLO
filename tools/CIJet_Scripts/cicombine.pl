#!/usr/bin/env perl
#
# created by K. Rabbertz: 16.10.2016
# last modified:
#
# modification log:
#
#-----------------------------------------------------------------------
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

my $ncount = 0;
my @massbins = ("MassBin-3600-4200","MassBin-4200-4800","MassBin-4800-5400","MassBin-5400-6000","MassBin-6000-6600","MassBin-6600-13000");
my @chibins = ("ChiBin-01-02","ChiBin-02-03","ChiBin-03-04","ChiBin-04-05","ChiBin-05-06","ChiBin-06-07",
               "ChiBin-07-08","ChiBin-08-09","ChiBin-09-10","ChiBin-10-12","ChiBin-12-14","ChiBin-14-16");

foreach my $massbin ( @massbins ) {
    foreach my $chibin ( @chibins ) {
        my $fglob = "CIJET_".$massbin."_".$chibin."_app1_"."*".".sum";
        my @files = glob("$fglob");

        my $nfiles = @files;
        if ( ! $nfiles ) {
            print "cicombine.pl: Could not find filename matching $fglob, skipped!\n";
            next;
        }
        if ( $ncount == 0 ) {
            $ncount = $nfiles;
        } else {
            if ( $ncount != $nfiles ) {die "cicombine.pl: Found unequal number of files per mass- and chi-bin. Cannot combine all files, please check your setup. Stopped!\n";}
        }
        print "cicombine.pl: Found $nfiles files matching $fglob\n";
        if ( $nfiles > 999 ) {die "cicombine.pl: Found too many files per mass- and chi-bin. Stopped!\n"}
        for ( my $i=0; $i<$nfiles; $i++ ) {
            my $no = sprintf("%03d",$i);
            my $infile  = $files[$i];
            open(INFILE,"< $infile") || die "cicombine.pl: Cannot open file $infile, aborted!\n";
            my $outfile = "CIJET_app1_".$no.".sum";
            open(OUTFILE,">> $outfile") || die "cicombine.pl: Cannot open file $outfile, aborted!\n";
            while ( my $in = <INFILE> ) {
                print OUTFILE $in;
            }
            close INFILE;
            close OUTFILE;
        }
    }
}

exit;
