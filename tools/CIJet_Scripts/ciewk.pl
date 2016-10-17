#!/usr/bin/env perl
#
# created by K. Rabbertz: 27.03.2014
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

my $fglob = "DijetAngularCMS-CT10nlo-8TeV_R0.5";
#my $fglob = "DijetAngularCMS-CT10nlo-8TeV_R0.7";
my @files = glob("${fglob}*");

my %kewk;
my @chibins;
my @chilows = ( 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 12, 14 );
my @chiupps = ( 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16 );
my $firstline;


foreach my $filnam ( @files ) {

# Skip already produced summary .txt files
    if ( $filnam =~ m/\.txt$/ ) {next;}
    my @frags = split("_",$filnam);
    my $chibin = $frags[2];
    push @chibins, $chibin;
    open(INFILE,"< $filnam") || die "ciewk.pl: Cannot open file $filnam, aborted!\n";
    my $nl = -1;
    while ( my $in = <INFILE> ) {
        $nl++;
        if ( $in =~ m/^\#/ ) {
            $firstline = $in;
            chomp $firstline;
        } else {
            my @numbers = split(" ",$in);
#            print "numbers @numbers\n";
            my $masslow = $numbers[0];
            my $massupp = $numbers[1];
            my $dtree   = $numbers[2];
            my $dloop   = $numbers[3];
            $masslow = sprintf("%4u",$masslow);
# Fill 2-dimensional hash with array content
            push(@{$kewk{$masslow}{$chibin}},$dtree);
            push(@{$kewk{$masslow}{$chibin}},$dloop);
        }
    }
    close INFILE;
}

foreach my $massbin ( sort keys %kewk ) {
    my $filout = $fglob."_MassBin-".$massbin."_AllChiBins.txt";
    open(OUTFILE,"> $filout") || die "ciewk.pl: Cannot open file $filout, aborted!\n";
    print OUTFILE "#\n";
    print OUTFILE "#  $fglob\n";
    print OUTFILE "#  Mass bin starting at $massbin GeV\n";
    print OUTFILE "#  All chi bins\n";
    print OUTFILE "#\n";
    print OUTFILE "$firstline\n";

    my $i = -1;
    foreach my $chibin ( sort keys %{$kewk{$massbin}} ) {
        $i++;
        my $line = sprintf("%6.1f     %6.1f     %12.5e     %12.5e",$chilows[$i],$chiupps[$i],$kewk{$massbin}{$chibin}[0],$kewk{$massbin}{$chibin}[1]);
        print OUTFILE "$line\n";
    }
    close OUTFILE;
}

exit (0);
