#!/usr/bin/env perl
#
# created by K. Rabbertz: 25.03.2014
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

my $massbin = shift;
my $lambda  = shift;
my $type    = shift;
if (! ($massbin || $lambda || $type) ) {die "cistats.pl: Please specify mass bin, lambda, and type. Stopped!\n";}

my @chibins = ("ChiBin-01-02","ChiBin-02-03","ChiBin-03-04","ChiBin-04-05","ChiBin-05-06","ChiBin-06-07",
               "ChiBin-07-08","ChiBin-08-09","ChiBin-09-10","ChiBin-10-12","ChiBin-12-14","ChiBin-14-16");

my @pdfs = ("cteq6ll","cteq66","CT10nlo");
foreach my $pdf ( @pdfs ) {

    my $filout = "CIJET_MassBin-".$massbin."_AllChiBins_".$pdf."_".$lambda."_".$type.".txt";
    open(OUTFILE,"> $filout") || die "cistats.pl: Cannot open file $filout, aborted!\n";

    foreach my $chibin ( @chibins ) {

        my $fglob = "CIJET_MassBin-".$massbin."*".$chibin."_app1_".$pdf."_".$lambda."_".$type.".xsc";
        my @files = glob("$fglob");

        my $nfiles = @files;
        if ( ! $nfiles ) {
            print "cistats.pl: Could not find filename matching $fglob, skipped!\n";
            next;
        }
        if ($nfiles > 1) {die "cistats.pl: Found more than one matching filename, please check your setup. Stopped!\n";}

        my $filnam = $files[0];
        open(INFILE,"< $filnam") || die "cistats.pl: Cannot open file $filnam, aborted!\n";

        my $nxs  = 0;
        my $xs   = 0;
        my $dxs  = 0;

        my $nl = -1;
        while ( my $in = <INFILE> ) {
            $nl++;
            if ($chibin eq "ChiBin-01-02" && $nl < 3) {
                my $out = $in;
                if ( $nl == 0 ) {
                    print OUTFILE "#\n";
                    my @frags = split("_",$filnam);
                    my $newnam = $frags[0]."_".$frags[1]."_".$frags[4]."_".$frags[5]."_".$frags[6];
                    print OUTFILE "#  $newnam\n";
                    print OUTFILE "#\n";
                    $out =~ s/0$//;
                    $out =~ s/ //g;
                    chomp $out;
                    if ( $filnam =~ m/$pdfs[0]/ ) {
                        $out = " ".$out."   LO\n";
                    } else {
                        $out = " ".$out."   NLO\n";
                    }
                }
                $out = "# ".$out;
                print OUTFILE $out;
            } elsif ( $in =~ m/^  1\.000  1\.000/ ) {
                $nxs++;
                my @numbers = split(" ",$in);
# LO PDF for LO result
                if ( $filnam =~ m/$pdfs[0]/ ) {
                    $xs  = $xs  + $numbers[2];
                    $dxs = $dxs + $numbers[2]*$numbers[2];
# NLO PDF for NLO result
                } elsif ( $filnam =~ m/$pdfs[1]/ || $filnam =~ m/$pdfs[2]/ ) {
                    $xs  = $xs  + $numbers[3];
                    $dxs = $dxs + $numbers[3]*$numbers[3];
# Undefined PDF
                } else {
                    die "cistats.pl: Check your PDF set selection for consistency with LO/NLO settings. Aborted!\n";
                }
            }
        }
        close INFILE;

        my $xsm  = $xs/$nxs;
        if ( $nxs < 2 ) {
            $dxs = $xsm;
        } else {
            $dxs  = 1/($nxs-1) * ($dxs - $xs*$xs/$nxs );
            if ($dxs < 0) {
                print "cistats.pl: Warning! Variance smaller than zero, taking absolute value!\n";
                print "cistats.pl: xsm = $xsm; dxs = $dxs; nxs = $nxs\n";
            }
            $dxs  = sqrt(abs($dxs));
        }
        my $dxsrel = 100.*$dxs/abs($xsm);

        if ( $chibin eq "ChiBin-01-02" ) {
            print OUTFILE "#\n";
            print OUTFILE "# Chi_low    Chi_upp   xsec[pb]         error[pb]        error[%]     N_jobs\n";
        } elsif ( $chibin eq "ChiBin-10-12" || $chibin eq "ChiBin-12-14" || $chibin eq "ChiBin-14-16" ) {
# Normalize to bin width 2 in Chi not done in CIJET
            $xsm = $xsm/2.;
            $dxs = $dxs/2.;
        }
        my @frags = split("-",$chibin);
        my $chilow = scalar $frags[1];
        my $chiupp = scalar $frags[2];
        my $line = sprintf("%6.1f     %6.1f     %12.5e     %12.5e     %8.3f      %2u",$chilow,$chiupp,$xsm,$dxs,$dxsrel,$nxs);
        print OUTFILE "$line\n";
    }
    close OUTFILE;
}

exit (0);
