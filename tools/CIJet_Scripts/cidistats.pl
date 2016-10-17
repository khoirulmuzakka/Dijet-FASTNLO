#!/usr/bin/env perl
#
# created by K. Rabbertz: 26.03.2014
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
if (! ($massbin || $lambda || $type) ) {die "cidistats.pl: Please specify mass bin, lambda, and type. Stopped!\n";}

my @chibins = ("ChiBin-01-02","ChiBin-02-03","ChiBin-03-04","ChiBin-04-05","ChiBin-05-06","ChiBin-06-07",
               "ChiBin-07-08","ChiBin-08-09","ChiBin-09-10","ChiBin-10-12","ChiBin-12-14","ChiBin-14-16");

my @orders = ("Order-0","Order-1");
my %pdfs;
$pdfs{$orders[0]}[0] = "cteq6ll";
$pdfs{$orders[1]}[0] = "cteq66";
my %types;
$types{CILHC}[0] = "LL-";
$types{DILHC}[0] = "LL+";

foreach my $ord ( @orders ) {

    my $filout = "CIDIJET_MassBin-".$massbin."_AllChiBins_".$pdfs{$ord}[0]."_".$lambda."_".$types{$type}[0].".txt";
    open(OUTFILE,"> $filout") || die "cidistats.pl: Cannot open file $filout, aborted!\n";

    my $fglob = "cidijet_DijetChi_".$type."_2012_Lambda-".$lambda."._MassBin-".$massbin."*".$ord."*.dat";
    my @files = glob("$fglob");

    my $nfiles = @files;
    if ( ! $nfiles ) {
        print "cidistats.pl: Could not find filename matching $fglob, skipped!\n";
        next;
    }
    if ($nfiles > 1) {die "cidistats.pl: Found more than one matching filename, please check your setup. Stopped!\n";}

    my $filnam = $files[0];
    open(INFILE,"< $filnam") || die "cidistats.pl: Cannot open file $filnam, aborted!\n";

    my $nxs  = 1;
    my $xs   = 0;
    my $dxs  = 0;
    my $chilow = 0;
    my $chiupp = 0;

    my $nl = -1;
    while ( my $in = <INFILE> ) {
        $nl++;
        if ($nl < 3) {
            my $out;
            if ( $nl == 0 ) {
                print OUTFILE "#\n";
                print OUTFILE "#  $filnam\n";
                print OUTFILE "#\n";
                $out = "# ".$in;
                chomp $out;
                if ( $filnam =~ m/$orders[0]/ ) {
                    $out = $out."   ".$pdfs{$ord}[0]."   LO\n";
                } else {
                    $out = $out."   ".$pdfs{$ord}[0]."   NLO\n";
                }
                print OUTFILE $out;
            } elsif ( $nl == 1 ) {
                $out = "# ".$in;
                print OUTFILE $out;
            } elsif ( $nl == 2 ) {
                $out = "# ".$in;
                print OUTFILE $out;
                print OUTFILE "#\n";
                print OUTFILE "# Chi_low    Chi_upp   xsec[pb]         error[pb]        error[%]     N_jobs\n";
            }
        } else {
            my @numbers = split(" ",$in);
            for (my $i=0; $i <= $#numbers; $i++) {
                $numbers[$i] =~ s/D/E/;
            }
            $chilow = $numbers[0];
            $chiupp = $numbers[1];
            $xs     = $numbers[2];
            $dxs    = $numbers[3];

            my $xsm = $xs/$nxs;
            $dxs = $dxs;
            my $dxsrel = 100.*$dxs/abs($xsm);

            if ( $nl > 11 ) {
# Normalize to bin width 2 in Chi not done in CIJET
                $xsm = $xsm/2.;
                $dxs = $dxs/2.;
            }

            my $line = sprintf("%6.1f     %6.1f     %12.5e     %12.5e     %8.3f      %2u",$chilow,$chiupp,$xsm,$dxs,$dxsrel,$nxs);
            print OUTFILE "$line\n";
        }
    }
    close INFILE;
    close OUTFILE;
}

exit (0);
