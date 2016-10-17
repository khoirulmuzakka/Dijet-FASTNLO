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



my @dirs = glob("job_*");
foreach my $jobdir ( @dirs ) {
    print "cijobrename.pl: Accessing jobdir $jobdir\n";
    my $jobno = $jobdir;
    $jobno =~ s/job_//;
    $jobno  = sprintf("%03u", $jobno);
#    print "jobno is $jobno\n";
    chdir $jobdir;
    my @cifiles = glob("CIJET*");
    foreach my $filnam ( @cifiles ) {
        my $newnam = $filnam;
        $newnam =~ s/\.dat$//;
        $newnam =~ s/\.sum$//;
        if (! ($newnam =~ m/$jobno$/)) {
            $newnam =~ s/\.//g;
            if ( ($newnam =~ m/ChiBin-[1-9]-[1-9]_/) ||
                 ($newnam =~ m/ChiBin-[1-9][0-9]-[1-9]_/) ||
                 ($newnam =~ m/ChiBin-[1-9]-[1-9][0-9]_/) ) {
                my @frags = split("_",$newnam);
                my $chibin = $frags[2];
                my @chifrags = split("-",$chibin);
                my $chilow = $chifrags[1];
                my $chiupp = $chifrags[2];
                $chilow = sprintf("%02u", $chilow);
                $chiupp = sprintf("%02u", $chiupp);
                my $chinew = $chifrags[0]."-".$chilow."-".$chiupp;
                $newnam =~ s/$chibin/$chinew/;
            }
            if ($filnam =~ m/dat/) {
                $newnam = $newnam."_".$jobno.".dat";
            } elsif ($filnam =~ m/sum/) {
                $newnam = $newnam."_".$jobno.".sum";
            }
            print "cijobrename.pl: Rename $filnam to $newnam\n";
#            print "cijobrename.pl: Test run, uncomment execution line below ...\n";
            rename $filnam, $newnam;
        } elsif ( ($newnam =~ m/ChiBin-[1-9]-[1-9]_/) ||
                  ($newnam =~ m/ChiBin-[1-9][0-9]-[1-9]_/) ||
                  ($newnam =~ m/ChiBin-[1-9]-[1-9][0-9]_/) ) {
            my @frags = split("_",$newnam);
            my $chibin = $frags[2];
            my @chifrags = split("-",$chibin);
            my $chilow = $chifrags[1];
            my $chiupp = $chifrags[2];
            $chilow = sprintf("%02u", $chilow);
            $chiupp = sprintf("%02u", $chiupp);
            my $chinew = $chifrags[0]."-".$chilow."-".$chiupp;
            $newnam =~ s/$chibin/$chinew/;
            if ($filnam =~ m/dat/) {
                $newnam = $newnam.".dat";
            } elsif ($filnam =~ m/sum/) {
                $newnam = $newnam.".sum";
            }
            print "cijobrename.pl: Rename $filnam to $newnam\n";
            rename $filnam, $newnam;
        }
    }
    chdir "..";
}

exit (0);
