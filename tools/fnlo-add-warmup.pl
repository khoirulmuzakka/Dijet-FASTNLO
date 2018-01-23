#!/usr/bin/env perl
#
# fastNLO warmup table addition script
#
# Version:
#
# created by K. Rabbertz: 05.03.2006
# last modified:
# 02.05.2011 KR: Added warmup table part
# adapted by K. Rabbertz from fastadd.pl: 29.07.2014
# generalised by KR to filename*_#####.wrm: 15.11.2017
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use List::Util qw[min max];
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n##################################################\n";
print "# fnlo-add-warmup.pl: Starting warmup table addition for fastNLO: WARMADD_$date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_e, $opt_h, $opt_o, $opt_v, $opt_w ) = ( "", "wrm", "", "", "2.4", "" );
getopts('de:ho:v:w:') or die "fnlo-add-warmup.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "fnlo-add-warmup.pl\n";
    print "Usage: fnlo-add-warmup.pl [switches/options] basefilename\n";
    print "       where filenames are basefilename*_nnnn.wrm\n";
    print "  -d              Verbose output\n";
    print "  -e              Warmup file extension. (def.=wrm)\n";
    print "  -h              Print this text\n";
    print "  -o outfile      Output filename, (def.=basefilename*.wrm)\n";
    print "  -v #            Choose between fastNLO versions, only 2.4 for now. (def.=2.4)\n";
    print "  -w dir          Directory for warmup tables, (def.=.)\n\n";
    exit;
}

#
# Parse arguments
#
if ( @ARGV != 1 ) {
    die "fnlo-add-warmup.pl: Error! Need exactly one basefilename!\n";}
my $base    = shift;
my $debug   = $opt_d;
my $ext     = $opt_e;
my $vers    = $opt_v;

#
# Initialization
#
if ( $vers != 2.4 ) {die "fnlo-add-warmup.pl: Error! Unsupported warmup file version: $vers\n"};
my $outfile = "${base}.${ext}";
if ( $opt_o ) {$outfile = $opt_o;}
my $wdir    = getcwd();
if ( $opt_w ) {$wdir = $opt_w;}
my $wrmglob = "${base}*_????.${ext}";

# Directory
my $sdir = getcwd();
chdir $wdir or die "fnlo-add-warmup.pl: ERROR! Could not cd to directory $wdir, aborted.\n";

#
# Analyze warmup tables
#
print "\nfnlo-add-warmup.pl: Analyzing warmup tables in directory $wdir\n";
my @files = glob $wrmglob;
chomp @files;
unless ( @files ) {
    die "fnlo-add-warmup.pl: Warning! No warmup files found for filename glob $wrmglob, stopped\n";
}
my $ifil = 0;
my $conti = "# This file has been calculated using";
my $contf = "contributions.";
my $nentf = "entries.";
my $contr = 0;
my $nentr = 0;
my $nent = 0;
my $nobs = 0;
my $stat = 0;
my $wrm0 = 0;
my $wrm1 = 0;
my $wrm2 = 0;
my $wrm3 = 0;
my $nscl = 0;
my @iobs;
my @xmin;
my @xmax;
my @pmin;
my @pmax;
my @qmin;
my @qmax;
my @ifill;

#
# Loop over all files determining min and max values, count contributions
#
foreach my $file ( @files ) {
    if ( $debug ) {print "fnlo-add-warmup.pl: Opening file: $file\n";}
    open(INFILE,"< $file") or die "fnlo-add-warmup.pl: Error! Could not open $file!\n";
    my $ient  = 0;
    my $ient2 = 0;
    if ( $debug ) {print "fnlo-add-warmup.pl: Analyzing file no.: $ifil\n";}
    while ( my $in = <INFILE> ) {
        if ( $debug ) {print "fnlo-add-warmup.pl: Line to analyze is: $in";}
        if ( $in =~ m/steerfile\:/ ) {
            my $tmp = $in;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
        } elsif ( $in =~ m/${contf}$/ ) {
            my $tmp = $in;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
            $contr = $contr + $tmps[1];
        } elsif ( $in =~ m/${nentf}$/ ) {
            my $tmp = $in;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
            $nentr = $nentr + $tmps[1];
        } elsif ( $in =~ "Warmup.Values" ) {
            $wrm0++;
        } elsif ( $in =~ "ObsBin" && $wrm0 ) {
            $wrm1++;
        } elsif ( $wrm1 && !($in =~ "}")) {
            my $tmp = $in;
            $tmp =~ s/^\s+//;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
            my $cols = scalar @tmps;
            if ( $cols < 5 ) {
                die "fnlo-add-warmup.pl: Error! Not enough columns found for iobs, x and mu1 limits: $cols\n";
            } elsif ( $cols == 5 || $cols == 7 ) {
                if (!$nscl) {$nscl = 1};
                $iobs[$ient] = $tmps[0];
                if ($tmps[1] > 1.e-6) {
                    $xmin[$ient] = defined $xmin[$ient] ? min($xmin[$ient],$tmps[1]) : min(1.0,$tmps[1]);
                }
                $xmax[$ient] = defined $xmax[$ient] ? max($xmax[$ient],$tmps[2]) : max(0.0,$tmps[2]);
                $pmin[$ient] = defined $pmin[$ient] ? min($pmin[$ient],$tmps[3]) : min(+1e10,$tmps[3]);
                $pmax[$ient] = defined $pmax[$ient] ? max($pmax[$ient],$tmps[4]) : max(-1e10,$tmps[4]);
                if ( $cols == 7 ) {
                    if (!$nscl || $nscl == 1) {$nscl = 2};
                    $qmin[$ient] = defined $qmin[$ient] ? min($qmin[$ient],$tmps[5]) : min(+1e10,$tmps[5]);
                    $qmax[$ient] = defined $qmax[$ient] ? max($qmax[$ient],$tmps[6]) : max(-1e10,$tmps[6]);
                }
                $ient++;
            } else {
                die "fnlo-add-warmup.pl: Error! Incorrect no. of columns found for iobs, x, mu1, and mu2 limits: $cols\n";
            }
        } elsif ( $in =~ "Warmup.Binning" ) {
            $wrm2++;
        } elsif ( $in =~ "ObsBin" && $wrm2 ) {
            $wrm3++;
        } elsif ( $wrm3 && !($in =~ "}")) {
            my $tmp = $in;
            $tmp =~ s/^\s+//;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
            my $cols = scalar @tmps;
            $ifill[$ient2] = defined $ifill[$ient2] ? $ifill[$ient2] + $tmps[$cols-1] : $tmps[$cols-1];
            $ient2++;
        } else {
            $wrm0 = 0;
            $wrm1 = 0;
            $wrm2 = 0;
            $wrm3 = 0;
        }
    }
    if ( !$ifil ) {
        $nent = $ient;
        $nobs = $iobs[$ient-1]+1;
    }
    if ( $debug ) {print "fnlo-add-warmup.pl: Number of entries: $nent\n";}
    if ( $debug ) {print "fnlo-add-warmup.pl: Number of observable bins: $nobs\n";}
    if ( $ifil && ($ient != $nent || $nent != $nobs) ) {
        print "fnlo-add-warmup.pl: Error! Inconsistent line numbers in warm-up files found \n";
        die "            (ient = $ient, nent = $nent, nobs = $nobs), stopped\n";
    }
    close INFILE;
    $ifil++;
}

#
# Loop over all lines of first file replacing limits with determined min and max values
#
open(INFILE,"< $files[0]") or die "fnlo-add-warmup.pl: Error! Could not open $files[0]!\n";
if ( $debug ) {print "fnlo-add-warmup.pl: Opening file $files[0] for replacing min/max values\n";}
open(OUTFILE,"> $outfile") or die "fnlo-add-warmup.pl: Error! Could not open $outfile!\n";
$wrm0 = 0;
$wrm1 = 0;
$wrm2 = 0;
$wrm3 = 0;
my $ient  = 0;
my $ient2 = 0;
while ( my $in = <INFILE> ) {
    if ( $in =~ $conti ) {
        print OUTFILE $in;
        if ( $debug ) {print "AAA: $in";}
    } elsif ( $in =~ m/${contf}$/ ) {
        printf(OUTFILE "#      %d $contf\n", $contr);
        if ( $debug ) {print "BBB\n";}
    } elsif ( $in =~ m/${nentf}$/ ) {
        printf(OUTFILE "#      %d $nentf\n", $nentr);
        if ( $debug ) {print "CCC\n";}
    } elsif ( $in =~ "Warmup.Values" ) {
        print OUTFILE $in;
        $wrm0++;
        if ( $debug ) {print "DDD wrm0 $wrm0\n";}
    } elsif ( $in =~ "ObsBin" && $wrm0 ) {
        print OUTFILE $in;
        $wrm1++;
        if ( $debug ) {print "EEE wrm1 $wrm1\n";}
    } elsif ( $wrm1 && !($in =~ "}")) {
        if ($nscl == 1) {
            printf(OUTFILE "   %4d     %9.2e  %9.2e  %16.4f  %16.4f\n",
                   $iobs[$ient],$xmin[$ient],$xmax[$ient],$pmin[$ient],$pmax[$ient]);
            if ( $debug ) {print "FFF\n";}
        } elsif ($nscl == 2) {
            printf(OUTFILE "   %4d    %9.2e  %9.2e  %14.6f  %14.6f  %14.6f  %14.6f\n",
                   $iobs[$ient],$xmin[$ient],$xmax[$ient],$pmin[$ient],$pmax[$ient],$qmin[$ient],$qmax[$ient]);
            if ( $debug ) {print "GGG\n";}
        }
        $ient++;
    } elsif ( $in =~ "Warmup.Binning" ) {
        print OUTFILE $in;
        $wrm2++;
        if ( $debug ) {print "HHH wrm2 $wrm2\n";}
    } elsif ( $in =~ "ObsBin" && $wrm2 ) {
        print OUTFILE $in;
        $wrm3++;
        if ( $debug ) {print "III wrm3 $wrm3\n";}
    } elsif ( $wrm3 && !($in =~ "}")) {
        if ( $debug ) {print "KKK\n";}
        my $tmp = $in;
        $tmp =~ s/^\s+//;
        chomp $tmp;
        my @tmps = split(/\s+/,$tmp);
        my $cols = scalar @tmps;
        $in =~ s/$tmps[$cols-1]$/$ifill[$ient2]/;
#        printf(OUTFILE "   %4d   %12.3f  %12.3f  %12.4f  %16d\n",
#               $tmps[0],$tmps[1],$tmps[2],$tmps[3],$ifill[$ient2]);
        print OUTFILE $in;
        $ient2++;
    } else {
        $wrm0 = 0;
        $wrm1 = 0;
        $wrm2 = 0;
        $wrm3 = 0;
        print OUTFILE $in;
    }
}
close OUTFILE;
close INFILE;
if ( ! -f "$sdir/$outfile" ) {
    my $ret = system("mv $outfile $sdir");
    if ( $ret ) {print "fnlo-add-warmup.pl: Couldn't move warmup summary file into ".
                     "work directory $sdir: $ret. Please look for file in $wdir.\n";}
} else {
    print "fnlo-add-warmup.pl: ATTENTION! Destination file $sdir/$outfile exists already!\n";
    print "                    Please remove this file first!\n";
}
chdir $sdir or die "fnlo-add-warmup.pl: ERROR! Could not cd back to directory $sdir, aborted.\n";

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-warmup.pl: Warmup table addition result file is: $outfile\n";
print "fnlo-add-warmup.pl: Finished warmup table addition: WARMEND_$date\n";
exit(0);
