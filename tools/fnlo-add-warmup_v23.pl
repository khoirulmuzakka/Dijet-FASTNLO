#!/usr/bin/env perl
#
# fastNLO warmup file addition script
#
# Version:
#
# created by K. Rabbertz: 05.03.2006
# last modified:
# 02.05.2011 KR: Added warmup table part
# adapted by K. Rabbertz from fastadd.pl: 29.07.2014
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
print "# fnlo-add-warmup.pl: Starting warmup file addition for fastNLO: WARMADD_$date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h, $opt_o, $opt_v, $opt_w ) = ( "", "", , "", "2.3", "" );
getopts('dho:v:w:') or die "fnlo-add-warmup.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "fnlo-add-warmup.pl\n";
    print "Usage: fnlo-add-warmup.pl [switches/options] glob (selects all files matching glob)\n";
    print "       Choose glob to match e.g. files for all orders but only one table type,\n";
    print "       e.g. for 2.4 the NNLOJET gridname: ZJtriple_yb0_ystar0_ptz\n";
    print "  -d              Verbose output\n";
    print "  -h              Print this text\n";
    print "  -o outfile      Output filename, (def.=glob_warmup.txt|wrm)\n";
    print "  -v #            Choose between fastNLO version 2.3 or 2.4 (def.=2.3)\n";
    print "  -w dir          Directory for warmup files, (def.=.)\n\n";
    exit;
}

#
# Parse arguments
#
if ( @ARGV != 1 ) {
    die "fnlo-add-warmup.pl: Error! Need exactly one filename glob!\n";
}
my $glob   = shift;
my $debug  = $opt_d;
my $vers   = $opt_v;

#
# Initialization
#
if (($vers != 2.3) && ($vers != 2.4)) {die "fnlo-add-warmup.pl: Error! Unsupported warmup file version: $vers\n"};
my $outfile = "${glob}_warmup.txt";
if ( $vers == 2.4 ) {$outfile = "${glob}_warmup.wrm";}
if ( $opt_o ) {$outfile = $opt_o;}
my $wdir    = ".";
if ( $opt_w ) {$wdir = $opt_w;}

my $wrmglob = "*${glob}*.txt";
if ( $vers == 2.4 ) {
    $wrmglob = "*${glob}*.wrm";
}

# Directory
my $sdir = getcwd();

#
# Analyze warmup files
#
print "\nfnlo-add-warmup.pl: Analyzing warmup files\n";
if ( -d "$wdir" ) {
    chdir $wdir;
    my @files = glob $wrmglob;
    chomp @files;
    unless ( @files ) {
        die "fnlo-add-warmup.pl: ERROR! No warm-up files found for filename glob $glob, aborted.\n";
    }
    my $ifil  = 0;
    my $conti = "# This file has been calculated using";
    my $contf = "contributions.";
    my $nentf = "entries.";
    my $contr = 0;
    my $nentr = 0;
    my $nent  = 0;
    my $nobs  = 0;
    my $nordi = "OrderInAlphas";
    my $nord  = 0;
    my $stat  = 0;
    my $wrm0  = 0;
    my $wrm1  = 0;
    my $wrm2  = 0;
    my $wrm3  = 0;
    my $nscl  = 0;
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
# Do not reuse already presummed files matching "warmup"
        if ( $file =~ m/warmup/ ) {next;}
        print "fnlo-add-warmup.pl: Opening file: $file\n";
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
            } elsif ( $vers < 2.4 && $in =~ $conti ) {
                my $tmp = $in;
                $tmp =~ s/^\s+//;
                chomp $tmp;
                my @tmps = split(/\s+/,$tmp);
                $contr = $contr + $tmps[7];
            } elsif ( $vers > 2.3 && $in =~ m/${contf}$/ ) {
                my $tmp = $in;
                chomp $tmp;
                my @tmps = split(/\s+/,$tmp);
                $contr = $contr + $tmps[1];
            } elsif ( $vers > 2.3 && $in =~ m/${nentf}$/ ) {
                my $tmp = $in;
                chomp $tmp;
                my @tmps = split(/\s+/,$tmp);
                $nentr = $nentr + $tmps[1];
            } elsif ( $in =~ m/${nordi}/ ) {
                my $tmp = $in;
                chomp $tmp;
                my @tmps = split(/\s+/,$tmp);
                $nord = max($nord,$tmps[1]);
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
                if ( $vers > 2.3 && $cols == 5 ) {
                    $ifill[$ient2] = defined $ifill[$ient2] ? $ifill[$ient2] + $tmps[4] : $tmps[4];
                    $ient2++;
                }
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
        if ( $vers < 2.4 && $in =~ $conti ) {
            printf(OUTFILE "%s %7.1e %s\n",$conti,$contr,$contf);
        } elsif ( $vers > 2.3 && $in =~ $conti ) {
            print OUTFILE $in;
        } elsif ( $in =~ m/${contf}$/ ) {
            printf(OUTFILE "#      %d $contf\n", $contr);
        } elsif ( $in =~ m/${nentf}$/ ) {
            printf(OUTFILE "#      %d $nentf\n", $nentr);
        } elsif ( $in =~ m/${nordi}/ ) {
            my $tmp = $in;
            chomp $tmp;
            my @tmps = split(/\s+/,$tmp);
            printf(OUTFILE "$tmps[0]   %4i           # Maximum value found in files\n", $nord);
        } elsif ( $in =~ "Warmup.Values" ) {
            print OUTFILE $in;
            $wrm0++;
        } elsif ( $in =~ "ObsBin" && $wrm0 ) {
            print OUTFILE $in;
            $wrm1++;
        } elsif ( $wrm1 && !($in =~ "}")) {
            if ($nscl == 1) {
                printf(OUTFILE "   %4d     %9.2e  %9.2e  %16.4f  %16.4f\n",
                       $iobs[$ient],$xmin[$ient],$xmax[$ient],$pmin[$ient],$pmax[$ient]);
            } elsif ($nscl == 2) {
                printf(OUTFILE "   %4d    %9.2e  %9.2e  %14.6f  %14.6f  %14.6f  %14.6f\n",
                       $iobs[$ient],$xmin[$ient],$xmax[$ient],$pmin[$ient],$pmax[$ient],$qmin[$ient],$qmax[$ient]);
            }
            $ient++;
        } elsif ( $in =~ "Warmup.Binning" ) {
            print OUTFILE $in;
            $wrm2++;
        } elsif ( $in =~ "ObsBin" && $wrm2 ) {
            print OUTFILE $in;
            $wrm3++;
        } elsif ( $wrm3 && !($in =~ "}")) {
            if ( $vers < 2.4 ) {
                print OUTFILE $in;
            } else {
                my $tmp = $in;
                $tmp =~ s/^\s+//;
                chomp $tmp;
                my @tmps = split(/\s+/,$tmp);
                my $cols = scalar @tmps;
                if ( $cols == 5 ) {
                    printf(OUTFILE "   %4d   %12.3f  %12.3f  %12.4f  %16d\n",
                           $tmps[0],$tmps[1],$tmps[2],$tmps[3],$ifill[$ient2]);
                }
                $ient2++;
            }
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
#    my $ret = system("mv $outfile $sdir");
#    if ( $ret ) {print "fnlo-add-warmup.pl: Couldn't move warmup summary file into ".
#                     "work directory $sdir: $ret. Please look for file in $wdir.\n";}
#    chdir $sdir;
}

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-warmup.pl: Warmup file addition result file is: $outfile\n";
print "fnlo-add-warmup.pl: Finished warmup file addition: WARMEND_$date\n";
exit(0);
