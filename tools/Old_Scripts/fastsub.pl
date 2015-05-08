#!/usr/bin/env perl
#
# fastNLO submit script
# Version:
#
# created by K. Rabbertz: 07.02.2006
# last modified:
#
#-----------------------------------------------------------------------
# Todo:
#

use Cwd;
use English;
#use FindBin qw($Bin);
#use lib "$Bin/modules";
use Getopt::Std;
use strict;
use warnings;

#
# Start
#
my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n##################################################\n";
print "# fastsub.pl: Starting submission of fastNLO jobs: $date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_f, $opt_h, $opt_i,
      $opt_m, $opt_n, $opt_o, $opt_p, $opt_r, $opt_t ) =
    ( "LCG", ".", "0", "187", "", "0001",
      "0", "1", "LO", "CTEQ", "", "" );
getopts('b:d:e:f:hi:m:n:o:p:q:rt:') or die "fastsub.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastsub.pl\n";
    print "Usage: fastsub.pl [switches/options] CE||queue ([scenario])\n";
    print "  -b batch        Batch system to use: LCG (def.) or PBS\n";
    print "  -d dir          Installation directory (def.=.)\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -f rev          fastNLO revision to use (def.=187)\n";
    print "  -h              Print this text\n";
    print "  -i nnnn         First job number (def.=0001)\n";
    print "  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only\n";
    print "  -n n            Number of jobs to submit (def.=1)\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -p pdf          CTEQ parton densities (def.) or LHAPDF\n";
    print "  -r              Reference calculation incl. pdf access\n";
    print "  -t dir          Output target directory: ".
        "(def.= {scen}{ref}_{jobnr} with\n                  ".
        "ref. to working directory in fastNLO installation)\n\n";
    exit;
}

unless ( $opt_b eq "LCG" || $opt_b eq "PBS" ) {
    die "fastsub.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fastsub.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( $opt_f =~ m/\d+/ && $opt_f !~ m/\D+/ ) {
    die "fastsub.pl: Error! Illegal fastNLO revision number: $opt_f, aborted.\n";
}
unless ( $opt_i =~ m/\d{4}/ && $opt_i !~ m/\D+/ ) {
    die "fastsub.pl: Error! Illegal starting number (nnnn): $opt_i, aborted.\n";
}
unless ( $opt_m == 0 || $opt_m == 1 || $opt_m == 2 || $opt_m == 3 ) {
    die "fastsub.pl: Error! Illegal job mode: $opt_m, aborted.\n";
}
unless ( $opt_n =~ m/\d{1,3}/ && $opt_n !~ m/\D+/ && $opt_n > 0 && $opt_n < 1000 ) {
    die "fastsub.pl: Error! Wrong number of jobs (1-999): $opt_n, aborted.\n";
}
unless ( $opt_o eq "LO" || $opt_o eq "NLO" ) {
    die "fastsub.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( $opt_p eq "CTEQ" || $opt_p eq "LHAPDF" ) {
    die "fastsub.pl: Error! Illegal option -o $opt_o, aborted.\n";
}

my $batch = $opt_b;
my $idir  = $opt_d;
my $nmax  = $opt_e;
my $frev  = $opt_f;
my $jobnr = $opt_i;
my $mode  = $opt_m;
my $njobs = $opt_n;
my $order = $opt_o;
my $pdf   = $opt_p;
my $ref   = "";
print "fastsub.pl: Using batch system $batch\n";
print "fastsub.pl: Installation directory $idir\n";
print "fastsub.pl: Maximal event number: $nmax\n";
print "fastsub.pl: Using fastNLO revision $frev\n";
print "fastsub.pl: Starting number of jobs to generate: $jobnr\n";
print "fastsub.pl: Job mode is $mode\n";
print "fastsub.pl: Number of jobs to generate: $njobs\n";
print "fastsub.pl: Submitting for order $order.\n";
print "fastsub.pl: Submitting for pdf $pdf.\n";
if ( $opt_r ) {
    $ref = "ref";
    print "fastsub.pl: Running in reference mode.\n";
}

#
# Parse arguments
#
unless ( ($mode == 1 && @ARGV == 1) || @ARGV == 2 ) {
    die "fastsub.pl: Error! Need one Compute Element or queue name plus a scenario name for mode 0,2,3!\n";
}
my $cequ = shift;
my $scen = shift;
my $tdir = "./";
if ( $opt_t ) {
    $tdir  = "${opt_t}/";
    print "fastsub.pl: Output reference directory changed to ${tdir}\n";
}

#
# Initialization
#
my $host = `hostname`;
chomp $host;
# Run mode settings: order, ordername, # events
my %runmode;
$runmode{LO}[0]  = "born";
$runmode{LO}[1]  = "100000000";
$runmode{NLO}[0] = "nlo";
$runmode{NLO}[1] = "10000000";

#
# Create and submit files
#
my @infiles = ( "cernlib-2003.tar.gz", "fastNLO-rev${frev}.tar.gz", "fastrun.pl",
                "nlojet++-2.0.1.tar.gz", "nlojet++-2.0.1-fix.tar.gz" );
if ( $pdf eq "LHAPDF" ) { push @infiles, "lhapdf-4.2.tar.gz"; }

my %sublines;
if ( $batch eq "LCG" ) {
    $sublines{VirtualOrganisation} = "\"cms\"";
    $sublines{RetryCount} = "0";
    $sublines{InputSandbox} = "{";
    foreach my $file ( @infiles ) {
        $sublines{InputSandbox} .= " \"$file\",";
    }
    chop $sublines{InputSandbox};
    $sublines{InputSandbox} .= " }";
    $sublines{Executable} = "\"fastrun.pl\"";
} elsif ( $batch eq "PBS" ) {
    $sublines{-j} = "oe";
#    $sublines{-W} = "stagein=";
#    foreach my $infile ( @infiles ) {
#       $sublines{-W} .= "${infile}\@${host}:${infile},";
#    }
}

for ( my $i = 0; $i < $njobs; $i++) {
    my $table = "${scen}${ref}/${scen}${ref}-hhc-$runmode{$order}[0]-2jet";
    my @outfiles = ( "${table}_${jobnr}", "${table}.raw_${jobnr}",
                     "${scen}${ref}_${order}.log_${jobnr}" );
    if ( $batch eq "LCG" ) {
        $sublines{OutputSandbox} = "{";
        foreach my $file ( @outfiles ) {
            $sublines{OutputSandbox} .= " \"$file\",";
        }
        chop $sublines{OutputSandbox};
        $sublines{OutputSandbox} .= " }";
        $sublines{Arguments} = "\"-f $frev -j $jobnr -o $order -p $pdf ";
        if ( $idir ne "." ) { $sublines{Arguments} .= "-d $idir ";}
        if ( $ref ) { $sublines{Arguments} .= "-r ";}
        $sublines{Arguments} .= "$scen\"";
        $sublines{StdOutput} = "\"${scen}${ref}_${order}.log_${jobnr}\"";
        $sublines{StdError} = "\"${scen}${ref}_${order}.log_${jobnr}\"";
        my $jdlfil = "${scen}${ref}_${order}.jdl_${jobnr}";
        open(NEWFILE,">$jdlfil") or
            die "fastsub.pl: Error! Could not open $jdlfil!";
        foreach my $key ( sort keys %sublines ) {
            print NEWFILE "$key = $sublines{$key};\n";
        }
        close NEWFILE;
        my $cmd = "edg-job-submit -r $cequ $jdlfil >& $jdlfil.http";
        print "fastsub.pl: Running command $cmd\n";
        system("$cmd");
    } elsif ( $batch eq "PBS" ) {
        my $sdir = $tdir."${scen}${ref}_${jobnr}";
        $sublines{-o} = "${sdir}/${scen}${ref}_${order}.log_${jobnr}";
        $sublines{-N} = "jobnr_${jobnr}";
        my $cshfil = "${scen}${ref}_${order}.csh_${jobnr}";
        open(NEWFILE,">$cshfil") or
            die "fastsub.pl: Error! Could not open $cshfil!";
        print NEWFILE "#!/bin/csh -f\n";
        foreach my $key ( sort keys %sublines ) {
            print NEWFILE "#PBS $key $sublines{$key}\n";
        }
        print NEWFILE "mkdir $sdir\n";
        print NEWFILE "cd $sdir\n";
        foreach my $infile ( @infiles ) {
            unless ( $mode == 3 && $infile =~ m/tar/ ) {
                print NEWFILE "cp -p \$PBS_O_WORKDIR/${infile} .\n";
            }
        }
        print NEWFILE "echo Actual directory: \`pwd\`:\n";
        print NEWFILE "ls -la\n";
        print NEWFILE "echo Working directory: \$PBS_O_WORKDIR\n";
        print NEWFILE "ls -la \$PBS_O_WORKDIR\n";
#       print NEWFILE "echo cd to working dir ...\n";
        my $args = "-v -b $batch -t $sdir -f $frev -j $jobnr -m $mode -o $order -p $pdf ";
        if ( $idir ne "." ) { $args .= "-d $idir ";}
        if ( $nmax ) { $args .= "-e $nmax ";}
        if ( $ref ) { $args .= "-r ";}
        if ( $scen ) { $args .= "$scen";}
        print NEWFILE "echo \"Running command ./fastrun.pl $args\\n\"\n";
        print NEWFILE "./fastrun.pl $args\n";
        print NEWFILE "echo Done!\\n\n";
        close NEWFILE;
        system("chmod 755 $cshfil");
        my $cmd = "qsub -q $cequ $cshfil";
        print "fastsub.pl: Running command $cmd\n";
        system("$cmd");
    }
    $jobnr++;
}

exit 0;
