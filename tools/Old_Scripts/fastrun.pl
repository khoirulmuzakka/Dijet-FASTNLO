#!/usr/bin/env perl
#
# fastNLO run script
# Version:
#
# created by K. Rabbertz: 29.01.2006
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
# Tee STDOUT and STDERR into file
my $gjobnr = "";
if ( defined $ENV{MY_JOBID} ) {
    $gjobnr = $ENV{MY_JOBID};
#    $gjobnr++;
    $gjobnr = substr("0000$gjobnr",-4);
    open STDOUT, "| tee fastrun_${gjobnr}.log" or die
        "fastrun.pl: ERROR! Can't tee STDOUT.\n";
    open STDERR, "| tee fastrun_${gjobnr}.err" or die
        "fastrun.pl: ERROR! Can't tee STDERR.\n";
#select STDERR; $| = 1;
#select STDOUT; $| = 1;
}

my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fastrun.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################\n";
print "# fastrun.pl: Starting run of fastNLO: FASTRUN0_$date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_f, $opt_g, $opt_h, $opt_i, $opt_j, $opt_m,
      $opt_n, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_t, $opt_v, $opt_w, $opt_x ) =
    ( "LOCAL", "", "0", "2.3.1pre-1804", "guc", "", ".", "0001", "0",
      "2jet", "LO", "LHAPDF", "none", "", ".", "", "2.1", "", "" );
getopts('b:de:f:g:hi:j:m:n:o:p:q:rs:t:v:wx:') or die "fastrun.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastrun.pl\n";
    print "Usage: fastrun.pl [switches/options] ([scenario])\n";
    print "  -b batch        Batch system used: LOCAL (def.),\n";
    print "                  GC (grid-control), GRID or PBS\n";
    print "  -d debug        Switch debug/verbose mode on\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -f rev|vers     fastNLO revision to use (v1.4, e.g.=542) or\n";
    print "                  fastNLO version to use (v2, def.=2.3.1pre-1804)\n";
    print "  -g prot         Grid storage protocol to use: guc (def.), srm\n";
    print "  -h              Print this text\n";
    print "  -i dir          Installation directory (def.=.)\n";
    print "  -j jobnr        Job number to attach (def.=0001)\n";
    print "  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only\n";
    print "  -n njet         NLOJet++ jet mode: 2jet (def.) for 2+, 3jet for 3+ final states\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -p pdf          CTEQ parton densities or LHAPDF (def.)\n";
    print "  -q rootversion  ROOT version to install (def.=none)\n";
    print "  -r              Reference calculation incl. pdf access\n";
    print "  -s dir          Archive source directory, abs. path (def.=`pwd`)\n";
    print "  -t dir          Output target directory: ".
        "(def.= {scen}{ref}_{jobnr} with\n                  ".
        "ref. to working directory in fastNLO installation)\n";
    print "  -v #            Choose between fastNLO version 1a, 1b, 2a, 2 and 2.1 (def.=2.1)\n";
    print "  -w              Warm-up run to determine x limits (def.=F; version 2 only)\n";
    print "  -x seed         Seed for random number generator (def.=F; version 2 only)\n";
    print "\n";
    print "Examples:\n";
    print "1) Install only (to install with LHAPDF use option -p):\n";
    print "   ./fastrun.pl [-i .|installdir] [-f rev|vers] -m 1 [-p CTEQ|LHAPDF] [-s .|sdir] [-v 1a|1b|2a|2|2.1]\n\n";
    print "2) Make only scenario (to make scenario for reference mode use option -r):\n";
    print "   ./fastrun.pl [-i .|installdir] [-f rev|vers] -m 2 [-p CTEQ|LHAPDF] [-r] [-v 1a|1b|2a|2|2.1] scenarioname\n\n";
    print "3) Run only (to run scenario in reference mode use option -r):\n";
    print "   ./fastrun.pl [-b LOCAL|GRID|batch] [-i .|installdir] [-e max-events] [-f rev|vers] -m 3 [-p CTEQ|LHAPDF] [-r] [-t ./{scen}{ref}_{jobnr}|tdir] [-v 1a|1b|2a|2|2.1] scenarioname\n\n";
    exit;
}

unless ( $opt_b eq "LOCAL" || $opt_b eq "GC" || $opt_b eq "GRID" || $opt_b eq "PBS" ) {
    die "fastrun.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( ($opt_f =~ m/\d+/ && $opt_f !~ m/\D+/) || $opt_f =~ m/\d\.\d\.\d/ ) {
    die "fastrun.pl: Error! Illegal fastNLO revision number: $opt_f, aborted.\n";
}
unless ( $opt_g eq "guc" || $opt_g eq "srm" ) {
    die "fastrun.pl: Error! No such grid storage protocol: $opt_g, aborted.\n";
}
unless ( -d $opt_i ) {
    die "fastrun.pl: Error! No such directory: $opt_i, aborted.\n";
}
unless ( $opt_j =~ m/\d{4}/ && $opt_j !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal job number (nnnn): $opt_j, aborted.\n";
}
unless ( $opt_m == 0 || $opt_m == 1 || $opt_m == 2 || $opt_m == 3 ) {
    die "fastrun.pl: Error! Illegal job mode: $opt_m, aborted.\n";
}
unless ( $opt_n eq "2jet" || $opt_n eq "3jet" ) {
    die "fastrun.pl: Error! Illegal njet mode: $opt_n, aborted.\n";
}
unless ( $opt_o eq "LO" || $opt_o eq "NLO" ) {
    die "fastrun.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( $opt_p eq "CTEQ" || $opt_p eq "LHAPDF" ) {
    die "fastrun.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( -d $opt_s || -l $opt_s ) {
    die "fastrun.pl: Error! Archive source directory $opt_s does not exist, aborted.\n";
}
unless ( $opt_v eq "1a" || $opt_v eq "1b" || $opt_v eq "2a" || $opt_v eq "2" || $opt_v eq "2.1" ) {
    die "fastrun.pl: Error! fastNLO version $opt_v does not exist, aborted.\n";
}

my $batch = $opt_b;
my $idir  = $opt_i;
#if ( $idir eq "." ) {
#    $idir = `pwd`;
#    chomp $idir;
#}
my $nmax  = $opt_e;
my $frev  = $opt_f;
my $prot  = $opt_g;
my $jobnr = $opt_j;
my $mode  = $opt_m;
my $njet  = $opt_n;
my $order = $opt_o;
my $pdf   = $opt_p;
my $rv    = $opt_q;
my $ref   = "";
if ( $opt_r ) { $ref = "ref";}
my $sdir  = $opt_s;
my $verb  = "";
my $vers  = $opt_v;
my $wrm   = "";
if ( $opt_w ) { $wrm = "wrm";}
if ( $opt_r && $opt_w ) {
    die "fastrun.pl: Error! Reference and warm-up mode are mutually exclusive!\n";
}
my $seed  = "";
if ( $opt_x ) { $seed = $opt_x ;}
print "fastrun.pl: Directory for/of installation is $idir\n";
print "fastrun.pl: Using fastNLO revision $frev\n";
print "fastrun.pl: Job mode is $mode\n";
unless ( $mode == 1 || $mode == 2 ) {
    print "fastrun.pl: Running on batch system $batch\n";
    print "fastrun.pl: Using grid storage protocol $prot\n";
    print "fastrun.pl: Maximal event number: $nmax\n";
    print "fastrun.pl: Attaching job number $jobnr\n";
    print "fastrun.pl: Running in order $order\n";
    if ( $ref eq "ref" ) {
        print "fastrun.pl: Running with pdf $pdf\n";
    }
    if ( $wrm eq "wrm" ) {
        print "fastrun.pl: Running in warm-up mode\n";
    }
}
if ( $mode != 1 && $ref ) {
    print "fastrun.pl: Making/running in reference mode\n";
}
if ( $mode == 0 || $mode == 1 ) {
    print "fastrun.pl: Looking for sources in $sdir\n";
}
if ( $opt_d ) {
    $verb = 1;
    print "fastrun.pl: Debug/verbose mode is active.\n";
}

#
# Parse arguments
#
my $scen = "def";
unless ( $mode == 1 || @ARGV == 1 ) {
    die "fastrun.pl: Error! Need one scenario name if mode != 1!\n";
}
if ( @ARGV > 0 ) {
    $scen   = shift;
}
my $tdir   = "${scen}${ref}_${jobnr}";
if ( $opt_t ) {
    $tdir  = $opt_t;
}
unless ( $mode == 1 || $mode == 2 ) {
    print "fastrun.pl: Output target directory $tdir\n";
}

#
# Set signal traps to store output table at end of allocated batch time
#
# Optimize batch running time, does not work reliably!
#if ( $batch eq "GRID" ) {
#    $SIG{INT}  = \&grid_storage;
#    $SIG{TERM} = \&grid_storage;
#}

#
# Initialization
#
# Switches for scenario cc file
my $refsw0 = " iref = 0";
my $refsw1 = " iref = 1";
# Makefile switches
my $maksw0 = "-o ${scen} ";
my $maksw1 = "-o ${scen}ref ";

# Run mode settings: order, ordername, # events
my %runmode;
$runmode{LO}[0]  = "born";
#$runmode{LO}[1]  = "100000000";
$runmode{LO}[1]  = "10000000";
if ( $nmax > 0 && $nmax < $runmode{LO}[1] ) {
    $runmode{LO}[1]  = "$nmax";
}
$runmode{NLO}[0] = "nlo";
$runmode{NLO}[1] = "10000000";
if ( $nmax > 0 && $nmax < $runmode{NLO}[1] ) {
    $runmode{NLO}[1]  = "$nmax";
}

# NLOJET++ table name
my $tabext;
if ( $vers eq "1a" || $vers eq "1b" ) {
    $tabext = "raw";
} else {
    $tabext = "tab";
}
my $tabnam = "${scen}${ref}${wrm}_${jobnr}-hhc-$runmode{$order}[0]-${njet}.${tabext}";

# Directories
my $pwdir = getcwd();
unless ( -d "$idir/src" ) {
    my $ret = system ("mkdir -p $idir/src");
    if ( $ret ) {die "fastrun.pl: Couldn't create unpacking ".
                     "directory $idir/src: $ret, aborted!\n";}
}
chdir "$idir";
my $aidir = getcwd();
chomp $aidir;
chdir "$pwdir";



# Define install hash
my %install;
# 0: Archive names
# 1: Unpacking directory of the archives
# 2: Final directories
# 3: Link shortcut
#
# Store used gcc version to add this to final directories
my $gccv1 = "3.3.6";
my $gccvers = `gcc -dumpversion`;
if ($ENV{CC}) {
    $gccvers = `$ENV{CC} -dumpversion`;
}

if ( $? ) {
    if ( $vers eq "1a" || $vers eq "1b" ) {
        print "fastrun.pl: Warning! System gcc command failed. ".
            " Will use version $gccv1 anyway.\n";
    } else {
        if ( $mode == 0 || $mode == 1 || $mode == 2 ) {
            die "fastrun.pl: Error! System gcc command failed, ".
                " aborted.\n";
        } else {
            print "fastrun.pl: Warning! System gcc command failed. ".
                " Will try to run anyway.\n";
        }
    }
}
if ( $vers eq "1a" || $vers eq "1b" ) {
    $gccvers = $gccv1;
} else {
    chomp $gccvers;
}
my $gccapp = "gcc$gccvers";
$gccapp =~ s/\.//g;
print "fastrun.pl: Using gcc compiler version $gccvers\n";

#$install{mcfm}[0]       = "mcfm-5.3";
#$install{mcfm}[1]       = "MCFM-5.3";
#$install{mcfmfix}[0]    = "mcfm-5.3-fix";
#$install{mcfmfix}[1]    = "MCFM-5.3";
if ( $vers eq "1a" || $vers eq "1b" ) {
    $install{fastNLO}[0]    = "fastNLO-rev${frev}";
    $install{fastNLO}[1]    = "fastNLO-rev${frev}";
    $install{gcccore}[0]    = "gcc-core-$gccvers";
    $install{gcccore}[1]    = "gcc-$gccvers";
    $install{gccgpp}[0]     = "gcc-g\+\+-$gccvers";
    $install{gccgpp}[1]     = "gcc-$gccvers";
    $install{gccg77}[0]     = "gcc-g77-$gccvers";
    $install{gccg77}[1]     = "gcc-$gccvers";
    $install{cernlib}[0]    = "cernlib-2003";
    $install{cernlib}[1]    = "cernlib-2003";
    $install{fastjet}[0]    = "fastjet-2.4.4";
    $install{fastjet}[1]    = "fastjet-2.4.4";
# Versions >= 5.4.0 don't work with gcc 3.3.6
    $install{lhapdf}[0]     = "lhapdf-5.3.1";
    $install{lhapdf}[1]     = "lhapdf-5.3.1";
    $install{nlojet}[0]     = "nlojet++-2.0.1";
    $install{nlojet}[1]     = "nlojet++-2.0.1";
    $install{nlojetfix}[0]  = "nlojet++-2.0.1-fix";
    $install{nlojetfix}[1]  = "nlojet++-2.0.1";
} elsif ( $vers eq "2a" || $vers eq "2" ) {
    if ( $gccvers lt "4.0.0" ) {
# 2006 seems to work on slc4
        $install{cernlib}[0]    = "cernlib-2006";
        $install{cernlib}[1]    = "cernlib-2006";
    } else {
# 2003 seems to work on SuSE 11.0 & Mandriva 10
        $install{cernlib}[0]    = "cernlib-2003";
        $install{cernlib}[1]    = "cernlib-2003";
    }
    $install{root}[0]       = "root-$rv";
    $install{root}[1]       = "root";
    $install{fastjet}[0]    = "fastjet-3.0.3";
    $install{fastjet}[1]    = "fastjet-3.0.3";
    $install{lhapdf}[0]     = "lhapdf-5.9.1";
    $install{lhapdf}[1]     = "lhapdf-5.9.1";
    $install{lhapdffix}[0]  = "";
    $install{lhapdffix}[1]  = "";
    $install{lhpdf}[0]      = "lhpdf-1.0.0";
    $install{lhpdf}[1]      = "lhpdf-1.0.0";
    $install{lhpdffix}[0]   = "lhpdf-1.0.0-fix";
    $install{lhpdffix}[1]   = "lhpdf-1.0.0";
    $install{fastNLO}[0]    = "fastnlo-${frev}";
    $install{fastNLO}[1]    = "fastnlo-${frev}";
    if ( $vers eq "2a" ) {
        $install{nlojet}[0]     = "nlojet++-4.0.1";
        $install{nlojet}[1]     = "nlojet++-4.0.1";
        $install{nlojetfix}[0]  = "nlojet++-4.0.1-fix";
        $install{nlojetfix}[1]  = "nlojet++-4.0.1";
        $install{znpdf}[0]      = "";
        $install{znpdf}[1]      = "";
        $install{znpdffix}[0]   = "";
        $install{znpdffix}[1]   = "";
        $install{fastNLOfix}[0] = "fastNLO-rev${frev}-fix";
        $install{fastNLOfix}[1] = "fastNLO-rev${frev}-fix";
    } else {
#        $install{nlojet}[0]     = "nlojet++-4.1.2";
#        $install{nlojet}[1]     = "nlojet++-4.1.2";
#        $install{nlojetfix}[0]  = "nlojet++-4.1.2-fix";
#        $install{nlojetfix}[1]  = "nlojet++-4.1.2";
        $install{nlojet}[0]     = "nlojet++-4.1.3";
        $install{nlojet}[1]     = "nlojet++-4.1.3";
        $install{nlojetfix}[0]  = "nlojet++-4.1.3-fix";
        $install{nlojetfix}[1]  = "nlojet++-4.1.3";
        $install{znpdf}[0]      = "cteq-pdf-1.0.4";
        $install{znpdf}[1]      = "cteq-pdf-1.0.4";
        $install{znpdffix}[0]   = "";
        $install{znpdffix}[1]   = "";
        $install{fastNLOfix}[0] = "";
        $install{fastNLOfix}[1] = "";
    }
} else {
    $install{fastjet}[0]       = "fastjet-3.0.6";
    $install{fastjet}[1]       = "fastjet-3.0.6";
    $install{lhapdf}[0]        = "lhapdf-5.9.1";
    $install{lhapdf}[1]        = "lhapdf-5.9.1";
    $install{lhapdffix}[0]     = "";
    $install{lhapdffix}[1]     = "";
    $install{nlojet}[0]        = "nlojet++-4.1.3";
    $install{nlojet}[1]        = "nlojet++-4.1.3";
    $install{nlojetfix}[0]     = "nlojet++-4.1.3-fix";
    $install{nlojetfix}[1]     = "nlojet++-4.1.3";
    $install{fastNLOtk}[0]     = "fastnlo_toolkit-${frev}";
    $install{fastNLOtk}[1]     = "fastnlo_toolkit-${frev}";
    $install{fastNLOnlojet}[0] = "fastnlo_interface_nlojet-${frev}";
    $install{fastNLOnlojet}[1] = "fastnlo_interface_nlojet-${frev}";
}

foreach my $comp ( keys %install ) {
    if ( $comp =~ m/gcc/ ) {
        $install{$comp}[2] = $install{$comp}[1];
        $install{$comp}[3] = "gcc";
#    } elsif ( $comp =~ m/cern/ ) {
#        $install{$comp}[2] = "cernlib-2006_slc4_ia32_gcc4";
#        $install{$comp}[3] = "cernlib-2006_slc4_ia32_gcc4-v2";
    } elsif ( $comp =~ m/fix/ ) {
        $install{$comp}[2] = "";
        $install{$comp}[3] = "";
    } else {
        $install{$comp}[2] = $install{$comp}[0]."-$gccapp";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            $install{$comp}[3] = $install{$comp}[0]."-v1";
        } else {
            $install{$comp}[3] = $install{$comp}[0]."-v2";
        }
    }
    $install{$comp}[0] .=".tar.gz";
    if ( $verb ) {
        print "fastrun.pl: Archive name: $install{$comp}[0]\n";
        print "fastrun.pl: Unpack dir  : $install{$comp}[1]\n";
        print "fastrun.pl: Install dir : $install{$comp}[2]\n";
        print "fastrun.pl: Final link  : $install{$comp}[3]\n";

    }
    unless ( $mode == 2 || $mode == 3 ||
             -d "$idir/$install{$comp}[2]" ||
             -f "$sdir/$install{$comp}[0]" ||
             ( $comp eq "lhapdf" && $pdf eq "CTEQ" ) ||
             ( $comp eq "root" && $rv eq "none" ) ) {
        die "fastrun.pl: Archive $install{$comp}[0] does not exist, ".
            "aborted!\n";
    }
}

#
# Print system info
#
if ( $verb ) {
    print "\n######################################################\n";
    print "fastrun.pl: System information for debugging purposes:\n";
    print "######################################################\n";
    my $host = `hostname`;
    if ( $? ) {
        print "fastrun.pl: Info: \"hostname\" command failed.\n\n";
    } else {
        chomp $host;
        print "fastrun.pl: The system's hostname is (hostname):\n$host\n\n";
    }
    my $osvers = `uname -a`;
    if ( $? ) {
        print "fastrun.pl: Info: \"uname -a\" command failed.\n\n";
    } else {
        print "fastrun.pl: Your operating system is (uname -a):\n$osvers\n";
    }
    my $procvers = `cat /proc/version`;
    if ( $? ) {
        print "fastrun.pl: Info: \"cat /proc/version\" command failed.\n\n";
    } else {
        print "fastrun.pl: Your linux version is (/proc/version):\n$procvers\n";
    }
    my $freemem = `free`;
    if ( $? ) {
        print "fastrun.pl: Info: \"free\" command failed.\n\n";
    } else {
        print "fastrun.pl: The available memory is:\n$freemem\n";
    }
    my $cpumod = `cat /proc/cpuinfo | grep \"model name\"`;
    if ( $? ) {
        print "fastrun.pl: Info: \"cat /proc/cpuinfo\" command failed.\n\n";
    } else {
        my $cpufrq = `cat /proc/cpuinfo | grep \"cpu MHz\"`;
        print "fastrun.pl: The processor type is:\n${cpumod}at\n${cpufrq}\n";
    }
    my $freedisk = `df -h`;
    if ( $? ) {
        print "fastrun.pl: Info: \"df -h\" command failed.\n\n";
    } else {
        print "fastrun.pl: The available disk space is:\n$freedisk\n";
    }
    my $freenode = `df -hi`;
    if ( $? ) {
        print "fastrun.pl: Info: \"df -hi\" command failed.\n\n";
    } else {
        print "fastrun.pl: The available inode space is:\n$freenode\n";
    }
    my $cwd = getcwd();
    if ( $? ) {
        print "fastrun.pl: Info: \"getcwd()\" command failed.\n\n";
    } else {
        print "fastrun.pl: The current working directory is:\n$cwd\n\n";
        print "fastrun.pl: The current working directorys content is:\n";
        system("ls -la");
    }
#print "fastrun.pl: Installation environment:\n";
#system("printenv");
    print "######################################################\n\n";
}



#
# Installation
#
if ( $mode == 0 || $mode == 1 ) {
    print "\nfastrun.pl: Setting fastNLO and gcc environment as ".
        "required by this installation of fastNLO version $vers!\n";
    print "fastrun.pl: Check fastNLO.(c)sh for reference and ".
        "later \"source\"ing.\n";
    open (FILE,">> fastNLO.csh");
    print FILE "#\n";
    print FILE "# fastrun.pl: Setting fastNLO and gcc environment as ".
        "required by this installation of fastNLO version $vers!\n";
    print FILE "#\n";



#
# 0) Install gcc
#
    if ( $vers eq "1a" || $vers eq "1b" ) {
# Standard environment to use system compiler still
        print FILE "# Add GCC 3.6.6 environment\n";
        unless ( -e "$idir/$install{gcccore}[2]" ) {
            $date = `date +%d%m%Y_%H%M%S`;
            chomp $date;
            print "\nfastrun.pl: Unpacking gcc-core sources in $install{gcccore}[0]: $date\n";
            my $ret = system("tar xz -C $idir -f $sdir/$install{gcccore}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gcccore}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            if ( -l "$idir/$install{gcccore}[3]" ) {system("rm -f $idir/$install{gcccore}[3]");}
            system("ln -s $install{gcccore}[2] $idir/$install{gcccore}[3]");
            chdir "$idir/$install{gcccore}[3]";
            my $cmd = "./configure ".
                "--prefix=$aidir/gcc ".
                "--bindir=$aidir/bin ".
                "--libdir=$aidir/lib ".
                "--with-gxx-include-dir=$aidir/include";
            print "fastrun.pl: Configuring gcc for first compile (system compiler): $cmd ...\n";
# Bugfix gcc
            $ENV{SHELL} = "/bin/sh";
            $ret = system("$cmd");
            if ( $ret ) {die "fastrun.pl: 1st configure step of gcc failed: $ret, aborted!\n";}
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: 1st make step of gcc failed: $ret, aborted!\n";}
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: 1st install step of gcc failed: $ret, aborted!\n";}
#        my @acmd = split(" ",$cmd);
#        system(@acmd);
            $ENV{PATH} = "$aidir/bin:$ENV{PATH}";
            if ( $ENV{LD_LIBRARY_PATH} ) {
                $ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64:".
                    "$ENV{LD_LIBRARY_PATH}";
            } else {
                $ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64";
            }
            $ENV{GCC_EXEC_PREFIX} = "$aidir/lib/gcc-lib/";
            chdir "..";
            print "fastrun.pl: Unpacking gcc-g++ and gcc-g77 sources in $install{gcccore}[1]: $date\n";
            $ret = system("tar xz -C $idir -f $sdir/$install{gccgpp}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gccgpp}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            $ret = system("tar xz -C $idir -f $sdir/$install{gccg77}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gccg77}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            chdir "$idir/gcc";
            system("make clean");
            print "fastrun.pl: Configuring gcc for recompile (incl. gcc): $cmd ...\n";
# In case a different "system's" compiler was chosen via CC, CXX and FC environment variables
            if ( $ENV{CC} ) {
                print "fastrun.pl: Changing system's CC compiler $ENV{CC} to $aidir/bin/gcc\n";
            }
            $ENV{CC}  = "$aidir/bin/gcc";
            if ( $ENV{CXX} ) {
                print "fastrun.pl: Changing system's CXX compiler $ENV{CXX} to $aidir/bin/g++\n";
            }
            $ENV{CXX}  = "$aidir/bin/g++";
            if ( $ENV{FC} ) {
                print "fastrun.pl: Changing system's FC compiler $ENV{FC} to $aidir/bin/g77\n";
            }
            $ENV{FC}  = "$aidir/bin/g77";
            $ret = system("$cmd");
            if ( $ret ) {die "fastrun.pl: 2nd configure step of gcc failed: $ret, aborted!\n";}
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: 2nd make step of gcc failed: $ret, aborted\n";}
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: 2nd install step of gcc failed: $ret, aborted\n";}
            chdir "$aidir";
        }
# Switch to proper gcc-3.3.6 compiler if not already done
        $ENV{SHELL} = "/bin/sh";
        $ENV{PATH} = "$aidir/bin:$ENV{PATH}";
        if ( $ENV{LD_LIBRARY_PATH} ) {
            $ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64:".
                "$ENV{LD_LIBRARY_PATH}";
        } else {
            $ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64";
        }
        print FILE "setenv GCC_EXEC_PREFIX $aidir/lib/gcc-lib/\n";
        print FILE "setenv CXXFLAGS \"-O3 -I .\"\n";
        print FILE "#\n";
        $ENV{GCC_EXEC_PREFIX} = "$aidir/lib/gcc-lib/";
        $ENV{CXXFLAGS} = "-O3 -I .";
    }



#
# 1) Install CERN libraries and ROOT (V2 only)
#
    unless ( (($vers eq "1a" || $vers eq "1b") && $ENV{CERNLIB}) ||
             (($vers eq "2") && $ENV{CERN_ROOT}) ||
             ($vers eq 2.1)
        ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        chdir "$aidir/src";
        print "\nfastrun.pl: Unpacking CERN libraries in $install{cernlib}[1]: $date\n";
        my $ret = system("tar xz -f $sdir/$install{cernlib}[0]");
        if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{cernlib}[0] ".
                         "in $aidir/src failed: $ret, aborted!\n";}
        $ret = system("mv $install{cernlib}[1] $install{cernlib}[2]");
        if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                         "$install{cernlib}[1] to ".
                         "$install{cernlib}[2]: $ret, aborted!\n";}
        if ( $vers eq "1a" || $vers eq "1b" ) {
            unless ( -d "$aidir/$install{cernlib}[2]" ) {
                my $ret = system("mv $install{cernlib}[2] $aidir");
                if ( $ret ) {die "fastrun.pl: Couldn't move cernlib ".
                                 "directory to $aidir: $ret, aborted!\n";}
            }
#            if ( -l "$aidir/$install{cernlib}[3]" ) {system("rm -f $aidir/$install{cernlib}[3]");}
#            chdir $aidir;
#            system("ln -s $install{cernlib}[2] $idir/$install{cernlib}[3]");
        } else {
            unless ( -d "$aidir/lib" ) {
                my $ret = system ("mkdir -p $aidir/lib");
                if ( $ret ) {die "fastrun.pl: Couldn't create lib ".
                                 "directory $aidir/lib: $ret, aborted!\n";}
            }
            $ret = system("cp -p $install{cernlib}[2]/\* $aidir/lib");
            if ( $ret ) {die "fastrun.pl: Couldn't copy CERNlibs into ".
                             "directory $aidir/lib: $ret, aborted!\n";}
        }
        print FILE "# Add CERNLIB environment\n";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            print FILE "setenv CERNLIB $aidir/$install{cernlib}[2]\n";
            $ENV{CERNLIB}  = "$aidir/$install{cernlib}[2]";
        } else {
            print FILE "setenv CERN_ROOT $aidir\n";
            print FILE "setenv CERNLIBPATH $aidir/lib\n";
            print FILE "setenv CERNLIBS \"-L\$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib\"\n";
            $ENV{CERN_ROOT}   = "$aidir";
            $ENV{CERNLIBPATH} = "$aidir/lib";
            $ENV{CERNLIBS}    = "-L$ENV{CERN_ROOT}/lib -lmathlib -lkernlib -lpacklib";
        }
        print FILE "#\n";
    }
    chdir "$pwdir";
# ROOT
    if ( ($vers eq "2") ) {
#        unless ( -e "$aidir/src/$install{root}[2]" ) {
        unless ( $ENV{ROOTSYS} || $rv eq "none") {
            $date = `date +%d%m%Y_%H%M%S`;
            chomp $date;
            chdir "$aidir/src";
            print "\nfastrun.pl: Installing ROOT from $install{root}[0]: $date\n";
            if ($verb) {print "fastrun.pl: Unpacking directory is $aidir/src\n";}
            print "\nfastrun.pl: Unpacking $install{root}[0] ...\n";
            my $ret = system("tar xz -f $sdir/$install{root}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive ".
                             "$sdir/$install{root}[0] ".
                             "into $aidir/src failed: $ret, aborted!\n";}
            $ret = system("mv $install{root}[1] $install{root}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{root}[1] to ".
                             "$install{root}[2]: $ret, aborted!\n";}
#            if ( -l "$idir/$install{root}[3]" ) {system("rm -f $idir/$install{root}[3]");}
#            system("ln -s $install{root}[2] $idir/$install{root}[3]");
            chdir "$install{root}[2]";
            print "\nfastrun.pl: Configuring ROOT ...\n";
            my $cmd = "./configure --enable-roofit --prefix=$aidir --etcdir=$aidir/etc --with-cxx=\"gcc -Wall\"\n";
            if ( $verb ) {print "$cmd"};
            $ret = system("./configure --enable-roofit --prefix=$aidir --etcdir=$aidir/etc --with-cxx=\"gcc -Wall\"");
            if ( $ret ) {die "fastrun.pl: Error $ret in ROOT configure step, aborted!\n";}
            print "\nfastrun.pl: Making ROOT ...\n";
# No multithreaded make for ROOT :-(
#            $ret = system("make -j2");
            $ret = system("make");
            if ( $ret ) {die "fastrun.pl: Error $ret in ROOT make step, aborted!\n";}
            print "\nfastrun.pl: Make install for ROOT ...\n";
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: Error $ret in ROOT make install step, aborted!\n";}
            print FILE "# Add ROOT environment\n";
            print FILE "setenv ROOTSYS $aidir\n";
            print FILE "setenv ROOTBINPATH $aidir/bin\n";
            print FILE "setenv ROOTLIBPATH $aidir/lib/root\n";
            print FILE "setenv ROOTINCLUDEPATH $aidir/include/root\n";
            print FILE "#\n";
            $ENV{ROOTSYS}         = "$aidir";
            $ENV{ROOTBINPATH}     = "$aidir/bin";
            $ENV{ROOTLIBPATH}     = "$aidir/lib/root";
            $ENV{ROOTINCLUDEPATH} = "$aidir/include/root";
        }
    }
    chdir "$pwdir";
#exit 1;



#
# 2) Install lhapdf
#
    unless ( $ENV{LHAPDF} || $pdf eq "CTEQ") {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installing lhapdf from $install{lhapdf}[0]: $date\n";
        print FILE "# Add LHAPDF environment\n";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            print "\nfastrun.pl: Unpacking $install{lhapdf}[0] ...\n";
            my $ret = system("tar xz -C $idir -f $sdir/$install{lhapdf}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{lhapdf}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            $ret = system("mv $install{lhapdf}[1] $idir/$install{lhapdf}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{lhapdf}[1] to ".
                             "$idir/$install{lhapdf}[2]: $ret, aborted!\n";}
#            if ( -l "$idir/$install{lhapdf}[3]" ) {system("rm -f $idir/$install{lhapdf}[3]");}
#            system("ln -s $install{lhapdf}[2] $idir/$install{lhapdf}[3]");
            chdir "$idir/$install{lhapdf}[2]";
            print "\nfastrun.pl: Configuring lhapdf ...\n";
#             system("./configure --prefix=`pwd` --exec-prefix=$aidir");
# Avoid Python Murks in LHAPDF >= 5.4.0
            $ret = system("./configure --prefix=`pwd` --disable-pyext");
#             $ret = system("./configure --prefix=`pwd`");
            if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF configure step, aborted!\n";}
            print "\nfastrun.pl: Making lhapdf ...\n";
# At least until LHAPDF 5.3.1 a race condition spoils multithreaded make!
            $ret = system("make");
# Works now, maybe not
#            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make step, aborted!\n";}
            print "\nfastrun.pl: Make install for lhapdf ...\n";
            $ret = system("make install");
# In addition, there is a double lhapdf-config creation => ignore this error
#             if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make install step, aborted!\n";}
            if ( $ret ) {print "fastrun.pl: Error $ret in LHAPDF make install step ignored!\n";}
            print FILE "setenv LHAPDF  $aidir/$install{lhapdf}[2]/lib\n";
            $ENV{LHAPDF}   = "$aidir/$install{lhapdf}[2]/lib";
        } else {
            chdir "$aidir/src";
            print "\nfastrun.pl: Unpacking $install{lhapdf}[0] ...\n";
            my $ret = system("tar xz -f $sdir/$install{lhapdf}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive ".
                             "$sdir/$install{lhapdf}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            unless ( $install{lhapdffix}[0] eq ".tar.gz" ) {
                print "\nfastrun.pl: Unpacking fix for $install{lhapdf}[0] ...\n";
                $ret = system("tar xz -f $sdir/$install{lhapdffix}[0]");
                if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{lhapdffix}[0] ".
                                 "in $aidir/src failed: $ret, aborted!\n";}
            }
            $ret = system("mv $install{lhapdf}[1] $install{lhapdf}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{lhapdf}[1] to ".
                             "$install{lhapdf}[2]: $ret, aborted!\n";}
            chdir "$install{lhapdf}[2]";
            print "\nfastrun.pl: Configuring lhapdf ...\n";
#             system("./configure --prefix=`pwd` --exec-prefix=$aidir");
# Avoid Python Murks in LHAPDF >= 5.4.0
            $ret = system("./configure --prefix=$aidir --disable-pyext");
            if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF configure step, aborted!\n";}
            print "\nfastrun.pl: Making lhapdf ...\n";
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make step, aborted!\n";}
            print "\nfastrun.pl: Make install for lhapdf ...\n";
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make install step, aborted!\n";}
            print FILE "setenv LHAPDF $aidir\n";
            print FILE "setenv LHAPDFBINPATH $aidir/bin\n";
            print FILE "setenv LHAPDFLIBPATH $aidir/lib\n";
            print FILE "setenv LHAPDFINCLUDEPATH $aidir/include/LHAPDF\n";
            $ENV{LHAPDF}            = "$aidir";
            $ENV{LHAPDFBINPATH}     = "$aidir/bin";
            $ENV{LHAPDFLIBPATH}     = "$aidir/lib";
            $ENV{LHAPDFINCLUDEPATH} = "$aidir/include/LHAPDF";
            $ret = `$aidir/bin/lhapdf-config --pdfsets-path`;
            chomp $ret;
            print FILE "setenv LHAPDFSETPATH $ret\n";
            $ENV{LHAPDFSETPATH} = "$ret";
        }
        print FILE "#\n";
    }
    chdir "$pwdir";
#exit 2;



#
# 3) Install fastjet
#
    unless ( $ENV{FASTJET} ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installing fastjet from $install{fastjet}[0]: $date\n";
        print FILE "# Add FASTJET environment\n";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            print "\nfastrun.pl: Unpacking $install{fastjet}[0] ...\n";
            my $ret = system("tar xz -C $idir -f $sdir/$install{fastjet}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastjet}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            $ret = system("mv $install{fastjet}[1] $idir/$install{fastjet}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastjet}[1] to ".
                             "$idir/$install{fastjet}[2]: $ret, aborted!\n";}
#            if ( -l "$idir/$install{fastjet}[3]" ) {system("rm -f $idir/$install{fastjet}[3]");}
#            system("ln -s $install{fastjet}[2] $idir/$install{fastjet}[3]");
            chdir "$idir/$install{fastjet}[2]";
            print "\nfastrun.pl: Configuring fastjet ...\n";
            $ret = system("./configure --enable-shared --enable-allplugins --prefix=`pwd` --bindir=$aidir/bin");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet configure step, aborted!\n";}
            print "\nfastrun.pl: Making fastjet ...\n";
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet make step, aborted!\n";}
            print "\nfastrun.pl: Make install for fastjet ...\n";
            $ret = system("make install");
            if ( $ret ) {print "fastrun.pl: Error $ret in fastjet make install step ignored!\n";}
            print "\nfastrun.pl: Checking fastjet ...\n";
            $ret = system("make check");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet check step, aborted!\n";}
            print FILE "setenv FASTJET $aidir/$install{fastjet}[2]\n";
            $ENV{FASTJET}  = "$aidir/$install{fastjet}[2]";
        } else {
            chdir "$aidir/src";
            print "\nfastrun.pl: Unpacking $install{fastjet}[0] ...\n";
            my $ret = system("tar xz -f $sdir/$install{fastjet}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastjet}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            $ret = system("mv $install{fastjet}[1] $install{fastjet}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastjet}[1] to ".
                             "$install{fastjet}[2]: $ret, aborted!\n";}
            chdir "$install{fastjet}[2]";
            print "\nfastrun.pl: Configuring fastjet ...\n";
            $ret = system("./configure --enable-shared --enable-allplugins --prefix=$aidir --bindir=$aidir/bin");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet configure step, aborted!\n";}
            print "\nfastrun.pl: Making fastjet ...\n";
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet make step, aborted!\n";}
            print "\nfastrun.pl: Make install for fastjet ...\n";
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet make install step, aborted!\n";}
            print "\nfastrun.pl: Checking fastjet ...\n";
            $ret = system("make check");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastjet check step, aborted!\n";}
            print FILE "setenv FASTJET $aidir\n";
            print FILE "setenv FASTJETBINPATH $aidir/bin\n";
            print FILE "setenv FASTJETLIBPATH $aidir/lib\n";
            print FILE "setenv FASTJETINCLUDEPATH $aidir/include/fastjet\n";
            $ENV{FASTJET}            = "$aidir";
            $ENV{FASTJETBINPATH}     = "$aidir/bin";
            $ENV{FASTJETLIBPATH}     = "$aidir/lib";
            $ENV{FASTJETINCLUDEPATH} = "$aidir/include/fastjet";
            $ret = `$aidir/bin/fastjet-config --libs`;
            chomp $ret;
            print FILE "setenv FASTJETLIBS \"$ret\"";
            print FILE "\n";
            $ENV{FASTJETLIBS} = "$ret";
        }
        print FILE "#\n";
    }
    chdir "$pwdir";
#exit 3;



#
# 4) Install Nlojet++
#
    unless ( $ENV{NLOJET} ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installing Nlojet++ from $install{nlojet}[0]: $date\n";
        print FILE "# Add NLOJET environment\n";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            print "\nfastrun.pl: Unpacking $install{nlojet}[0] ...\n";
            my $ret = system("tar xz -C $idir -f $sdir/$install{nlojet}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojet}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            unless ( $install{nlojetfix}[0] eq ".tar.gz" ) {
                print "\nfastrun.pl: Unpacking fix for $install{nlojet}[0] ...\n";
                $ret = system("tar xzv -C $idir -f $sdir/$install{nlojetfix}[0]");
                if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojetfix}[0] ".
                                 "in $idir failed: $ret, aborted!\n";}
            }
            $ret = system("mv $install{nlojet}[1] $idir/$install{nlojet}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{nlojet}[1] to ".
                             "$idir/$install{nlojet}[2]: $ret, aborted!\n";}
#            if ( -l "$idir/$install{nlojet}[3]" ) {system("rm -f $idir/$install{nlojet}[3]");}
#            system("ln -s $install{nlojet}[2] $idir/$install{nlojet}[3]");
            chdir "$idir/$install{nlojet}[2]";
            print "\nfastrun.pl: Configuring Nlojet++ ...\n";
#        system("./configure --prefix=`pwd` --exec-prefix=$aidir");
            $ret = system("./configure --prefix=`pwd`");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ configure step, aborted!\n";}
            print "\nfastrun.pl: Making Nlojet++ ...\n";
            $ret = system("make -j2 CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make step, aborted!\n";}
            print "\nfastrun.pl: Make install for Nlojet++ ...\n";
            $ret = system("make install CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make install step, aborted!\n";}
            print FILE "setenv NLOJET  $aidir/$install{nlojet}[2]\n";
            $ENV{NLOJET}   = "$aidir/$install{nlojet}[2]";
        } else {
            chdir "$aidir/src";
            print "\nfastrun.pl: Unpacking $install{nlojet}[0] ...\n";
            my $ret = system("tar xz -f $sdir/$install{nlojet}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojet}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            unless ( $install{nlojetfix}[0] eq ".tar.gz" ) {
                print "\nfastrun.pl: Unpacking fix for $install{nlojet}[0] ...\n";
                $ret = system("tar xz -f $sdir/$install{nlojetfix}[0]");
                if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojetfix}[0] ".
                                 "in $aidir/src failed: $ret, aborted!\n";}
            }
            $ret = system("mv $install{nlojet}[1] $install{nlojet}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{nlojet}[1] to ".
                             "$install{nlojet}[2]: $ret, aborted!\n";}
            chdir "$install{nlojet}[2]";
            print "\nfastrun.pl: Configuring Nlojet++ ...\n";
            $ret = system("./configure --prefix=$aidir");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ configure step, aborted!\n";}
            print "\nfastrun.pl: Making Nlojet++ ...\n";
            $ret = system("make -j2 CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make step, aborted!\n";}
            print "\nfastrun.pl: Make install for Nlojet++ ...\n";
            $ret = system("make install CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
            if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make install step, aborted!\n";}
            print FILE "setenv NLOJET $aidir\n";
            print FILE "setenv NLOJETBINPATH $aidir/bin\n";
            print FILE "setenv NLOJETLIBPATH $aidir/lib\n";
            print FILE "setenv NLOJETINCLUDEPATH $aidir/include/nlo++\n";
            $ENV{NLOJET}            = "$aidir";
            $ENV{NLOJETBINPATH}     = "$aidir/bin";
            $ENV{NLOJETLIBPATH}     = "$aidir/lib";
            $ENV{NLOJETINCLUDEPATH} = "$aidir/include/nlo++";
        }
        print FILE "#\n";
    }
    chdir "$pwdir";
#exit 4;



#
# 4a) Install lhpdf for Nlojet++, only for fastNLO version 2
#     (This is used within fastNLO v2 reference jobs via the pdf-cteq6.h include file)
#
    if ( ($vers eq "2") && ! $ENV{LHPDF} ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installing lhpdf from $install{lhpdf}[0]: $date\n";
        chdir "$aidir/src";
        print "\nfastrun.pl: Unpacking $install{lhpdf}[0] ...\n";
        my $ret = system("tar xz -f $sdir/$install{lhpdf}[0]");
        if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{lhpdf}[0] ".
                         "in $aidir/src failed: $ret, aborted!\n";}
        unless ( $install{lhpdffix}[0] eq ".tar.gz" ) {
            print "\nfastrun.pl: Unpacking fix for $install{lhpdf}[0] ...\n";
            $ret = system("tar xz -f $sdir/$install{lhpdffix}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{lhpdffix}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
        }
        $ret = system("mv $install{lhpdf}[1] $install{lhpdf}[2]");
        if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                         "$install{lhpdf}[1] to ".
                         "$install{lhpdf}[2]: $ret, aborted!\n";}
        chdir "$install{lhpdf}[2]";
        print "\nfastrun.pl: Configuring lhpdf ...\n";
        $ret = system("./configure --prefix=$aidir");
        if ( $ret ) {die "fastrun.pl: Error $ret in lhpdf configure step, aborted!\n";}
        print "\nfastrun.pl: Making lhpdf ...\n";
        $ret = system("make -j2");
        if ( $ret ) {die "fastrun.pl: Error $ret in lhpdf make step, aborted!\n";}
        print "\nfastrun.pl: Make install for lhpdf ...\n";
        $ret = system("make install");
        if ( $ret ) {die "fastrun.pl: Error $ret in lhpdf make install step, aborted!\n";}
        print FILE "# Add LHPDF environment\n";
        print FILE "setenv LHPDF $aidir\n";
        print FILE "setenv LHPDFLIBPATH $aidir/lib\n";
        print FILE "setenv LHPDFINCLUDEPATH $aidir/include/lhpdf\n";
        print FILE "#\n";
        $ENV{LHPDF}            = "$aidir";
        $ENV{LHPDFLIBPATH}     = "$aidir/lib";
        $ENV{LHPDFINCLUDEPATH} = "$aidir/include/lhpdf";
    }
    chdir "$pwdir";
#exit 4;



#
# 4b) Install znpdf for Nlojet++, only for fastNLO version 2
#     (NLOJet++ and fastNLO v2 compile and run without this. Drop ???)
#
    if ( ($vers eq "2") && ! $ENV{ZNPDF} ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installing znpdf from $install{znpdf}[0]: $date\n";
        chdir "$aidir/src";
        print "\nfastrun.pl: Unpacking $install{znpdf}[0] ...\n";
        my $ret = system("tar xz -f $sdir/$install{znpdf}[0]");
        if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{znpdf}[0] ".
                         "in $aidir/src failed: $ret, aborted!\n";}
        unless ( $install{znpdffix}[0] eq ".tar.gz" ) {
            print "\nfastrun.pl: Unpacking fix for $install{znpdf}[0] ...\n";
            $ret = system("tar xz -f $sdir/$install{znpdffix}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{znpdffix}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
        }
        $ret = system("mv $install{znpdf}[1] $install{znpdf}[2]");
        if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                         "$install{znpdf}[1] to ".
                         "$install{znpdf}[2]: $ret, aborted!\n";}
        chdir "$install{znpdf}[2]";
        print "\nfastrun.pl: Configuring znpdf ...\n";
        $ret = system("./configure --prefix=$aidir");
        if ( $ret ) {die "fastrun.pl: Error $ret in znpdf configure step, aborted!\n";}
        print "\nfastrun.pl: Making znpdf ...\n";
        $ret = system("make -j2");
        if ( $ret ) {die "fastrun.pl: Error $ret in znpdf make step, aborted!\n";}
        print "\nfastrun.pl: Make install for znpdf ...\n";
        $ret = system("make install");
        if ( $ret ) {die "fastrun.pl: Error $ret in znpdf make install step, aborted!\n";}
        print FILE "# Add ZNPDF environment\n";
        print FILE "setenv ZNPDF $aidir\n";
        print FILE "setenv ZNPDFLIBPATH $aidir/lib\n";
        print FILE "setenv ZNPDFINCLUDEPATH $aidir/include/cteq\n";
        print FILE "#\n";
        $ENV{ZNPDF}            = "$aidir";
        $ENV{ZNPDFLIBPATH}     = "$aidir/lib";
        $ENV{ZNPDFINCLUDEPATH} = "$aidir/include/cteq";
    }
    chdir "$pwdir";
#exit 4;



#
# 5) Install mcfm
#
# CERNLIB var already needed for mcfm
#     if ( $vers eq 3 ) {
#     print "\nfastrun.pl: Setting environment variable CERNLIB for mcfm to:\n";
#     my $cwdir = getcwd();
#     $ENV{CERNLIB} = "$ENV{CERN_ROOT}/lib";
#     print "\nfastrun.pl: $ENV{CERNLIB}\n";
#     unless ( -e "$idir/$install{mcfm}[2]" ) {
#         $date = `date +%d%m%Y_%H%M%S`;
#         chomp $date;
#         print "\nfastrun.pl: Installing mcfm from $install{mcfm}[0]: $date\n";
#         print "\nfastrun.pl: Unpacking $install{mcfm}[0] ...\n";
# # Old mcfm version 5
# #        system ("mkdir $idir/$install{mcfm}[0]");
# #        my $ret = system("tar xz -C $idir/$install{mcfm}[0] -f $sdir/$install{mcfm}[1]");
#         my $ret = system("tar xz -C $idir -f $sdir/$install{mcfm}[0]");
#         if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{mcfm}[0] ".
#                          "in $idir failed: $ret, aborted!\n";}
#         unless ( $install{mcfmfix}[0] eq ".tar.gz" ) {
#             print "\nfastrun.pl: Unpacking fix for $install{mcfm}[0] ...\n";
#             $ret = system("tar xzv -C $idir -f $sdir/$install{mcfmfix}[0]");
#             if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{mcfmfix}[0] ".
#                              "in $idir failed: $ret, aborted!\n";}
#         }
#         $ret = system("mv $install{mcfm}[1] $idir/$install{mcfm}[2]");
#         if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
#                          "$install{mcfm}[1] to ".
#                          "$idir/$install{mcfm}[2]: $ret, aborted!\n";}
# #        if ( -l "$idir/$install{mcfm}[3]" ) {system("rm -f $idir/$install{mcfm}[3]");}
# #        system("ln -s $install{mcfm}[2] $idir/$install{mcfm}[3]");
#         chdir "$idir/$install{mcfm}[2]";
#         print "\nfastrun.pl: Configuring mcfm ...\n";
#         $ret = system("./Install");
#         if ( $ret ) {die "fastrun.pl: Error $ret in mcfm install script, aborted!\n";}
#         $ret = system("make -j2");
#         if ( $ret ) {die "fastrun.pl: Error $ret in mcfm make step, aborted!\n";}
#         chdir "..";
#     }
#     }
#exit 5;



#
# 6) Install fastNLO
#
    unless ( $ENV{FASTNLO} ) {
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        if ( $vers ne "2.1" ) {
            print "\nfastrun.pl: Installing fastNLO from $install{fastNLO}[0]: $date\n";
        } else {
            print "\nfastrun.pl: Installing fastNLO from $install{fastNLOtk}[0] and $install{fastNLOnlojet}[0]: $date\n";
        }
        print FILE "# Add FASTNLO environment\n";
        if ( $vers eq "1a" || $vers eq "1b" ) {
            print "\nfastrun.pl: Unpacking $install{fastNLO}[0] ...\n";
            my $ret = system("tar xz -C $idir -f $sdir/$install{fastNLO}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLO}[0] ".
                             "in $idir failed: $ret, aborted!\n";}
            $ret = system("mv $install{fastNLO}[1] $idir/$install{fastNLO}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastNLO}[1] to ".
                             "$idir/$install{fastNLO}[2]: $ret, aborted!\n";}
#            if ( -l "$idir/$install{fastNLO}[3]" ) {system("rm -f $idir/$install{fastNLO}[3]");}
#            system("ln -s $install{fastNLO}[2] $idir/$install{fastNLO}[3]");
            print FILE "setenv FASTNLO $aidir/$install{fastNLO}[2]\n";
            $ENV{FASTNLO}  = "$aidir/$install{fastNLO}[2]";
        } elsif ( $vers eq "2a" || $vers eq "2" ) {
            chdir "$aidir/src";
            print "\nfastrun.pl: Unpacking $install{fastNLO}[0] ...\n";
            my $ret = system("tar xz -f $sdir/$install{fastNLO}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLO}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            unless ( $install{fastNLOfix}[0] eq ".tar.gz" ) {
                print "\nfastrun.pl: Unpacking fix for $install{fastNLO}[0] ...\n";
                $ret = system("tar xz -f $sdir/$install{fastNLOfix}[0]");
                if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLOfix}[0] ".
                                 "in $aidir/src failed: $ret, aborted!\n";}
            }
            $ret = system("mv $install{fastNLO}[1] $install{fastNLO}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastNLO}[1] to ".
                             "$install{fastNLO}[2]: $ret, aborted!\n";}
            chdir "$install{fastNLO}[2]";
            $ret = system("./configure --prefix=$aidir");
            if ( $ret ) {die "fastrun.pl: configure of fastNLO V2 failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: make -j2 of fastNLO V2 failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: make install of fastNLO V2 failed: ".
                             "$ret, aborted!\n";}
            print FILE "setenv FASTNLO $aidir\n";
            print FILE "setenv FASTNLOBINPATH $aidir/bin\n";
            print FILE "setenv FASTNLOLIBPATH $aidir/lib/fastnlo\n";
            print FILE "setenv FASTNLOSRCPATH $aidir/src/$install{fastNLO}[2]\n";
            $ENV{FASTNLO}            = "$aidir";
            $ENV{FASTNLOBINPATH}     = "$aidir/bin";
            $ENV{FASTNLOLIBPATH}     = "$aidir/lib/fastnlo";
            $ENV{FASTNLOSRCPATH}     = "$aidir/src/$install{fastNLO}[2]";
        } else {
            chdir "$aidir/src";
            print "\nfastrun.pl: Unpacking $install{fastNLOtk}[0] ...\n\n";
            my $ret = system("tar xz -f $sdir/$install{fastNLOtk}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLOtk}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            $ret = system("mv $install{fastNLOtk}[1] $install{fastNLOtk}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastNLOtk}[1] to ".
                             "$install{fastNLOtk}[2]: $ret, aborted!\n";}
            chdir "$install{fastNLOtk}[2]";
            $ret = system("./configure --prefix=$aidir");
            if ( $ret ) {die "fastrun.pl: configure of fastNLO_toolkit failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: make -j2 of fastNLO_toolkit failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: make install of fastNLO_toolkit failed: ".
                             "$ret, aborted!\n";}
            print "\nfastrun.pl: Unpacking $install{fastNLOnlojet}[0] ...\n\n";
            $ret = system("tar xz -f $sdir/$install{fastNLOnlojet}[0]");
            if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLOnlojet}[0] ".
                             "in $aidir/src failed: $ret, aborted!\n";}
            $ret = system("mv $install{fastNLOnlojet}[1] $install{fastNLOnlojet}[2]");
            if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
                             "$install{fastNLOnlojet}[1] to ".
                             "$install{fastNLOnlojet}[2]: $ret, aborted!\n";}
            chdir "$install{fastNLOnlojet}[2]";
            $ret = system("./configure --prefix=$aidir");
            if ( $ret ) {die "fastrun.pl: configure of fastNLO_interface_nlojet failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make -j2");
            if ( $ret ) {die "fastrun.pl: make -j2 of fastNLO_interface_nlojet failed: ".
                             "$ret, aborted!\n";}
            $ret = system("make install");
            if ( $ret ) {die "fastrun.pl: make install of fastNLO_interface_nlojet failed: ".
                             "$ret, aborted!\n";}
            print FILE "setenv FASTNLO $aidir\n";
            print FILE "setenv FASTNLOBINPATH $aidir/bin\n";
            print FILE "setenv FASTNLOLIBPATH $aidir/lib\n";
            print FILE "setenv FASTNLOSRCPATH $aidir/src/$install{fastNLOtk}[2]\n";
            $ENV{FASTNLO}            = "$aidir";
            $ENV{FASTNLOBINPATH}     = "$aidir/bin";
            $ENV{FASTNLOLIBPATH}     = "$aidir/lib";
            $ENV{FASTNLOSRCPATH}     = "$aidir/src/$install{fastNLOtk}[2]";
        }
        print FILE "#\n";
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Installation of fastNLO and all required ".
            "components finished: $date\n";
    }
}
#exit 6;



#
# Set PATH environment
#
#    my $cwdir = getcwd();
if ( $mode == 0 || $mode == 1 ) {
    print FILE "# Add to system paths PATH and LD_LIBRARY_PATH\n";}
if ( $vers eq "1a" || $vers eq "1b" ) {
    if ( $mode == 2 || $mode == 3 ) {
# For old v1 installation (1a)
        if ( $vers eq "1a" ) {
            $ENV{FASTJET}  = "$aidir/fastjet";
            $ENV{NLOJET}   = "$aidir/nlojet";
            $ENV{LHAPDF}   = "$aidir/lhapdf/lib";
            $ENV{FASTNLO}  = "$aidir/fastNLO";
# For new v1 installation (1b)
        } elsif ( $vers eq "1b" ) {
            $ENV{FASTJET}  = "$aidir/$install{fastjet}[2]";
            $ENV{NLOJET}   = "$aidir/$install{nlojet}[2]";
            $ENV{LHAPDF}   = "$aidir/$install{lhapdf}[2]/lib";
            $ENV{FASTNLO}  = "$aidir/$install{fastNLO}[2]";
        }
    }
    if ( $ENV{PATH} ) {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv PATH $aidir/bin:$ENV{FASTJET}:".
                "$ENV{NLOJET}/bin:\${PATH}\n";}
        $ENV{PATH} = "$aidir/bin:$ENV{FASTJET}:".
            "$ENV{NLOJET}/bin:$ENV{PATH}";
    } else {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv PATH $aidir/bin:$ENV{FASTJET}:".
                "$ENV{NLOJET}/bin\n";}
        $ENV{PATH} = "$aidir/bin:$ENV{FASTJET}:".
            "$ENV{NLOJET}/bin";
    }
    if ( $ENV{LD_LIBRARY_PATH} ) {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv LD_LIBRARY_PATH $aidir/lib:$aidir/lib64:".
                "$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
                "$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
                "$ENV{FASTJET}/plugins/CDFCones/.libs:".
                "$ENV{NLOJET}/lib:".
                "$ENV{LHAPDF}:\${LD_LIBRARY_PATH}\n";}
        $ENV{LD_LIBRARY_PATH} ="$aidir/lib:$aidir/lib64:".
            "$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
            "$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
            "$ENV{FASTJET}/plugins/CDFCones/.libs:".
            "$ENV{NLOJET}/lib:".
            "$ENV{LHAPDF}:$ENV{LD_LIBRARY_PATH}";
    } else {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv LD_LIBRARY_PATH $aidir/lib:$aidir/lib64:".
                "$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
                "$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
                "$ENV{FASTJET}/plugins/CDFCones/.libs:".
                "$ENV{NLOJET}/lib:".
                "$ENV{LHAPDF}\n";}
        $ENV{LD_LIBRARY_PATH} ="$aidir/lib:$aidir/lib64:".
            "$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
            "$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
            "$ENV{FASTJET}/plugins/CDFCones/.libs:".
            "$ENV{NLOJET}/lib:".
            "$ENV{LHAPDF}";
    }
} else {
    if ( $mode == 2 || $mode == 3 ) {
        $ENV{FASTJET} = "$aidir";
        $ENV{NLOJET}  = "$aidir";
        $ENV{LHAPDF}  = "$aidir";
    }
    if ( $ENV{PATH} ) {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv PATH $aidir/bin:\${PATH}\n";}
        $ENV{PATH} = "$aidir/bin:$ENV{PATH}";
    } else {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv PATH $aidir/bin\n";}
        $ENV{PATH} = "$aidir/bin";
    }
    if ( $ENV{LD_LIBRARY_PATH} ) {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv LD_LIBRARY_PATH $aidir/lib:$aidir/lib/fastnlo:".
                "$aidir/lib64:\${LD_LIBRARY_PATH}\n";}
        $ENV{LD_LIBRARY_PATH} ="$aidir/lib:$aidir/lib/fastnlo:$aidir/lib64:".
            "$ENV{LD_LIBRARY_PATH}";
    } else {
        if ( $mode == 0 || $mode == 1 ) {
            print FILE "setenv LD_LIBRARY_PATH $aidir/lib:$aidir/lib/fastnlo:".
                "$aidir/lib64\n";}
        $ENV{LD_LIBRARY_PATH} ="$aidir/lib:$aidir/lib/fastnlo:$aidir/lib64";
    }
}
if ( $mode == 0 || $mode == 1 ) {close (FILE);}
#exit 7;



#
# 7) Make fastNLO scenario
#
my $scendir;
if ( $mode == 0 || $mode == 2 ) {
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Making fastNLO scenario for version $vers: $date\n";

    if ( $vers eq "1a" || $vers eq "1b" ) {
        chdir $idir;
# Structure change in fastNLO following change in revision 212!
        if ( $frev < 212 ) {
            $scendir = "$ENV{FASTNLO}/author1c/hadron";
        } else {
            $scendir = "$ENV{FASTNLO}/trunk/v1.4/author1c/hadron";
        }
        chdir "$scendir" or die
            "fastrun.pl: Could not cd to dir $scendir!\n";

# Adapting scenario cc file for ref or noref calculation
        my $scenfil = "${scen}.cc";
        my $scenlib = "${scen}${ref}.la";
        my $deb = `pwd`;
        chomp $deb;
        if ( $verb ) {
            print "fastrun.pl: DEBUG! Current directory to run NLOJET++: $deb\n";
            print "fastrun.pl: DEBUG! ls -la:\n\n";
            system("ls -la");
            print "\nfastrun.pl: DEBUG! Check on scenario library: $scenlib\n";
        }
        if ( -f $scenlib ) {
            print "\nfastrun.pl: Scenario already compiled, make step skipped!\n";
        } else {
            print "\nfastrun.pl: Checking scenario file $scenfil in $scendir for ".
                "reference setting: $ref ...\n";
            my $oldfil = "$scendir/$scenfil";
            my $test = "";
            if ( $ref ) {
                my $cmd = "grep \"$refsw1\" $oldfil";
                if ( $verb ) {print "fastrun.pl: DEBUG! Running command $cmd\n";}
                $test = `$cmd`;
            } else {
                my $cmd = "grep \"$refsw0\" $oldfil";
                if ( $verb ) {print "fastrun.pl: DEBUG! Running command $cmd\n";}
                $test = `$cmd`;
            }
            if ( $test ) {
                print "\nfastrun.pl: Scenario file $scenfil already adapted, skipped!\n";
            } else {
                print "\nfastrun.pl: Adapting scenario file $scenfil in $scendir for ".
                    "reference setting: $ref ...\n";
                open(OLDFILE,"<$oldfil") or die "fastrun.pl: Error! Could not open $oldfil!\n";
                my $newfil = $oldfil.".new";
                open(NEWFILE,">$newfil");
                while ( my $in = <OLDFILE> ) {
                    my $out = $in;
                    if ( $ref ) {
                        $out =~ s/$refsw0/$refsw1/;
                    } else {
                        $out =~ s/$refsw1/$refsw0/;
                    }
                    print NEWFILE $out;
                }
                close(OLDFILE);
                close(NEWFILE);
                system("mv -f $oldfil ${oldfil}.old");
                system("mv -f $newfil $oldfil");
            }

# Adapting scenario Makefile for ref or noref calculation
            $scenfil = "Makefile";
            $oldfil  = "$scendir/$scenfil";
            $test = "";
            if ( $ref ) {
                my $cmd = "grep \"\\$maksw1\" $oldfil";
                if ( $verb ) {print "fastrun.pl: DEBUG! Running command $cmd\n";}
                $test = `$cmd`;
            } else {
                my $cmd = "grep \"\\$maksw0\" $oldfil";
                if ( $verb ) {print "fastrun.pl: DEBUG! Running command $cmd\n";}
                $test = `$cmd`;
            }
            if ( $test ) {
                print "\nfastrun.pl: $scenfil already adapted, skipped!\n";
            } else {
                print "\nfastrun.pl: Adapting Makefile in $scendir for ".
                    "reference setting: $ref ...\n";
                open(OLDFILE,"<$oldfil") or die "fastrun.pl: Error! Could not open $oldfil!\n";
                my $newfil = $oldfil.".new";
                open(NEWFILE,">$newfil");
                while ( my $in = <OLDFILE> ) {
                    my $out = $in;
                    if ( $ref ) {
                        $out =~ s/$maksw0/$maksw1/;
                    } else {
                        $out =~ s/$maksw1/$maksw0/;
                    }
                    print NEWFILE $out;
                }
                close(OLDFILE);
                close(NEWFILE);
                system("mv -f $oldfil ${oldfil}.old");
                system("mv -f $newfil $oldfil");
            }
            $date = `date +%d%m%Y_%H%M%S`;
            chomp $date;
            print "\nfastrun.pl: Making nlofast-add for fastNLO: $date\n";
            my $ret = system("make -j2 nlofast-add");
            if ( $ret ) {die "fastrun.pl: Error $ret in nlofast-add make step, aborted!\n";}
            print "\nfastrun.pl: Making scenario $scen of fastNLO: $date\n";
            chdir $scendir;
            $ret = system("make -j2 $scen");
            if ( $ret ) {die "fastrun.pl: Error $ret in fastNLO make step, aborted!\n";}
        }
    } else {
#        print "\nfastrun.pl: (Re-)making scenario $scen of fastNLO: $date\n";
        print "\nfastrun.pl: (Re-)making scenarios of fastNLO: $date\n";
        chdir "$ENV{FASTNLOSRCPATH}/trunk/v2.0";
#        my $ret = system("make -j2 lib${scen}.la");
        my $ret = system("make -j2");
#        if ( $ret ) {die "fastrun.pl: (re-)make -j2 lib${scen}.la ".
        if ( $ret ) {die "fastrun.pl: (re-)make -j2 ".
                         "of fastNLO V2 failed: $ret, aborted!\n";}
#        $ret = system("make install lib${scen}.la");
        $ret = system("make install");
#        if ( $ret ) {die "fastrun.pl: (re-)make install lib${scen}.la ".
        if ( $ret ) {die "fastrun.pl: (re-)make install ".
                         "of fastNLO V2 failed: $ret, aborted!\n";}
    }
}



#
# 8) Run fastNLO
#
if ( $mode == 0 || $mode == 3 ) {
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Running fastNLO version $vers scenario on batch system $batch: $date\n";

# Run without complete installation by source
# In case of non-local running need to unpack binary fastNLO distribution
    if ( $mode == 3 && $batch ne "LOCAL" ) {
# Set minimal environment with respect to current working directory
#        my $cwd = `pwd`;
#        chomp $cwd;
#        if ( $vers eq "1a" || $vers eq "1b" ) {
#            $ENV{NLOJET} = "$cwd/nlojet";
#            $ENV{FASTNLO} = "$cwd/fastNLO";
#        } else {
#            $ENV{NLOJET}  = $cwd;
#            $ENV{FASTNLO} = $cwd;
#            if (
#            $ENV{LD_LIBRARY_PATH} =
#        }

# Fetching and unpacking of fastNLO binary archive, use version incl. CTEQ PDFs for reference
        my $file = "fastNLO-bin";
#        if ( $ref ) {
        if ( $vers eq "1a" || $vers eq "1b" ) {
            $file .= "-${pdf}";
#        }
        } else {
            $file .= "-v${vers}";
        }

        $file .= ".tgz";

        if ( ! -f $file ) {
            if ( $batch eq "GC" ) {
                die "fastrun.pl: ERROR! Could not find binary tgz of fastNLO $file working dir, aborted!\nCheck your grid storage options for grid-control!\n";
            } elsif ( $batch eq "GRID" ) {
                grid_storage("FETCH","fastNLO_archives",$file,".",$file,$prot);
            } elsif ( $batch eq "PBS" ) {
                die "fastrun.pl: ERROR! Could not find binary tgz of fastNLO $file in PBS working dir, aborted!\n";
            }
        }
        if ( -f $file ) {
            system ("tar xfz $file");
        } else {
            die "fastrun.pl: ERROR! Could not find binary tgz of fastNLO $file, aborted!\n";
        }
    }

# Structure change in fastNLO following change in revision 212!
    if ( ($vers eq "1a" || $vers eq "1b") && $frev < 212 ) {
        $scendir = "$ENV{FASTNLO}/author1c/hadron";
    } elsif ( $vers eq "1a" || $vers eq "1b" ) {
        $scendir = "$ENV{FASTNLO}/trunk/v1.4/author1c/hadron";
    } else {
        $scendir = "$aidir";
    }

    system("pwd");
    system("ls -la");

    chdir "$scendir" or die
        "fastrun.pl: ERROR! Could not cd to dir $scendir!\n";

    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Running fastNLO: $date\n";
    my $cmd;
    if ( $vers eq "1a" || $vers eq "1b" ) {
        $cmd = "$ENV{NLOJET}/bin/nlojet++ -P dipole ".
            "--save-after $runmode{$order}[1] ".
            "-c$runmode{$order}[0] ".
            "-d $tdir ".
            "-n ${scen}${ref}_${jobnr} ".
            "-u ${scen}${ref}.la ";
        if ( $nmax ) { $cmd .= "--max-event $nmax"; }
    } else {
        $cmd = "$ENV{NLOJET}/bin/nlojet++ --calculate ".
            "--save-after $runmode{$order}[1] ".
            "-c$runmode{$order}[0] ".
            "-d $tdir ".
            "-n ${scen}${ref}${wrm}_${jobnr} ".
            "-u lib/lib${scen}.la ";
        if ( $nmax ) { $cmd .= "--max-event $nmax "; }
        if ( $seed ) { $cmd .= "-s $seed "; }
    }
# Do not try to maximize CPU time yet, too unstable
    if ( $batch eq "MAX" ) {
# Fork NLO calculation
        print "fastrun.pl: Forking command ((time $cmd) 2>&1)& in background\n";
        system("((time $cmd) 2>&1)&");
    } else {
# Run NLO calculation
        my $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "fastrun.pl: Starting calculation: FASTCAL0_$date\n";
        print "fastrun.pl: Running command (time $cmd) 2>&1 in foreground\n";
        my $ret = system("(time $cmd) 2>&1");
        if ( $ret ) {die "fastrun.pl: ERROR! Error $ret in fastNLO run step, aborted!\n";}
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\nfastrun.pl: Calculation finished: FASTCAL1_$date\n";
# Copy table to grid storage
        if ( $batch ne "LOCAL" ) {
            my $spath = "${scendir}/${tdir}";
            my $sfile = "${tabnam}";
# Remove ref from table path
#            my $tpath = "fastNLO_tables/${scen}${ref}";
            my $tpath = "fastNLO_tables/${scen}";
            if ( $tpath ) {
                $tpath =~ s/\/\.\//\//g;
            } else {
                $tpath = "";
            }
# Change job numbering according to grid-control
            my $tfile  = $sfile;
            $tfile =~ s/_0001//;
            if ( $batch eq "GC" ) {
                $tfile =~ s/\.${tabext}/\.${tabext}/;
            } else {
                $tfile =~ s/\.${tabext}/_${gjobnr}\.${tabext}/;
            }
# GC does not support my file rename scheme with leading zeros in job numbers
            if ( $batch eq "GC" ) {
                print "fastrun.pl: INFO: Batch mode $batch: Grid storage done by grid-control.\n";
                print "                  Copy table to current directory for storage.\n";
                my $ret = system("cp -p $spath/$sfile $pwdir/$tfile");
                if ( $ret ) {die "fastrun.pl: Couldn't copy table into ".
                                 "current directory $aidir: $ret, aborted!\n";}
            } elsif ( $batch eq "GRID" ) {
                grid_storage("TABSAV","$spath","$sfile","$tpath","$tfile",$prot);
            } else {
                my $ret = system("cp -p $spath/$sfile $pwdir/$tfile");
                if ( $ret ) {die "fastrun.pl: Couldn't copy table into ".
                                 "current directory $aidir: $ret, aborted!\n";}
            }
            unless ( $batch eq "GC" ) {
                my $date = `date +%d%m%Y_%H%M%S`;
                chomp $date;
                print "fastrun.pl: Table stored: TABSAV1_$date\n";
            }
        }
        $date = `date +%d%m%Y_%H%M%S`;
        chomp $date;
        print "\n###############################\n";
        print "# fastrun.pl: fastNLO finished: FASTRUN1_$date\n";
        print "###############################\n\n";
#        close STDERR;
#        close STDOUT;
# Copy log files to grid storage, but rename first to go with the job
# If size = 0, as should be normal for job.stderr, do NOT rename => no grid storage.
        if ( $batch ne "LOCAL" ) {
#            my $tpath = "fastNLO_tables/${scen}${ref}";
            my $tpath = "fastNLO_tables/${scen}";
            if ( $tpath ) {
                $tpath =~ s/\/\.\//\//g;
            } else {
                $tpath = "";
            }
            my $tfile = "${scen}${ref}-hhc-$runmode{$order}[0]-${njet}_${gjobnr}";
            if ( $batch eq "GC" ) {
                $tfile = "${scen}${ref}-hhc-$runmode{$order}[0]-${njet}";
            }
            chdir $pwdir or die "fastrun.pl: ERROR! Couldn't cd to $pwdir, aborted!\n";
#            if ( -z "job.stderr" ) {
#                my $ret = system("rm -f job.stderr");
#                if ( $ret ) {die "fastrun.pl: ERROR! Couldn't remove empty file job.stderr: ".
#                                 "$ret, aborted!\n";}
#            }
#            if ( -z "job.stdout" ) {
#                my $ret = system("rm -f job.stdout");
#                if ( $ret ) {die "fastrun.pl: ERROR! Couldn't remove empty file job.stdout: ".
#                                 "$ret, aborted!\n";}
#            }
# GC does not support my file rename scheme with leading zeros in job numbers
            if ( $batch eq "GC" ) {
# GC will copy renamed files to SE
                print "fastrun.pl: INFO: Batch mode $batch: Grid storage done by grid-control.\n";
                print "                  Move log files to correct file names for storage.\n\n";
                if ( -f "job.stderr" && ! -z "job.stderr" ) {
                    my $ret = system("cp -p job.stderr ${tfile}.err");
                    if ( $ret ) {die "fastrun.pl: ERROR! Couldn't copy job.stderr ".
                                     "fastrun_${gjobnr}.err: $ret, aborted!\n";}
                }
                if ( -f "job.stdout" && ! -z "job.stdout" ) {
                    $ret = system("cp -p job.stdout ${tfile}.log");
                    if ( $ret ) {die "fastrun.pl: ERROR! Couldn't copy job.stdout ".
                                     "fastrun_${gjobnr}.log: $ret, aborted!\n";}
                }
                $ret = system("pwd");
                if ( $ret ) {print "fastrun.pl: WARNING! Couldn't print cwd!\n";}
                $ret = system("ls -la");
                if ( $ret ) {print "fastrun.pl: WARNING! Couldn't list cwd!\n";}
            } elsif ( $batch eq "GRID" ) {
                grid_storage("LOGSAV","$aidir","fastrun_${gjobnr}.err","$tpath","${tfile}.err",$prot);
                grid_storage("LOGSAV","$aidir","fastrun_${gjobnr}.log","$tpath","${tfile}.log",$prot);
            } else {
# For PBS with grid-control log files are in $pwdir as job.stdout and .stderr
# ... only copy
                if ( -f "job.stderr" ) {
                    my $ret = system("cp -p job.stderr ${tfile}.err");
                    if ( $ret ) {die "fastrun.pl: Couldn't copy job.stderr ".
                                     "fastrun_${gjobnr}.err: $ret, aborted!\n";}
                }
                if ( -f "job.stdout" ) {
                    $ret = system("cp -p job.stdout ${tfile}.log");
                    if ( $ret ) {die "fastrun.pl: Couldn't copy job.stdout ".
                                     "fastrun_${gjobnr}.log: $ret, aborted!\n";}
                }
            }
            unless ( $batch eq "GC" ) {
                my $date = `date +%d%m%Y_%H%M%S`;
                chomp $date;
                print "fastrun.pl: Log files stored: LOGSAV1_$date\n";
            }
        }
        exit 0;
    }
# Check on process before entering loop
    my $user  = `id -un`;
    chomp $user;
    my $fpid  = `/bin/ps -o pid,args --no-headers | grep \"nlojet++\" | grep \"${jobnr}\" | grep -v \"grep\"`;
    chomp $fpid;
    if ( $verb ) {
        print "fastrun.pl: DEBUG! user $user\n";
        print "fastrun.pl: DEBUG! fpid1 $fpid\n";
    }
    my @tmp = split(" ",$fpid);
    $fpid = $tmp[0];
    my $keep = 0;
    if ( $fpid ) {
        if ( $verb ) {print "fastrun.pl: DEBUG! fpid2 $fpid\n";}
    } else {
        if ( $verb ) {
            print "fastrun.pl: DEBUG! fpid2 not yet there (race condidtion)?\n";
            print "fastrun.pl: DEBUG! Setting fpid anyway to enter while loop.\n";
        }
        $keep = 1;
    }

    while ( 1 ) {
#    while ( $fpid || $keep ) {
        if ( $batch eq "GRID" ) {
            sleep 1000;
            $date = `date +%d%m%Y_%H%M%S`;
            print "\nfastrun.pl: Saving table: $date\n";
            chomp $date;
            grid_storage("SAVE");
        } else {
            sleep 500;
            if ( $keep ) {
                $keep  = 0;
                $fpid  = `/bin/ps -o pid,args --no-headers | grep \"nlojet++\" | grep \"${jobnr}\" | grep -v \"grep\"`;
                chomp $fpid;
                my @tmp = split(" ",$fpid);
                $fpid = $tmp[0];
            }
            if ( $fpid ) {
                my $fptest = `/bin/ps --no-headers -p $fpid`;
                chomp $fptest;
                if ( $fptest ) {
                    if ( $verb ) {
                        print "fastrun.pl: DEBUG! fpid3  $fpid\n";
                        print "fastrun.pl: DEBUG! fptest $fptest\n";
                    }
                } else {
                    $fpid = 0;
                    $date = `date +%d%m%Y_%H%M%S`;
                    chomp $date;
                    print "\nfastrun.pl: fastNLO finished: $date\n";
                }
            } else {
                $date = `date +%d%m%Y_%H%M%S`;
                chomp $date;
                print "\nfastrun.pl: fastNLO finished: $date\n";
            }
        }
    }
}

exit 0;



#
# Signal handling and grid storage: Save table and logs on grid storage
#
# Original table naming example: fnl0002kt10_0001-hhc-born-2jet.raw/tab
# Grid table naming example:     fnl0002kt10-hhc-born-2jet_0000.raw/tab
#
sub grid_storage {
    my $signam   = shift;
    my $spath    = shift;
    my $sfile    = shift;
    my $tpath    = shift;
    my $tfile    = shift;
    my $protocol = shift;
    print "fastrun.pl: spath $spath\n";
    print "fastrun.pl: sfile $sfile\n";
    print "fastrun.pl: tpath $tpath\n";
    print "fastrun.pl: tfile $tfile\n";
    print "fastrun.pl: protocol $protocol\n";

#    my $sename = "ekp-lcg-se.physik.uni-karlsruhe.de";
#    my $sename = "dcache-se-cms.desy.de";
#    my $sename = "dcache-door-cms02.desy.de";
#    my $sepath = "/pnfs/desy.de/cms/analysis/qcd/rabbertz/";

# KITEKP via globus-url-copy
    my $gucname = "ic-kit-lcgse.rz.uni-karlsruhe.de";
    my $gucpath = "/wlcg/data/users/cms/rabbertz/";

# Aachen T2 via srm
#    my $srmname = "dcache-se-cms.desy.de:8443/srm/managerv2?SFN=";
    my $srmname = "grid-srm.physik.rwth-aachen.de:8443";
    my $srmpath = "/pnfs/physik.rwth-aachen.de/cms/store/user/krabbert/";

# DESY T2 via srm
#    my $srmname = "dcache-se-cms.desy.de:8443/srm/managerv2?SFN=";
#    my $srmname = "dcache-se-cms.desy.de:8443";
#    my $srmpath = "/pnfs/desy.de/cms/tier2/store/user/krabbert/";

    my $gucopt = "-cd";
    if ($protocol eq "guc") {
# Check globus-url-copy version
        my $gcmd = "globus-url-copy -version 2>&1";
        my $ret = `$gcmd`;
        chomp $ret;
        my @tmp = split(" ",$ret);
        my $guccmd = $tmp[0];
        my $gucver = $tmp[1];
        unless ( $gucver =~ m/\d*\.\d*/ ) {
            print "fastrun.pl: ERROR! ".
                "globus-url-copy not found: $ret!\n";}
        print "fastrun.pl: Found command $guccmd version $gucver\n";
        $gucopt = "-cd";
        if ( $gucver <= 2.9 ) {
            print "fastrun.pl: WARNING! Version of globus-url-copy ".
                "too old: $gucver. Cannot create directories on the fly!\n";
            $gucopt = "";
        }
    } else {
# Check srmcp version
        my $gcmd = "srmcp -version 2>&1 | grep -i version";
        my $ret = `$gcmd`;
        chomp $ret;
        my @tmp = split(" ",$ret);
        my $srmver = $tmp[$#tmp];
        chomp $srmver;
        print "fastrun.pl: Found srmcp command version $srmver\n";
    }

    my $sename = $gucname;
    my $sepath = $gucpath;
    if ($protocol eq "srm" ) {
        $sename = $srmname;
        $sepath = $srmpath;
    }

    my $source = "${spath}/${sfile}";
    my $target = "${sepath}${tpath}/${tfile}";

    if ( $signam eq "FETCH" ) {
        $source = "${sepath}${spath}/${sfile}";
        $target = "${tpath}/${tfile}";
    }
    $source =~ s/\/\.\//\//g;
    $target =~ s/\/\.\//\//g;

    print "fastrun.pl: Source: $source\n";
    print "fastrun.pl: Target: $target\n";

    my $gjobnr = "";
    if ( defined $ENV{MY_JOBID} ) {
        $gjobnr = $ENV{MY_JOBID};
#        $gjobnr++;
        $gjobnr = substr("0000$gjobnr",-4);
    }

    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    if ( $signam eq "TABSAV" || $signam eq "LOGSAV" ) {
        print "\nfastrun.pl: Received signal $signam: ${signam}0_$date\n";
        print "fastrun.pl: Saving source $source to\n";
        print "fastrun.pl: target $target for\n";
        print "fastrun.pl: fastNLO job no. $jobnr and\n";
        print "fastrun.pl: grid-control job no. $gjobnr\n";
        print "fastrun.pl: Trying to find source $source ...\n";
        system("pwd");
        system("ls -la $spath");
        unless ( -f "$source" ) {
            die "fastrun.pl: ERROR! Could not find source $source\n";
        }
# Create target directory if necessary
# Deprecated due to grid chaos! What does work? srm? edg/glite-gridftp-mkdir?
#                my $gcmd = "edg-gridftp-mkdir ".
#                    "gsiftp://${sename}/${sepath}".
#                    "fastNLO_tables/${tabdir}";
#                print "Command $gcmd\n";
#                my $ret = system("$gcmd");
#                if ( $ret ) {print "fastrun.pl: WARNING! Creation of grid storage directory failed: $ret!\n";}
# Create target directory directly with globus-url-copy if possible
# Don't bother to copy empty error files ...
        if ( -z "$source" ) {
            print "fastrun.pl: WARNING! Skipping storage of empty file $source\n";
        } else {
            if ($protocol eq "guc") {
                my $gcmd = "globus-url-copy $gucopt ".
                    "file://${source} ".
                    "gsiftp://${sename}/${target}";
                print "Command $gcmd\n";
                my $ret = system("$gcmd");
                if ( $ret ) {die "fastrun.pl: ERROR! Grid storage of ".
                                 "source $source to target $target failed: $ret!\n";}
            } else {
                $ENV{SRM_PATH} = "";
                my $gcmd = "srmcp ".
                    "file:////${source} ".
                    "srm://${sename}/${target}";
                print "Command $gcmd\n";
                my $ret = system("$gcmd");
                if ( $ret ) {die "fastrun.pl: ERROR! Grid storage of ".
                                 "source $source to target $target failed: $ret!\n";}
            }
        }
    } elsif ( $signam eq "FETCH" ) {
        print "fastrun.pl: Fetching binary fastNLO archive $sfile from\n";
        print "fastrun.pl: SE $sename in path $sepath\n";
        print "fastrun.pl: Trying to fetch binary fastNLO archive $source on\n";
        print "fastrun.pl: SE $sename in path $sepath\n";
        my $gcmd = "globus-url-copy ".
            "gsiftp://${sename}/${source} ".
            "file://`pwd`/${target}";
        if ($protocol eq "srm") {
            $ENV{SRM_PATH} = "";
            $gcmd = "srmcp -streams_num=1 ".
                "srm://${sename}/${source} ".
                "file:///`pwd`/${target}";
        }
        print "Command $gcmd\n";
        system("$gcmd");
        system("pwd");
        system("ls -la");
    } else {
        die "fastrun.pl: ERROR! Caught unsupported signal $signam, aborted!\n";
    }

    return 0;
}
