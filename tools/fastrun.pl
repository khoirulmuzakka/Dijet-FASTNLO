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
chomp $date;
print "\n######################################\n";
print "# fastrun.pl: Starting run of fastNLO: FASTRUN0_$date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_f, $opt_h, $opt_i, $opt_j,
      $opt_m, $opt_o, $opt_p, $opt_r, $opt_s, $opt_t, $opt_v ) =
    ( "LOCAL", "", "0", "187", "", ".", "0001",
      "0", "LO", "CTEQ", "", ".", "", "1" );
getopts('b:de:f:hi:j:m:o:p:rs:t:v:') or die "fastrun.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastrun.pl\n";
    print "Usage: fastrun.pl [switches/options] ([scenario])\n";
    print "  -b batch        Batch system used: LOCAL (def.), GRID or PBS\n";
    print "  -d debug        Switch debug/verbose mode on\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -f rev          fastNLO revision to use (def.=187)\n";
    print "  -h              Print this text\n";
    print "  -i dir          Installation directory (def.=.)\n";
    print "  -j jobnr        Job number to attach (def.=0001)\n";
    print "  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -p pdf          CTEQ parton densities (def.) or LHAPDF\n";
    print "  -r              Reference calculation incl. pdf access\n";
    print "  -s dir          Archive source directory (def.=.)\n";
    print "  -t dir          Output target directory: ".
	"(def.= {scen}{ref}_{jobnr} with\n                  ".
	"ref. to working directory in fastNLO installation)\n";
    print "  -v #            Choose between fastNLO version 1 and 2 (def.=1)\n";
    print "\n";
    print "Examples:\n";
    print "1) Install only (to install with LHAPDF use option -p):\n";
    print "   ./fastrun.pl [-i .|installdir] [-f 187|rev] -m 1 [-p CTEQ|LHAPDF] [-s .|sdir] [-v 1|2]\n\n";
    print "2) Make only scenario (to make scenario for reference mode use option -r):\n";
    print "   ./fastrun.pl [-i .|installdir] [-f 187|rev] -m 2 [-p CTEQ|LHAPDF] [-r] scenarioname\n\n";
    print "3) Run only (to run scenario in reference mode use option -r):\n";
    print "   ./fastrun.pl [-b LOCAL|GRID|batch] [-i .|installdir] [-e max-events] [-f 187|rev] -m 3 [-p CTEQ|LHAPDF] [-r] [-t ./{scen}{ref}_{jobnr}|tdir] scenarioname\n\n";
    exit;
}

unless ( $opt_b eq "LOCAL" || $opt_b eq "GRID" || $opt_b eq "PBS" ) {
    die "fastrun.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( $opt_f =~ m/\d+/ && $opt_f !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal fastNLO revision number: $opt_f, aborted.\n";
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
unless ( $opt_o eq "LO" || $opt_o eq "NLO" ) {
    die "fastrun.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( $opt_p eq "CTEQ" || $opt_p eq "LHAPDF" ) {
    die "fastrun.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( -d $opt_s || -l $opt_s ) {
    die "fastrun.pl: Error! Archive source directory $opt_s does not exist, aborted.\n";
}
unless ( $opt_v == 1 || $opt_v == 2) {
    die "fastrun.pl: Error! fastNLO version $opt_v does not exist, aborted.\n";
}

my $batch = $opt_b;
my $idir  = $opt_i;
my $nmax  = $opt_e;
my $frev  = $opt_f;
my $jobnr = $opt_j;
my $mode  = $opt_m;
my $order = $opt_o;
my $pdf   = $opt_p;
my $ref   = "";
if ( $opt_r ) { $ref = "ref";}
my $sdir  = $opt_s;
my $verb  = "";
my $vers  = $opt_v;
print "fastrun.pl: Directory for/of installation is $idir\n";
print "fastrun.pl: Using fastNLO revision $frev\n";
print "fastrun.pl: Job mode is $mode\n";
unless ( $mode == 1 || $mode == 2 ) {
    print "fastrun.pl: Running on batch system $batch\n";
    print "fastrun.pl: Maximal event number: $nmax\n";
    print "fastrun.pl: Attaching job number $jobnr\n";
    print "fastrun.pl: Running in order $order\n";
    if ( $ref eq "ref" ) {
	print "fastrun.pl: Running with pdf $pdf\n";
    }
}
if ( $mode != 1 && $ref ) {
    print "fastrun.pl: Making/running in reference mode\n";
}
if ( $mode == 0 || $mode ==1 ) {
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
my $tabnam = "${scen}${ref}_${jobnr}-hhc-$runmode{$order}[0]-2jet.raw";

# Directories
my $pwdir = getcwd();
chdir $idir;
my $aidir = getcwd();
chdir $pwdir;

# Define install hash
my %install;
# 0: Archive names
# 1: Unpacking directory of the archives
# 2: Final directories
# 3: Link shortcut
#
# Store used gcc version to add this to final directories
my $gccvers = "3.3.6";
if ( $vers == 2 ) {
    $gccvers = `gcc -dumpversion`;
    chomp $gccvers;
}
print "fastrun.pl: Using gcc compiler version $gccvers\n"; 

$install{fastjet}[0]    = "fastjet-2.3.4";
$install{fastjet}[1]    = "fastjet-2.3.4";
$install{mcfm}[0]       = "mcfm-5.3";
$install{mcfm}[1]       = "MCFM-5.3";
$install{fastNLO}[0]    = "fastNLO-rev${frev}";
$install{fastNLO}[1]    = "fastNLO-rev${frev}";
if ( $vers == 1 ) { 
    $install{gcccore}[0]    = "gcc-core-$gccvers";
    $install{gcccore}[1]    = "gcc-$gccvers";
    $install{gccgpp}[0]     = "gcc-g\+\+-$gccvers";
    $install{gccgpp}[1]     = "gcc-$gccvers";
    $install{gccg77}[0]     = "gcc-g77-$gccvers";
    $install{gccg77}[1]     = "gcc-$gccvers";
    $install{cernlib}[0]    = "cernlib-2003";
    $install{cernlib}[1]    = "cernlib-2003";
# Versions >= 5.4.0 don't work with gcc 3.3.6
    $install{lhapdf}[0]     = "lhapdf-5.3.1";
    $install{lhapdf}[1]     = "lhapdf-5.3.1";
    $install{nlojet}[0]     = "nlojet++-2.0.1";
    $install{nlojet}[1]     = "nlojet++-2.0.1";
    $install{nlojetfix}[0]  = "nlojet++-2.0.1-fix";
    $install{nlojetfix}[1]  = "nlojet++-2.0.1";
} else {
    $install{root}[0]       = "root-5.22";
    $install{root}[1]       = "root";
    $install{lhapdf}[0]     = "lhapdf-5.6.0";
    $install{lhapdf}[1]     = "lhapdf-5.6.0";
    $install{nlojet}[0]     = "nlojet++-4.0.1";
    $install{nlojet}[1]     = "nlojet++-4.0.1";
}

foreach my $comp ( keys %install ) {
    $install{$comp}[0] .=".tar.gz";
    $install{$comp}[2] = $install{$comp}[1]."-$gccvers";
    $install{$comp}[3] = $install{$comp}[1];"

    unless ( $mode == 2 || $mode == 3 ||
	     -d "$idir/$install{$comp}[0]" ||
	     -f "$sdir/$install{$comp}[1]" ||
	     ( $comp eq "lhapdf" && $pdf eq "CTEQ" ) ) {
	die "fastrun.pl: Archive $install{$comp}[1] does not exist, ".
	    "aborted!\n";
    }
}
exit 0;

#
# Print system info
#
if ( $verb ) {
    print "\n######################################################\n";
    print "fastrun.pl: System information for debugging purposes:\n";
    print "######################################################\n";
    my $host = `hostname`;
    chomp $host;
    print "fastrun.pl: The system's hostname is (hostname):\n$host\n\n"; 
    my $osvers = `uname -a`;
    print "fastrun.pl: Your operating system is (uname -a):\n$osvers\n";
    my $procvers = `cat /proc/version`;
    print "fastrun.pl: Your linux version is (/proc/version):\n$procvers\n";
    my $freemem = `free`;
    print "fastrun.pl: The available memory is:\n$freemem\n";
    my $cpumod = `cat /proc/cpuinfo | grep \"model name\"`;
    my $cpufrq = `cat /proc/cpuinfo | grep \"cpu MHz\"`;
    print "fastrun.pl: The processor type is:\n${cpumod}at\n${cpufrq}\n";
    my $freedisk = `df -h`;
    print "fastrun.pl: The available disk space is:\n$freedisk\n";
    my $freenode = `df -hi`;
    print "fastrun.pl: The available inode space is:\n$freenode\n";
#print "fastrun.pl: Installation environment:\n";
#system("printenv");
    print "######################################################\n";
}

#
# 0) Install gcc
#
if ( $mode == 0 || $mode == 1 ) {
    if ( $vers == 1 ) {
# Standard environment to use system compiler still
	unless ( -e "$idir/$install{gcccore}[0]" ) {
	    $date = `date +%d%m%Y_%H%M%S`;
	    chomp $date;
	    print "\nfastrun.pl: Unpacking gcc-core sources in $install{gcccore}[1]: $date\n";
	    my $ret = system("tar xz -C $idir -f $sdir/$install{gcccore}[1]");
	    if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gcccore}[1] ".
			     "in $idir failed: $ret, aborted!\n";}
	    if ( -l "$idir/gcc" ) {system("rm -f $idir/gcc");}
	    system("ln -s $install{gcccore}[0] $idir/gcc");
	    chdir "$idir/gcc";
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
#	my @acmd = split(" ",$cmd);
#	system(@acmd);
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
	    $ret = system("tar xz -C $idir -f $sdir/$install{gccgpp}[1]");
	    if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gccgpp}[1] ".
			     "in $idir failed: $ret, aborted!\n";}
	    $ret = system("tar xz -C $idir -f $sdir/$install{gccg77}[1]");
	    if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{gccg77}[1] ".
			     "in $idir failed: $ret, aborted!\n";}
	    chdir "$idir/gcc";
	    system("make clean");
	    print "fastrun.pl: Configuring gcc for recompile (incl. gcc): $cmd ...\n";
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
	$ENV{GCC_EXEC_PREFIX} = "$aidir/lib/gcc-lib/";
    }	

#
# 1) Install CERN libraries (V1 only) resp. ROOT (V2)
#
    if ( $vers == 1 ) {
	unless ( -e "$idir/$install{cernlib}[0]" ) {
	    $date = `date +%d%m%Y_%H%M%S`;
	    chomp $date;
	    print "\nfastrun.pl: Unpacking CERN libraries in $install{cernlib}[1]: $date\n";
	    my $ret =system("tar xz -C $idir -f $sdir/$install{cernlib}[1]");
	    if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{cernlib}[1] ".
			     "in $idir failed: $ret, aborted!\n";}
	    if ( -l "$idir/cernlib" ) {system("rm -f $idir/cernlib");}
	    system("ln -s $install{cernlib}[0] $idir/cernlib");
	}
    } else {
	unless ( -e "$idir/$install{root}[0]" ) {
	    $date = `date +%d%m%Y_%H%M%S`;
	    chomp $date;
	    print "\nfastrun.pl: Installing ROOT from $install{root}[1]: $date\n";
	    print "\nfastrun.pl: Unpacking $install{root}[1] ...\n";
	    my $ret = system("tar xz -C $idir -f $sdir/$install{root}[1]");
	    if ( $ret ) {die "fastrun.pl: Unpacking of archive ".
			     "$sdir/$install{root}[1] ".
			     "in $idir failed: $ret, aborted!\n";}
	    $ret = system("mv $install{root}[2] $idir/$install{root}[0]");
	    if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
			     "$install{root}[2] to ".
			     "$idir/$install{root}[0]: $ret, aborted!\n";}
	    if ( -l "$idir/root" ) {system("rm -f $idir/root");}
	    system("ln -s $install{root}[0] $idir/root");
	    chdir "$idir/$install{root}[0]";
	    print "\nfastrun.pl: Configuring ROOT ...\n";
	    my $cmd = "./configure --enable-roofit --prefix=`pwd` --etcdir=`pwd`/etc --with-cxx=\"gcc -Wall\"\n";
	    if ( $verb ) {print "$cmd"};
	    $ret = system("./configure --enable-roofit --prefix=`pwd` --etcdir=`pwd`/etc --with-cxx=\"gcc -Wall\"");
	    if ( $ret ) {die "fastrun.pl: Error $ret in ROOT configure step, aborted!\n";}
	    print "\nfastrun.pl: Making ROOT ...\n";
# No multithreaded make for ROOT :-( 
#	    $ret = system("make -j2");
	    $ret = system("make");
	    if ( $ret ) {die "fastrun.pl: Error $ret in ROOT make step, aborted!\n";}
	    print "\nfastrun.pl: Make install for ROOT ...\n";
	    $ret = system("make install");
	    if ( $ret ) {die "fastrun.pl: Error $ret in ROOT make install step, aborted!\n";}
	    chdir "$pwdir";
	}
    }

#
# 2) Install lhapdf
#
    unless ( -e "$idir/$install{lhapdf}[0]" || $pdf eq "CTEQ" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Installing lhapdf from $install{lhapdf}[1]: $date\n";
	print "\nfastrun.pl: Unpacking $install{lhapdf}[1] ...\n";
	my $ret = system("tar xz -C $idir -f $sdir/$install{lhapdf}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{lhapdf}[1] ".
			 "in $idir failed: $ret, aborted!\n";}
	$ret = system("mv $install{lhapdf}[2] $idir/$install{lhapdf}[0]");
	if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
			 "$install{lhapdf}[2] to ".
			 "$idir/$install{root}[0]: $ret, aborted!\n";}
	if ( -l "$idir/lhapdf" ) {system("rm -f $idir/lhapdf");}
	system("ln -s $install{lhapdf}[0] $idir/lhapdf");
	chdir "$idir/$install{lhapdf}[0]";
	print "\nfastrun.pl: Configuring lhapdf ...\n";
#	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
# Avoid Python Murks in LHAPDF >= 5.4.0
	$ret = system("./configure --prefix=`pwd` --disable-pyext");
#	$ret = system("./configure --prefix=`pwd`");
	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF configure step, aborted!\n";}
	print "\nfastrun.pl: Making lhapdf ...\n";
# At least until LHAPDF 5.3.1 a race condition spoils multithreaded make!
#	$ret = system("make");
# Works now
	$ret = system("make -j2");
	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make step, aborted!\n";}
	print "\nfastrun.pl: Make install for lhapdf ...\n";
	$ret = system("make install");
# In addition, there is a double lhapdf-config creation => ignore this error
#	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make install step, aborted!\n";}
	if ( $ret ) {print "fastrun.pl: Error $ret in LHAPDF make install step ignored!\n";}
	chdir "$pwdir";
    }

#
# 3) Install fastjet
#
    unless ( -e "$idir/$install{fastjet}[0]" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Installing fastjet from $install{fastjet}[1]: $date\n";
	print "\nfastrun.pl: Unpacking $install{fastjet}[1] ...\n";
	my $ret = system("tar xz -C $idir -f $sdir/$install{fastjet}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastjet}[1] ".
			 "in $idir failed: $ret, aborted!\n";}
	$ret = system("mv $install{fastjet}[2] $idir/$install{fastjet}[0]");
	if ( $ret ) {die "fastrun.pl: Couldn't move unpacking dir ".
			 "$install{fastjet}[2] to ".
			 "$idir/$install{fastjet}[0]: $ret, aborted!\n";}
	if ( -l "$idir/fastjet" ) {system("rm -f $idir/fastjet");}
	system("ln -s $install{fastjet}[0] $idir/fastjet");
	chdir "$idir/$install{fastjet}[0]";
	print "\nfastrun.pl: Configuring fastjet ...\n";
	$ret = system("./configure --enable-shared --prefix=`pwd` --bindir=$aidir/bin");
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
	chdir "$pwdir";
    }


#
# 4) Install Nlojet++
#
# NLOJet++ Version 2
    unless ( -e "$idir/$install{nlojet}[0]" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Installing Nlojet++ from $install{nlojet}[1]: $date\n";
	print "\nfastrun.pl: Unpacking $install{nlojet}[1] ...\n";
	my $ret = system("tar xz -C $idir -f $sdir/$install{nlojet}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojet}[1] ".
			 "in $idir failed: $ret, aborted!\n";}
	print "\nfastrun.pl: Unpacking fix for $install{nlojet}[1] ...\n";
	$ret = system("tar xzv -C $idir -f $sdir/$install{nlojetfix}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojetfix}[1] ".
			 "in $idir failed: $ret, aborted!\n";}
	if ( -l "$idir/nlojet" ) {system("rm -f $idir/nlojet");}
	system("ln -s  $install{nlojet}[0] $idir/nlojet");
	print "0: $install{nlojet}[0]\n";
	print "1: $install{nlojet}[1]\n";
	print "2: $install{nlojet}[2]\n";
	chdir "$idir/$install{nlojet}[2]";
	print "\nfastrun.pl: Configuring Nlojet++ ...\n";
#	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
	$ret = system("./configure --prefix=`pwd`");
	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ configure step, aborted!\n";}
	print "\nfastrun.pl: Making Nlojet++ ...\n";
	$ret = system("make -j2 CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make step, aborted!\n";}
	print "\nfastrun.pl: Make install for Nlojet++ ...\n";
	$ret = system("make install CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ make install step, aborted!\n";}
	chdir "$pwdir";
    }
# NLOJet++ Version 4
#    unless ( -e "$idir/$install{nlojet4}[0]" ) {
#	$date = `date +%d%m%Y_%H%M%S`;
#	chomp $date;
#	print "\nfastrun.pl: Installing Nlojet++ version 4 from $install{nlojet4}[1]: $date\n";
#	print "\nfastrun.pl: Unpacking $install{nlojet4}[1] ...\n";
#	my $ret = system("tar xz -C $idir -f $sdir/$install{nlojet4}[1]");
#	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{nlojet4}[1] ".
#			 "in $idir failed: $ret, aborted!\n";}
#	if ( -l "$idir/nlojet4" ) {system("rm -f $idir/nlojet4");}
#	system("ln -s  $install{nlojet4}[0] $idir/nlojet4");
#	chdir "$idir/$install{nlojet4}[0]";
#	print "\nfastrun.pl: Configuring Nlojet++ version 4 ...\n";
##	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
#	$ret = system("./configure --prefix=`pwd`");
#	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ version 4 configure step, aborted!\n";}
#	print "\nfastrun.pl: Making Nlojet++ version 4 ...\n";
#	$ret = system("make -j2 CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
#	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ version 4 make step, aborted!\n";}
#	print "\nfastrun.pl: Make install for Nlojet++ version 4 ...\n";
#	$ret = system("make install CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
#	if ( $ret ) {die "fastrun.pl: Error $ret in NLOJET++ version 4 make install step, aborted!\n";}
#	chdir "$pwdir";
#    }

#
# 5) Install mcfm
#
# CERNLIB var already needed for mcfm    
    print "\nfastrun.pl: Setting environment variable CERNLIB for mcfm:\n";
    my $cwdir = getcwd();
    $ENV{CERNLIB}  = "$cwdir/cernlib"; 
    unless ( -e "$idir/$install{mcfm}[0]" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Installing mcfm from $install{mcfm}[1]: $date\n";
	print "\nfastrun.pl: Unpacking $install{mcfm}[1] ...\n";
	system ("mkdir $idir/$install{mcfm}[0]");
	my $ret = system("tar xz -C $idir/$install{mcfm}[0] -f $sdir/$install{mcfm}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{mcfm}[1] ".
			 "in $idir/$install{mcfm}[0] failed: $ret, aborted!\n";}
	print "\nfastrun.pl: Unpacking fix for $install{mcfm}[1] ...\n";
	$ret = system("tar xzv -C $idir/$install{mcfm}[0] -f $sdir/$install{mcfmfix}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{mcfmfix}[1] ".
			 "in $idir/$install{mcfm}[0] failed: $ret, aborted!\n";}
	if ( -l "$idir/mcfm" ) {system("rm -f $idir/mcfm");}
	system("ln -s  $install{mcfm}[0] $idir/mcfm");
	chdir "$idir/$install{mcfm}[0]";
	print "\nfastrun.pl: Configuring mcfm ...\n";
	$ret = system("./Install");
	if ( $ret ) {die "fastrun.pl: Error $ret in mcfm install script, aborted!\n";}
	$ret = system("make -j2");
	if ( $ret ) {die "fastrun.pl: Error $ret in mcfm make step, aborted!\n";}
        chdir "..";
    }
    exit 0;

#
# 6) Install fastNLO
#
    unless ( -e "$idir/$install{fastNLO}[0]" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Installing fastNLO from $install{fastNLO}[1]: $date\n";
	print "\nfastrun.pl: Unpacking $install{fastNLO}[1] ...\n";
	my $ret = system("tar xz -C $idir -f $sdir/$install{fastNLO}[1]");
	if ( $ret ) {die "fastrun.pl: Unpacking of archive $sdir/$install{fastNLO}[1] ".
			 "in $idir failed: $ret, aborted!\n";}
	if ( -l "$idir/fastNLO" ) {system("rm -f $idir/fastNLO");}
	system("ln -s $install{fastNLO}[0] $idir/fastNLO");
    }
}

#
# Set environment that is expected to have been set up and
# save it to file fastnlo.(c)sh for later "source"ing.
#
print "\nfastrun.pl: Setting fastNLO and gcc environment to be used ".
    "with this installation!\n";
print "fastrun.pl: Check fastNLO.(c)sh for reference and ".
    "later \"source\"ing.\n";
open (FILE,"> fastNLO.csh");
my $cwdir = getcwd();
print FILE "#\n";
print FILE "# fastrun.pl: Setting fastNLO and gcc environment to be used ".
    "with this installation!\n";
print FILE "#\n";
print FILE "setenv CERNLIB $cwdir/cernlib\n";
print FILE "setenv FASTJET $cwdir/fastjet\n";
print FILE "setenv FASTNLO $cwdir/fastNLO\n";
print FILE "setenv LHAPDF  $cwdir/lhapdf/lib\n";
print FILE "setenv NLOJET  $cwdir/nlojet\n";
print FILE "setenv NLOJET4 $cwdir/nlojet4\n";
$ENV{CERNLIB}  = "$cwdir/cernlib"; 
$ENV{FASTJET}  = "$cwdir/fastjet";
$ENV{FASTNLO}  = "$cwdir/fastNLO";
$ENV{LHAPDF}   = "$cwdir/lhapdf/lib";
$ENV{NLOJET}   = "$cwdir/nlojet";
$ENV{NLOJET4}  = "$cwdir/nlojet4";
print FILE "setenv PATH $cwdir/bin:$ENV{FASTJET}:".
    "$ENV{NLOJET4}/bin:\${PATH}".
    "$ENV{NLOJET}/bin:\${PATH}\n";
if ( $ENV{LD_LIBRARY_PATH} ) {
    print FILE "setenv LD_LIBRARY_PATH $cwdir/lib:$cwdir/lib64:".
	"$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
	"$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
	"$ENV{FASTJET}/plugins/CDFCones/.libs:".
	"$ENV{NLOJET}/lib:$ENV{NLOJET4}/lib:".
	"$ENV{LHAPDF}:\${LD_LIBRARY_PATH}\n";
    $ENV{LD_LIBRARY_PATH} ="$cwdir/lib:$cwdir/lib64:".
	"$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
	"$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
	"$ENV{FASTJET}/plugins/CDFCones/.libs:".
	"$ENV{NLOJET}/lib:$ENV{NLOJET4}/lib:".
	"$ENV{LHAPDF}:$ENV{LD_LIBRARY_PATH}";
} else {
    print FILE "setenv LD_LIBRARY_PATH $cwdir/lib:$cwdir/lib64:".
	"$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
	"$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
	"$ENV{FASTJET}/plugins/CDFCones/.libs:".
	"$ENV{NLOJET}/lib:$ENV{NLOJET4}/lib:".
	"$ENV{LHAPDF}\n";
    $ENV{LD_LIBRARY_PATH} ="$cwdir/lib:$cwdir/lib64:".
	"$ENV{FASTJET}/lib:$ENV{FASTJET}/plugins/SISCone/.libs:".
	"$ENV{FASTJET}/plugins/SISCone/siscone/siscone/.libs:".
	"$ENV{FASTJET}/plugins/CDFCones/.libs:".
	"$ENV{NLOJET}/lib:$ENV{NLOJET4}/lib:".
	"$ENV{LHAPDF}";
}
print FILE "setenv GCC_EXEC_PREFIX $cwdir/lib/gcc-lib/\n";
print FILE "setenv CXXFLAGS \"-O3 -I .\"\n";
$ENV{PATH}            = "$cwdir/bin:$ENV{FASTJET}:".
    "$ENV{NLOJET4}/bin:$ENV{NLOJET}/bin:$ENV{PATH}";
$ENV{GCC_EXEC_PREFIX} ="$cwdir/lib/gcc-lib/";
$ENV{CXXFLAGS} = "-O3 -I .";
close (FILE);
if ( $verb ) {
    my $tmp = `which gcc`;
    print "fastrun.pl: DEBUG! gcc executable used: $tmp\n";
}

#
# 7) Make fastNLO scenario
#
my $scendir;
if ( $mode == 0 || $mode == 2 ) {
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Making fastNLO scenario: $date\n";
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
}
    
#
# 8) Run fastNLO
#
if ( $mode == 0 || $mode == 3 ) {
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Running fastNLO scenario on batch system $batch: $date\n";
    
# a) Run on grid without complete installation by source
    if ( $batch eq "GRID" && $mode == 3 ) {
# Fetching and unpacking of fastNLO binary archive
	my $file = "fastNLO-bin";
	if ( $ref ) {
	    $file .= "-${pdf}";
	}
	$file .= ".tgz";
	if ( ! -f $file ) {
	    grid_storage("FETCH","fastNLO_archives",$file,".",$file);
	}
	if ( -f $file ) {
	    system ("tar xfz $file");
	} else {
	    die "fastrun.pl: ERROR! Could not find binary tgz of fastNLO $file, aborted!\n";
	}
    }
    
# Structure change in fastNLO following change in revision 212!
    if ( $frev < 212 ) { 
	$scendir = "$ENV{FASTNLO}/author1c/hadron";
    } else {
	$scendir = "$ENV{FASTNLO}/trunk/v1.4/author1c/hadron";
    }
    chdir "$scendir" or die
	"fastrun.pl: ERROR! Could not cd to dir $scendir!\n";
    
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Running fastNLO: $date\n";
    my $cmd = "$ENV{NLOJET}/bin/nlojet++ -P dipole ".
	"--save-after $runmode{$order}[1] ".
	"-c$runmode{$order}[0] ".
	"-d $tdir ".
	"-n ${scen}${ref}_${jobnr} ".
	"-u ${scen}${ref}.la ";
    if ( $nmax ) { $cmd .= "--max-event $nmax"; }
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
	if ( $batch eq "GRID" ) {
	    my $spath = "${scendir}/${tdir}";
	    my $sfile = "${tabnam}";
	    my $tpath = "fastNLO_tables/${scen}${ref}";
	    if ( $tpath ) {
		$tpath =~ s/\/\.\//\//g;
	    } else {
		$tpath = "";
	    }
# Change job numbering according to grid-control
	    my $tfile  = $sfile;
	    $tfile =~ s/_0001//;
	    $tfile =~ s/\.raw/_${gjobnr}\.raw/;
	    grid_storage("TABSAV","$spath","$sfile","$tpath","$tfile");
	    my $date = `date +%d%m%Y_%H%M%S`;
	    chomp $date;
	    print "fastrun.pl: Table stored: TABSAV1_$date\n";
	}
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\n###############################\n";
	print "# fastrun.pl: fastNLO finished: FASTRUN1_$date\n";
	print "###############################\n\n";
#	close STDERR;
#	close STDOUT;
# Copy log files to grid storage
	if ( $batch eq "GRID" ) {
	    my $tpath = "fastNLO_tables/${scen}${ref}";
	    if ( $tpath ) {
		$tpath =~ s/\/\.\//\//g;
	    } else {
		$tpath = "";
	    }
	    my $tfile = "${scen}${ref}-hhc-$runmode{$order}[0]-2jet_${gjobnr}";
	    grid_storage("LOGSAV","$cwdir","fastrun_${gjobnr}.err","$tpath","${tfile}.err");
	    grid_storage("LOGSAV","$cwdir","fastrun_${gjobnr}.log","$tpath","${tfile}.log");
	    my $date = `date +%d%m%Y_%H%M%S`;
	    chomp $date;
	    print "fastrun.pl: Log files stored: LOGSAV1_$date\n";
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
# Original table naming example: fnl0002kt10_0001-hhc-born-2jet.raw
# Grid table naming example:     fnl0002kt10-hhc-born-2jet_0000.raw
#
sub grid_storage {
    my $signam = shift;
    my $spath  = shift;
    my $sfile  = shift;
    my $tpath  = shift;
    my $tfile  = shift;
    print "fastrun.pl: spath $spath\n";
    print "fastrun.pl: sfile $sfile\n";
    print "fastrun.pl: tpath $tpath\n";
    print "fastrun.pl: tfile $tfile\n";

#    my $sename = "ekp-lcg-se.physik.uni-karlsruhe.de";
#    my $sename = "dcache-se-cms.desy.de";
#    my $sename = "dcache-door-cms02.desy.de";
#    my $sepath = "/pnfs/desy.de/cms/analysis/qcd/rabbertz/";
    my $sename = "ic-kit-lcgse.rz.uni-karlsruhe.de";
    my $sepath = "/wlcg/data/users/cms/rabbertz/";
    
    my $source = "${spath}/${sfile}";
    my $target = "${sepath}${tpath}/${tfile}";
    if ( $signam eq "FETCH" ) {
	$source = "${sepath}${spath}/${sfile}";
	$target = "${tpath}/${tfile}";
    }
    $source =~ s/\/\.\//\//g;
    $target =~ s/\/\.\//\//g;
    
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
	my $gucopt = "-cd";
	if ( $gucver <= 2.9 ) {
	    print "fastrun.pl: WARNING! Version of globus-url-copy ".
		"too old: $gucver. Cannot create directories on the fly!\n";
	    $gucopt = "";
	}
# Create target directory if necessary
# Deprecated due to grid chaos! What does work? srm? edg/glite-gridftp-mkdir?
#		my $gcmd = "edg-gridftp-mkdir ".
#		    "gsiftp://${sename}/${sepath}".
#		    "fastNLO_tables/${tabdir}";
#		print "Command $gcmd\n";
#		my $ret = system("$gcmd");
#		if ( $ret ) {print "fastrun.pl: WARNING! Creation of grid storage directory failed: $ret!\n";}
# Create target directory directly with globus-url-copy if possible
	$gcmd = "globus-url-copy $gucopt ".
	    "file://${source} ".
	    "gsiftp://${sename}/${target}";
	print "Command $gcmd\n";
	$ret = system("$gcmd");
	if ( $ret ) {die "fastrun.pl: ERROR! Grid storage of ".
			 "source $source to target $target failed: $ret!\n";}
    } elsif ( $signam eq "FETCH" ) {
	print "fastrun.pl: Fetching binary fastNLO archive $sfile from\n";
	print "fastrun.pl: SE $sename in path $sepath\n";
	print "fastrun.pl: Trying to fetch binary fastNLO archive $source on\n";
	print "fastrun.pl: SE $sename in path $sepath\n";
	my $gcmd = "globus-url-copy ".
	    "gsiftp://${sename}/${source} ". 
	    "file://`pwd`/${target}";
	print "Command $gcmd\n";
	system("$gcmd");
	system("pwd");
	system("ls -la");
    } else {
	die "fastrun.pl: ERROR! Caught unsupported signal $signam, aborted!\n";
    }
    
    return 0;
}
