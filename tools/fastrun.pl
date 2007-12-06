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
my $date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n######################################\n";
print "# fastrun.pl: Starting run of fastNLO: $date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_f, $opt_h, $opt_j,
      $opt_m, $opt_o, $opt_p, $opt_r, $opt_s, $opt_t, $opt_v ) =
    ( "LOCAL", ".", "0", "187", "", "0001",
      "0", "LO", "CTEQ", "", ".", "", "" );
getopts('b:d:e:f:hj:m:o:p:rs:t:v') or die "fastrun.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastrun.pl\n";
    print "Usage: fastrun.pl [switches/options] ([scenario])\n";
    print "  -b batch        Batch system used: LOCAL (def.), GRID or PBS\n";
    print "  -d dir          Installation directory (def.=.)\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -f rev          fastNLO revision to use (def.=187)\n";
    print "  -h              Print this text\n";
    print "  -j jobnr        Job number to attach (def.=0001)\n";
    print "  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -p pdf          CTEQ parton densities (def.) or LHAPDF\n";
    print "  -r              Reference calculation incl. pdf access\n";
    print "  -s dir          Archive source directory (def.=.)\n";
    print "  -t dir          Output target directory: ".
	"(def.= {scen}{ref}_{jobnr} with\n                  ".
	"ref. to working directory in fastNLO installation)\n";
    print "  -v              Switch verbose mode on\n";
    print "\n";
    print "Examples:\n";
    print "1) Install only (to install with LHAPDF use option -p):\n";
    print "   ./fastrun.pl [-d .|installdir] [-f 187|rev] -m 1 [-p CTEQ|LHAPDF] [-s .|sdir]\n\n";
    print "2) Make only scenario (to make scenario for reference mode use option -r):\n";
    print "   ./fastrun.pl [-d .|installdir] [-f 187|rev] -m 2 [-p CTEQ|LHAPDF] [-r] scenarioname\n\n";
    print "3) Run only (to run scenario in reference mode use option -r):\n";
    print "   ./fastrun.pl [-b LOCAL|GRID|batch] [-d .|installdir] [-e max-events] [-f 187|rev] -m 3 [-p CTEQ|LHAPDF] [-r] [-t ./{scen}{ref}_{jobnr}|tdir] scenarioname\n\n";
    exit;
}

unless ( $opt_b eq "LOCAL" || $opt_b eq "GRID" || $opt_b eq "PBS" ) {
    die "fastrun.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( -d $opt_d ) {
    die "fastrun.pl: Error! No such directory: $opt_d, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( $opt_f =~ m/\d+/ && $opt_f !~ m/\D+/ ) {
    die "fastrun.pl: Error! Illegal fastNLO revision number: $opt_f, aborted.\n";
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

my $batch = $opt_b;
my $idir  = $opt_d;
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
if ( $opt_v ) {
    $verb = 1;
    print "fastrun.pl: Verbose mode is active.\n";
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
# First entry (index 0!): Subdirecory name into which the archive is unpacked!
$install{cernlib}[0]    = "cernlib-2003";
$install{lhapdf}[0]     = "lhapdf-5.3.1";
$install{nlojet}[0]     = "nlojet++-2.0.1";
$install{nlojetfix}[0]  = "nlojet++-2.0.1";
$install{fastNLO}[0]    = "fastNLO-rev${frev}";
$install{gcccore}[0]    = "gcc-3.3.6";
$install{gccgpp}[0]     = "gcc-3.3.6";
$install{gccg77}[0]     = "gcc-3.3.6";
# Second: Archive filenames 
foreach my $comp ( keys %install ) {
    my $tmp = $install{$comp}[0].".tar.gz";
    if ( $comp eq "nlojetfix" ) {
	$install{$comp}[1] = "nlojet++-2.0.1-fix.tar.gz";
    } elsif ( $comp eq "gcccore" ) {
	$install{$comp}[1] = "gcc-core-3.3.6.tar.gz";
    } elsif ( $comp eq "gccgpp" ) {
	$install{$comp}[1] = "gcc-g\+\+-3.3.6.tar.gz";
    } elsif ( $comp eq "gccg77" ) {
	$install{$comp}[1] = "gcc-g77-3.3.6.tar.gz";
    } else {
	$install{$comp}[1] = $tmp;
    }
    unless ( $mode == 2 || $mode == 3 ||
	     -d "$idir/$install{$comp}[0]" ||
	     -f "$sdir/$install{$comp}[1]" ||
	     ( $comp eq "lhapdf" && $pdf eq "CTEQ" ) ) {
	die "fastrun.pl: Archive $install{$comp}[1] does not exist, ".
	    "aborted!\n";
    }
}

#
# Print system info
#
if ( $verb ) {
    print "\nfastrun.pl: System information for debugging purposes:\n";
    print "\n############################################################\n";
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
    print "#**********************************************************#\n";
}

#
# 0) Install gcc
#
if ( $mode == 0 || $mode == 1) {
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
	$ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64:".
	    "$ENV{LD_LIBRARY_PATH}";
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
    $ENV{LD_LIBRARY_PATH} = "$aidir/lib:$aidir/lib64:".
	"$ENV{LD_LIBRARY_PATH}";
    $ENV{GCC_EXEC_PREFIX} = "$aidir/lib/gcc-lib/";
    
#
# 1) Unpack CERN libraries
#
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
	if ( -l "$idir/lhapdf" ) {system("rm -f $idir/lhapdf");}
	system("ln -s $install{lhapdf}[0] $idir/lhapdf");
	chdir "$idir/$install{lhapdf}[0]";
	print "\nfastrun.pl: Configuring lhapdf ...\n";
#	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
	$ret = system("./configure --prefix=`pwd`");
	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF configure step, aborted!\n";}
	print "\nfastrun.pl: Making lhapdf ...\n";
# At least until LHAPDF 5.3.1 a race condition spoils multithreaded make!
#	$ret = system("make -j2");
	$ret = system("make");
	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make step, aborted!\n";}
	print "\nfastrun.pl: Make install for lhapdf ...\n";
	$ret = system("make install");
# In addition, there is a double lhapdf-config creation => ignore this error
#	if ( $ret ) {die "fastrun.pl: Error $ret in LHAPDF make install step, aborted!\n";}
	if ( $ret ) {print "fastrun.pl: Error $ret in LHAPDF make install step ignored!\n";}
	chdir "$pwdir";
    }

#
# 3) Install Nlojet++
#
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
	chdir "$idir/$install{nlojet}[0]";
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

#
# 4) Install fastNLO
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
# Set environment that is expected to have been set up!   
#
print "\nfastrun.pl: Setting environment variables for fastNLO:\n";
my $cwdir = getcwd();
print "fastrun.pl: fastNLO and gcc environment used with this installation:\n";
print "setenv CERNLIB $cwdir/cernlib\n";
print "setenv FASTNLO $cwdir/fastNLO\n";
print "setenv LHAPDF  $cwdir/lhapdf/lib\n";
print "setenv NLOJET  $cwdir/nlojet\n";
$ENV{CERNLIB}  = "$cwdir/cernlib"; 
$ENV{FASTNLO}  = "$cwdir/fastNLO";
$ENV{LHAPDF}   = "$cwdir/lhapdf/lib";
$ENV{NLOJET}   = "$cwdir/nlojet";
print "setenv PATH $cwdir/bin:$ENV{NLOJET}/bin:\${PATH}\n";
print "setenv LD_LIBRARY_PATH $cwdir/lib:$cwdir/lib64:$ENV{NLOJET}/lib:\${LD_LIBRARY_PATH}\n";
print "setenv GCC_EXEC_PREFIX $cwdir/lib/gcc-lib/\n";
print "setenv CXXFLAGS \"-O3 -I .\"\n";
$ENV{PATH}            = "$cwdir/bin:$ENV{PATH}";
$ENV{LD_LIBRARY_PATH} ="$cwdir/lib:$cwdir/lib64:".
    "$ENV{NLOJET}/lib:$ENV{LD_LIBRARY_PATH}";
$ENV{GCC_EXEC_PREFIX} ="$cwdir/lib/gcc-lib/";
$ENV{CXXFLAGS} = "-O3 -I .";
if ( $verb ) {
    my $tmp = `which gcc`;
    print "fastrun.pl: DEBUG! gcc executable used: $tmp\n";
}

#
# 5) Make fastNLO scenario
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
# 6) Run fastNLO
#
if ( $mode == 0 || $mode == 3 ) {
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Running fastNLO scenario on batch system $batch: $date\n";
    
# a) Run on grid without complete installation by source
    if ( $batch eq "GRID" && $mode == 3 ) {
# Fetching and unpacking of fastNLO binary archive
	my $file = "fastNLO-bin.tgz";
	if ( ! -f $file ) {
	    grid_storage("FETCH",$file);
	}
	if ( -f $file ) {
	    system ("tar xfz $file");
	} else {
	    die "fastrun.pl: Could not find binary tgz of fastNLO $file, aborted!\n";
	}
    }
    
# Structure change in fastNLO following change in revision 212!
    if ( $frev < 212 ) { 
	$scendir = "$ENV{FASTNLO}/author1c/hadron";
    } else {
	$scendir = "$ENV{FASTNLO}/trunk/v1.4/author1c/hadron";
    }
    chdir "$scendir" or die
	"fastrun.pl: Could not cd to dir $scendir!\n";
    
    if ( $frev < 212 ) { 
	$scendir = "$ENV{FASTNLO}/author1c/hadron";
    } else {
	$scendir = "$ENV{FASTNLO}/trunk/v1.4/author1c/hadron";
    }

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
	print "fastrun.pl: Forking command $cmd in background\n";
	system("$cmd &");
    } else {
# Run NLO calculation
	print "fastrun.pl: Running command $cmd in foreground\n";
	my $ret = system("$cmd");
	if ( $ret ) {die "fastrun.pl: Error $ret in fastNLO run step, aborted!\n";}
# Copy table to grid storage
	if ( $batch eq "GRID" ) {
	    grid_storage("SAVE","${scendir}/${tdir}/${tabnam}","$scen$ref");
	}
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: fastNLO finished: $date\n";
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
# Signal handling and grid storage: Save table on grid storage
#
# Original table naming example: fnl0002kt10_0001-hhc-born-2jet.raw
# Grid table naming example:     fnl0002kt10-hhc-born-2jet_0000.raw
#
sub grid_storage {
    my $signam = shift;
    my $trfile = shift;
    my $tabdir = shift;
    $trfile =~ s/\/\.\//\//g;
    if ( $tabdir ) {
	$tabdir =~ s/\/\.\//\//g;
    } else {
	$tabdir = "";
    }
    
    my $sename = "ekp-lcg-se.physik.uni-karlsruhe.de";
    my $sepath = "/wlcg/data/users/cms/rabbertz/";
    my @files = ( "$trfile" );
    
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    my $gjobnr = $ENV{MY_JOB};
    print "\nfastrun.pl: Received signal $signam: $date\n";
    unless ( $signam eq "FETCH" ) {
	print "fastrun.pl: Saving table $trfile for\n";
	print "fastrun.pl: fastNLO job no. $jobnr and\n";
	print "fastrun.pl: grid-control job no. $gjobnr\n";
    } else {
	print "fastrun.pl: Fetching binary fastNLO archive $trfile from\n";
	print "fastrun.pl: SE $sename in path $sepath\n";
    }	

    foreach my $file ( @files ) {
	unless ( $signam eq "FETCH" ) {
	    print "fastrun.pl: Trying to find table file $file ...\n";
	    system("pwd");
	    system("ls -la");
	    if ( -f "${file}" ) {
# Create target directory if necessary
		my $gcmd = "edg-gridftp-mkdir ".
		    "gsiftp://${sename}/${sepath}".
		    "fastNLO_tables/${tabdir}";
		print "Command $gcmd\n";
		my $ret = system("$gcmd");
		if ( $ret ) {print "fastrun.pl: WARNING! Creation of grid storage directory failed: $ret!\n";}
# Change job numbering according to grid-control
		my $newnum = substr("0000$gjobnr",-4);	    
		my @tmp = split("/",$file);
		my $tabnam = pop(@tmp);
		$tabnam =~ s/_0001//;
		$tabnam =~ s/\.raw/_${newnum}\.raw/;
		$gcmd = "globus-url-copy ".
		    "file://${file} ".
		    "gsiftp://${sename}/${sepath}".
		    "fastNLO_tables/${tabdir}/${tabnam}";
		print "Command $gcmd\n";
		system("$gcmd");
	    } else {
		system("pwd");
		system("ls -la");
		die "fastrun.pl: Could not find table file $file\n";
	    }
	} else {
	    print "fastrun.pl: Trying to fetch binary fastNLO archive $file on\n";
	    print "fastrun.pl: SE $sename in path $sepath\n";
	    my $gcmd = "globus-url-copy gsiftp://${sename}/${sepath}". 
		"fastNLO_archives/${file} file://`pwd`/${file}";
	    print "Command $gcmd\n";
	    system("$gcmd");
	    system("pwd");
	    system("ls -la");
	}
    }
    return 0;
}
