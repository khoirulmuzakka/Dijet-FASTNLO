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
    ( "LCG", ".", "0", "187", "", "0001",
      "0", "LO", "CTEQ", "", ".", "", "" );
getopts('b:d:e:f:hj:m:o:p:rs:t:v') or die "fastrun.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastrun.pl\n";
    print "Usage: fastrun.pl [switches/options] ([scenario])\n";
    print "  -b batch        Batch system used: LCG (def.) or PBS\n";
    print "  -d dir          Installation directory (def.=.)\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -f rev          fastNLO revision to use (def.=187)\n";
    print "  -h              Print this text\n";
    print "  -j jobnr        Job number to attach\n";
    print "  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -p pdf          CTEQ parton densities (def.) or LHAPDF\n";
    print "  -r              Reference calculation incl. pdf access\n";
    print "  -s dir          Archive source directory (def.=.)\n";
    print "  -t dir          Output target directory: ".
	"(def.= {scen}{ref}_{jobnr} with\n                  ".
	"ref. to working directory in fastNLO installation)\n";
    print "  -v              Switch verbose mode on\n\n";
    exit;
}

unless ( $opt_b eq "LCG" || $opt_b eq "PBS" ) {
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
my $sdir  = $opt_s;
my $verb  = "";
print "fastrun.pl: Job mode is $mode\n";
print "fastrun.pl: Using fastNLO revision $frev\n";
print "fastrun.pl: Directory for/of installation is $idir\n";
unless ( $mode == 1 || $mode == 2 ) {
    print "fastrun.pl: Running on batch system $batch\n";
    print "fastrun.pl: Attaching job number $jobnr\n";
    print "fastrun.pl: Running in order $order\n";
    print "fastrun.pl: Running with pdf $pdf\n";
    print "fastrun.pl: Maximal event number: $nmax\n";
}
if ( $opt_r ) {
    $ref = "ref";
    unless ( $mode == 2 ) {
	print "fastrun.pl: Making/running in reference mode.\n";
    }
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
# Set signal traps to store output at end of allocated batch time
#
# LCG only for the moment
if ( $batch eq "LCG" ) { 
    $SIG{INT}  = \&grid_store;
    $SIG{TERM} = \&grid_store;
}

#
# Initialization
#
# Global constants
# Switches for scenario cc file
our $refsw0 = " iref = 0"; 
our $refsw1 = " iref = 1"; 
# Makefile switches
our $maksw0 = "-o $scen ";
our $maksw1 = "-o ${scen}ref ";

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

# Global vars for output storage on signal catch, jobnr not yet attached!
our $tabdir  = "${scen}${ref}";
our $tabnam0 = "${scen}${ref}-hhc-$runmode{$order}[0]-2jet";
our $tabnam1 = "${tabnam0}.raw";
our $lognam  = "${scen}${ref}_${order}.log";

# Directories
my $pwdir = getcwd();
chdir $idir;
my $aidir = getcwd();
chdir $pwdir;

# Define install hash
my %install;
# First entry (index 0!): Subdirecory name into which the archive is unpacked!
$install{cernlib}[0]    = "cernlib-2003";
$install{lhapdf}[0]     = "lhapdf-5.3.0";
$install{nlojet}[0]     = "nlojet++-2.0.1";
$install{nlojetfix}[0]  = "nlojet++-2.0.1";
$install{fastNLO}[0]    = "fastNLO-rev${frev}";
# Second: Archive filenames 
foreach my $comp ( keys %install ) {
    my $tmp = $install{$comp}[0].".tar.gz";
    if ( $comp eq "nlojetfix" ) {
	$install{$comp}[1] = "nlojet++-2.0.1-fix.tar.gz";
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
# 1) Unpack CERN libraries
#
if ( $mode == 0 || $mode == 1) {
    unless ( -e "$idir/$install{cernlib}[0]" ) {
	$date = `date +%d%m%Y_%H%M%S`;
	chomp $date;
	print "\nfastrun.pl: Unpacking CERN libraries in $install{cernlib}[1]: $date\n";
	system("tar xz -C $idir -f $sdir/$install{cernlib}[1]");
#    system("rm -f $install{cernlib}[1]");
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
	system("tar xz -C $idir -f $sdir/$install{lhapdf}[1]");
#    system("rm -f $install{lhapdf}[1]");
	system("ln -s $install{lhapdf}[0] $idir/lhapdf");
	chdir "$idir/$install{lhapdf}[0]";
	print "\nfastrun.pl: Configuring lhapdf ...\n";
#	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
	system("./configure --prefix=`pwd`");
	print "\nfastrun.pl: Making lhapdf ...\n";
	system("make");
	print "\nfastrun.pl: Make install for lhapdf ...\n";
	system("make install");
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
	system("tar xz -C $idir -f $sdir/$install{nlojet}[1]");
#    system("rm -f $install{nlojet}[1]");
	print "\nfastrun.pl: Unpacking fix for $install{nlojet}[1] ...\n";
	system("tar xzv -C $idir -f $sdir/$install{nlojetfix}[1]");
#    system("rm -f $install{nlojetfix}[1]");
	system("ln -s  $install{nlojet}[0] $idir/nlojet");
	chdir "$idir/$install{nlojet}[0]";
	print "\nfastrun.pl: Configuring Nlojet++ ...\n";
#	system("./configure --prefix=`pwd` --exec-prefix=$aidir");
	system("./configure --prefix=`pwd`");
	print "\nfastrun.pl: Making Nlojet++ ...\n";
	system("make CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
	print "\nfastrun.pl: Make install for Nlojet++ ...\n";
	system("make install CFLAGS=\"-O3 -Wall\" CXXFLAGS=\"-O3 -Wall\"");
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
	system("tar xz -C $idir -f $sdir/$install{fastNLO}[1]");
#    system("rm -f $install{fastNLO}[1]");
    }
}

#
# 5) Make fastNLO scenario
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastrun.pl: Setting environment for fastNLO: $date\n";
chdir $idir;
my $cwdir = getcwd();
$ENV{CERNLIB}  = "$cwdir/$install{cernlib}[0]"; 
$ENV{LHAPDF}   = "$cwdir/$install{lhapdf}[0]/lib";
$ENV{NLOJET}   = "$cwdir/$install{nlojet}[0]";
$ENV{fastNLO}  = "$cwdir/$install{fastNLO}[0]";
$ENV{CXXFLAGS} = "-O3 -I .";
print "CERNLIB: $ENV{CERNLIB}\n";
print "LHAPDF: $ENV{LHAPDF}\n";
print "NLOJET: $ENV{NLOJET}\n";
print "fastNLO: $ENV{fastNLO}\n";
print "CXXFLAGS: $ENV{CXXFLAGS}\n";
# Structure change in fastNLO following change in revision 212!
my $scendir;
if ( $frev < 212 ) { 
    $scendir = "$ENV{fastNLO}/author1c/hadron";
} else {
    $scendir = "$ENV{fastNLO}/trunk/v1.4/author1c/hadron";
}
chdir "$scendir" or die
    "fastrun.pl: Could not cd to dir $scendir!\n";
    
if ( $mode == 0 || $mode == 2 ) {
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
	print "\nfastrun.pl: Making scenario $scen of fastNLO: $date\n";
	chdir $scendir;
	system("make $scen");
    }
}
    
#
# 6) Run fastNLO
#
chdir "$scendir" or die
    "fastrun.pl: Could not cd to dir $scendir!\n";
if ( $mode == 0 || $mode == 3 ) {
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
    if ( $batch eq "LCG" ) {
# Fork NLO calculation
	print "fastrun.pl: Forking command $cmd in background\n";
	system("$cmd &");
    } else {
# Run NLO calculation
	print "fastrun.pl: Running command $cmd in foreground\n";
	system("$cmd");
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
	if ( $batch eq "LCG" ) {
	    sleep 1000;
	    $date = `date +%d%m%Y_%H%M%S`;
	    print "\nfastrun.pl: Saving table: $date\n";
	    chomp $date;
	    grid_store("SAVE");
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
# Signal handling: Save table and log file on grid storage
#
sub grid_store {
    my $signam = shift;
    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfastrun.pl: Received signal $signam: $date\n";
    print "fastrun.pl: Saving table and log files for job nr. $jobnr\n";

    my $gcmd = "globus-url-copy file:`pwd`/${lognam}_${jobnr} ".
	"gsiftp://ekp-lcg-se.physik.uni-karlsruhe.de/".
	"grid/users/cms/cmssgm/${lognam}_${jobnr}";
    print "Command $gcmd\n";
    system("$gcmd");
    my @files = ( "$tabnam0", "$tabnam1" );
    foreach my $file ( @files ) {
#	my $lcmd = "lcg-cr --vo cms ".
#	    "-l lfn:${file}_${jobnr} file:`pwd`/${file}"; 
	if ( -f "${tabdir}/${file}" ) {
	    my $gcmd = "globus-url-copy file:`pwd`/${tabdir}/${file} ".
		"gsiftp://ekp-lcg-se.physik.uni-karlsruhe.de/".
		"grid/users/cms/cmssgm/${file}_${jobnr}";
	    print "Command $gcmd\n";
	    system("$gcmd");
	}
    }
    return 0;
}
