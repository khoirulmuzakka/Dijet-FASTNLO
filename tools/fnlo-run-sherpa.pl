#!/usr/bin/env perl
#
# fastNLO run script for SHERPA
# Version:
#
# created by K. Rabbertz: 29.01.2006
# adapted by K. Rabbertz from fnlo-run-nlojet.pl: 06.10.2015
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
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
    $gjobnr = substr("0000$gjobnr",-4);
    open STDOUT, "| tee fnlo-run-sherpa_${gjobnr}.log" or die
        "fnlo-run-sherpa.pl: ERROR! Can't tee STDOUT.\n";
    open STDERR, "| tee fnlo-run-sherpa_${gjobnr}.err" or die
        "fnlo-run-sherpa.pl: ERROR! Can't tee STDERR.\n";
}

my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-run-sherpa.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################\n";
print "# fnlo-run-sherpa.pl: Starting run of Sherpa with fastNLO: FASTRUN0_$date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_h, $opt_j,
      $opt_t, $opt_v, $opt_w, $opt_x ) =
    ( "GC", "", "0", "", "0001",
      "", "2.3", "", "" );
getopts('b:de:hj:t:v:wx:') or die "fnlo-run-sherpa.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-run-sherpa.pl\n";
    print "Usage: fnlo-run-sherpa.pl [switches/options] ([ScenarioName])\n";
    print "  -b batch        Batch system used: GC (def.),\n";
    print "                  GC (grid-control), or LOCAL\n";
    print "  -d debug        Switch debug/verbose mode on\n";
    print "  -e nevents      Number of events (def.=from Run_prod.dat file)\n";
    print "  -h              Print this text\n";
    print "  -j jobnr        Job number to attach (def.=0001)\n";
    print "  -t dir          Output target directory: ".
        "(def.= mcgrid/MCgrid_\${scenname} with ref. to working directory)\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n";
    print "  -w              Warm-up run to determine x limits (def.=F; version 2.3+ only)\n";
    print "  -x seed         Seed for random number generator (def.=F; version 2.3+ only)\n";
    print "\n";
    print "Example:\n";
    print "Run SHERPA scenario:\n";
    print "   ./fnlo-run-sherpa.pl [-b GC|LOCAL] [-e nevents] [-x seed] ATLAS_2012_I1082936\n\n";
    exit;
}

unless ( $opt_b eq "LOCAL" || $opt_b eq "GC" || $opt_b eq "GRID" || $opt_b eq "PBS" ) {
    die "fnlo-run-sherpa.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( $opt_b eq "GC" ) {
    die "fnlo-run-sherpa.pl: Error! Batch system other than GC (grid-control) not updated/supported: $opt_b, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fnlo-run-sherpa.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( $opt_j =~ m/\d{4}/ && $opt_j !~ m/\D+/ ) {
    die "fnlo-run-sherpa.pl: Error! Illegal job number (nnnn): $opt_j, aborted.\n";
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-run-sherpa.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}

my $batch = $opt_b;
my $nmax  = $opt_e;
my $jobnr = $opt_j;
my $verb  = "";
my $vers  = $opt_v;
my $wrm   = "";
if ( $opt_w ) { $wrm = "wrm";}
my $seed  = "";
if ( $opt_x ) { $seed = $opt_x ;}
print "fnlo-run-sherpa.pl: Running on batch system $batch\n";
print "fnlo-run-sherpa.pl: Maximal event number: $nmax\n";
print "fnlo-run-sherpa.pl: Attaching job number $jobnr\n";
if ( $wrm eq "wrm" ) {
    print "fnlo-run-sherpa.pl: Running in warm-up mode\n";
}
if ( $opt_d ) {
    $verb = 1;
    print "fnlo-run-sherpa.pl: Debug/verbose mode is active.\n";
}

#
# Parse arguments
#
my $scenname = "def";
unless ( @ARGV == 1 ) {
    die "fnlo-run-sherpa.pl: Error! Need one scenario description in the form of ScenarioName=RivetAnalysisName!\n";
}
if ( @ARGV > 0 ) {
    $scenname = shift;
    chomp $scenname;
}
print "fnlo-run-sherpa.pl: Running scenario ${scenname}.\n";
my $tdir = "mcgrid/MCgrid_${scenname}";
if ( $opt_t ) {
    $tdir  = $opt_t;
}
print "fnlo-run-sherpa.pl: Output target directory $tdir\n";

#
# Initialization
#

# Directories
my $rundir = getcwd();
chomp $rundir;

#
# Print system info (some parts in debug mode only; df commands can get stuck ...)
#
print "\n######################################################\n";
print "fnlo-run-sherpa.pl: System information for debugging purposes:\n";
print "######################################################\n";
my $host = `hostname`;
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"hostname\" command failed.\n\n";
} else {
    chomp $host;
    print "fnlo-run-sherpa.pl: The system's hostname is (hostname):\n$host\n\n";
}
my $osvers = `uname -a`;
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"uname -a\" command failed.\n\n";
} else {
    print "fnlo-run-sherpa.pl: Your operating system is (uname -a):\n$osvers\n";
}
my $procvers = `cat /proc/version`;
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"cat /proc/version\" command failed.\n\n";
} else {
    print "fnlo-run-sherpa.pl: Your linux version is (/proc/version):\n$procvers\n";
}
my $freemem = `free`;
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"free\" command failed.\n\n";
} else {
    print "fnlo-run-sherpa.pl: The available memory is:\n$freemem\n";
}
my $cpumod = `cat /proc/cpuinfo | grep \"model name\"`;
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"cat /proc/cpuinfo\" command failed.\n\n";
} else {
    my $cpufrq = `cat /proc/cpuinfo | grep \"cpu MHz\"`;
    print "fnlo-run-sherpa.pl: The processor type is:\n${cpumod}at\n${cpufrq}\n";
}
if ( $verb ) {
    my $freedisk = `df -h`;
    if ( $? ) {
        print "fnlo-run-sherpa.pl: Info: \"df -h\" command failed.\n\n";
    } else {
        print "fnlo-run-sherpa.pl: The available disk space is:\n$freedisk\n";
    }
    my $freenode = `df -hi`;
    if ( $? ) {
        print "fnlo-run-sherpa.pl: Info: \"df -hi\" command failed.\n\n";
    } else {
        print "fnlo-run-sherpa.pl: The available inode space is:\n$freenode\n";
    }
}
my $cwd = getcwd();
if ( $? ) {
    print "fnlo-run-sherpa.pl: Info: \"getcwd()\" command failed.\n\n";
} else {
    print "fnlo-run-sherpa.pl: The current working directory is:\n$cwd\n\n";
    print "fnlo-run-sherpa.pl: The current working directory's content is:\n";
    my $ret = system("ls -laR");
    if ( $ret ) {print "fnlo-run-sherpa.pl: Couldn't list current directory: $ret, skipped!\n";}
}

#
# Set system paths environment
#
$ENV{FASTNLO} = "$ENV{HOME}/local";
if ( $ENV{PATH} ) {
    $ENV{PATH} .= ":$ENV{FASTNLO}/bin";
} else {
    $ENV{PATH}  = "$ENV{FASTNLO}/bin";
}
if ( $ENV{LD_LIBRARY_PATH} ) {
    $ENV{LD_LIBRARY_PATH} .= ":$ENV{FASTNLO}/lib";
} else {
    $ENV{LD_LIBRARY_PATH}  = "$ENV{FASTNLO}/lib";
}
if ( $ENV{PYTHONPATH} ) {
    $ENV{PYTHONPATH} .= ":$ENV{HOME}/local/lib64/python2.6/site-packages:$ENV{HOME}/local/lib/python2.6/site-packages";
} else {
    $ENV{PYTHONPATH}  = "$ENV{HOME}/local/lib64/python2.6/site-packages:$ENV{HOME}/local/lib/python2.6/site-packages";
}
$ENV{RIVET_ANALYSIS_PATH} = "$ENV{HOME}/local/share/Rivet:$ENV{HOME}/local/share/fastnlo_interface_nlojet/RivetAdditions:$ENV{HOME}/local/share/mcgrid-2.0-examples/${scenname}/analysis";
$ENV{RIVET_INFO_PATH}     = "$ENV{HOME}/local/share/Rivet:$ENV{HOME}/local/share/fastnlo_interface_nlojet/RivetAdditions";
$ENV{RIVET_REF_PATH}      = "$ENV{HOME}/local/share/Rivet:$ENV{HOME}/local/share/fastnlo_interface_nlojet/RivetAdditions";
$ENV{PKG_CONFIG_PATH}     = "$ENV{HOME}/local/lib/pkgconfig";
print "\nfnlo-run-sherpa.pl: The environment settings are:\n\n";
my $ret = system("printenv");
if ( $ret ) {print "fnlo-run-sherpa.pl: Couldn't print environment settings: $ret, skipped!\n";}
print "######################################################\n\n";

#
# Run fastNLO
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-sherpa.pl: Running fastNLO version $vers scenario on batch system $batch: $date\n";

# In case of non-local running need to unpack binary fastNLO distribution
if ( $batch ne "LOCAL" ) {

# Fetching and unpacking of fastNLO binary archive
    my $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-sherpa.pl: Get/unpack MCgrid/Rivet analysis package: BINGET0_$date\n";

    my $file = "${scenname}.tgz";
    if ( ! -f $file ) {
        if ( $batch eq "GC" ) {
            die "fnlo-run-sherpa.pl: ERROR! Could not find $file in working dir, aborted!\nCheck your grid storage options for grid-control!\n";
        }
    }
    if ( -f $file ) {
        system ("tar xfz $file");
    } else {
        die "fnlo-run-sherpa.pl: ERROR! Could not find $file, aborted!\n";
    }

    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-sherpa.pl: Unpacked MCgrid/Rivet analysis package: BINGET1_$date\n";
}






my $scendir = "$rundir/${scenname}";
chdir "$scendir" or die
    "fnlo-run-sherpa.pl: ERROR! Could not cd to dir $scendir!\n";

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-sherpa.pl: Running fastNLO: $date\n";
my $cmd;
if ( $vers eq "2.3" ) {
    $cmd = "$ENV{FASTNLO}/bin/Sherpa -f Run_prod.dat ";
    if ( $nmax ) { $cmd .= "-e $nmax "; }
    if ( $seed ) { $cmd .= "-R $seed "; }
} else {
    die "fnlo-run-sherpa.pl: ERROR! Unsupported fastNLO version $vers requested, aborted!\n";
}

# Run NLO calculation
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-sherpa.pl: Starting calculation: FASTCAL0_$date\n";
print "\nfnlo-run-sherpa.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "fnlo-run-sherpa.pl: ERROR! Error $ret in fastNLO run step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-sherpa.pl: Calculation finished: FASTCAL1_$date\n";

# Copy/rename results for grid storage
print "\nfnlo-run-sherpa.pl: The final working directory's content is:\n";
$ret = system("ls -laR");
if ( $ret ) {print "fnlo-run-sherpa.pl: Couldn't list current directory: $ret, skipped!\n";}
my $refext   = "yoda";
my $tabext   = "tab";
my @reffiles = `find . -follow -name \*.${refext} -printf \"%f\\n\"`;
chomp @reffiles;
my @tabfiles = `find . -follow -name \*.${tabext} -printf \"%f\\n\"`;
chomp @tabfiles;
my $refpath  = ".";
my $tabpath  = "mcgrid/MCgrid_${scenname}";
print "\nfnlo-run-sherpa.pl: Sherpa product files @reffiles @tabfiles\n";
print "\nfnlo-run-sherpa.pl: Source paths: $refpath $tabpath\n";
my $trefpath = "${scenname}ref";
my $ttabpath = "${scenname}";
print "fnlo-run-sherpa.pl: Target paths: $trefpath $ttabpath\n";
my @treffiles = @reffiles;
foreach my $file ( @treffiles ) {
# Eliminate ending
    $file =~ s/\.yoda//;
# Add scenario name to front and job number to the end
    $file = "${scenname}_${gjobnr}";
# Add back ending
    $file .= ".yoda";
}
my @ttabfiles = @tabfiles;
foreach my $file ( @ttabfiles ) {
# Eliminate ending
    $file =~ s/\.${tabext}//;
# Add scenario name to front and job number to the end
    $file = "${scenname}_${file}_${gjobnr}";
# Add back ending
    $file .= ".${tabext}";
}

print "fnlo-run-sherpa.pl: Target files: @treffiles @ttabfiles\n";

if ( $batch ne "LOCAL" ) {
# Copy/rename for grid storage via GC
    if ( $batch eq "GC" ) {
        print "fnlo-run-sherpa.pl: Info: Batch mode $batch: Grid storage done by grid-control.\n";
# Table files	
	foreach my $i (0 .. $#tabfiles) {
	    my $spath = $tabpath;
	    my $sfile = $tabfiles[$i];
	    my $tfile = $ttabfiles[$i];
	    if ( ! -f "$rundir/$sfile" ) {
		print "fnlo-run-sherpa.pl: Copy product to current directory for storage.\n";
		my $ret = system("cp -p $spath/$sfile $rundir/$sfile");
		if ( $ret ) {die "fnlo-run-sherpa.pl: Couldn't copy product into ".
				 "current directory $rundir: $ret, aborted!\n";}
	    }
	    print "fnlo-run-sherpa.pl: Rename product $sfile to expected file name $tfile for storage.\n";
	    my $ret = system("mv -f $rundir/$sfile $rundir/$tfile");
	    if ( $ret ) {die "fnlo-run-sherpa.pl: Couldn't rename product ${sfile} into ".
			     "${tfile}: $ret, aborted!\n";}

	}
# Reference files (yoda)	
	foreach my $i (0 .. $#reffiles) {
	    my $spath = $refpath;
	    my $sfile = $reffiles[$i];
	    my $tfile = $treffiles[$i];
	    if ( ! -f "$rundir/$sfile" ) {
		print "fnlo-run-sherpa.pl: Copy product to current directory for storage.\n";
		my $ret = system("cp -p $spath/$sfile $rundir/$sfile");
		if ( $ret ) {die "fnlo-run-sherpa.pl: Couldn't copy product into ".
				 "current directory $rundir: $ret, aborted!\n";}
	    }
	    print "fnlo-run-sherpa.pl: Rename product $sfile to expected file name $tfile for storage.\n";
	    my $ret = system("mv -f $rundir/$sfile $rundir/$tfile");
	    if ( $ret ) {die "fnlo-run-sherpa.pl: Couldn't rename product ${sfile} into ".
			     "${tfile}: $ret, aborted!\n";}
	}
# Also rename related log and err files in cwd for grid storage, if file size larger than zero
# Remove product extension to rename log/err files
        chdir $rundir or die "fnlo-run-nlojet.pl: ERROR! Couldn't cd to $rundir, aborted!\n";
	my $tfile = "${scenname}_${gjobnr}";
	print "fnlo-run-sherpa.pl: Copy log files to correct file names for storage.\n\n";
	if ( -f "job.stderr" && ! -z "job.stderr" ) {
	    my $ret = system("cp -p job.stderr ${tfile}.err");
	    if ( $ret ) {die "fnlo-run-sherpa.pl: ERROR! Couldn't copy job.stderr ".
			     "to ${tfile}.err: $ret, aborted!\n";}
	}
	if ( -f "job.stdout" && ! -z "job.stdout" ) {
	    $ret = system("cp -p job.stdout ${tfile}.log");
	    if ( $ret ) {die "fnlo-run-sherpa.pl: ERROR! Couldn't copy job.stdout ".
			     "to ${tfile}.log: $ret, aborted!\n";}
	}

	print "\n";
	$ret = system("pwd");
	if ( $ret ) {print "fnlo-run-sherpa.pl: WARNING! Couldn't print cwd!\n";}
	$ret = system("ls -la");
	if ( $ret ) {print "fnlo-run-sherpa.pl: WARNING! Couldn't list cwd!\n";}
    }
    my $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-sherpa.pl: Results stored: TABSAV1_$date\n";
}

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n###############################\n";
print "# fnlo-run-sherpa.pl: fastNLO finished: FASTRUN1_$date\n";
print "###############################\n\n";
exit 0;
