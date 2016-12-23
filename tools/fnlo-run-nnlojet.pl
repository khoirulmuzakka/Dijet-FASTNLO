#!/usr/bin/env perl
#
# fastNLO run script for NNLOJet
# Version:
#
# created by K. Rabbertz: 13.12.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use File::Basename;
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
    open STDOUT, "| tee fnlo-run-nlojet_${gjobnr}.log" or die
        "fnlo-run-nnlojet.pl: ERROR! Can't tee STDOUT.\n";
    open STDERR, "| tee fnlo-run-nlojet_${gjobnr}.err" or die
        "fnlo-run-nnlojet.pl: ERROR! Can't tee STDERR.\n";
}

my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-run-nnlojet.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################\n";
print "# fnlo-run-nnlojet.pl: Starting run of fastNLO: FASTRUN0_$date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h ) =
    ( "", "" );
getopts('dh') or die "fnlo-run-nnlojet.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-run-nnlojet.pl\n";
    print "Usage: fnlo-run-nnlojet.pl [switches/options] (ScenarioType) (ScenarioOrder) (ScenarioMode)\n";
    print "  -d debug        Switch debug/verbose mode on\n";
    print "  -h              Print this text\n";
    print "\n";
    print "Example:\n";
    print "Run NNLOJET scenario:\n";
    print "   ./fnlo-run-nnlojet.pl ZJ LO fastwarmup\n\n";
    exit;
}
my $verb = 0;
if ( $opt_d ) {
    $verb = 1;
    print "fnlo-run-nnlojet.pl: Debug/verbose mode is active.\n";
}

#
# Parse arguments
#
my $scentype = "def";
my $scenord  = "def";
my $scenmode = "def";
unless ( @ARGV == 3 ) {
    die "fnlo-run-nnlojet.pl: Error! Need three argument of scenario type, order, and mode!\n";
}
if ( @ARGV > 0 ) {
    $scentype = shift;
    $scenord  = shift;
    $scenmode = shift;
}
print "fnlo-run-nnlojet.pl: The scenario type, order, and mode are ${scentype}, ${scenord}, and ${scenmode}.\n";

#
# Print system info (some parts in debug mode only; df commands can get stuck ...)
#
print "\n######################################################\n";
print "fnlo-run-nnlojet.pl: System information for debugging purposes:\n";
print "######################################################\n";
my $host = `hostname`;
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"hostname\" command failed.\n\n";
} else {
    chomp $host;
    print "fnlo-run-nnlojet.pl: The system's hostname is (hostname):\n$host\n\n";
}
my $osvers = `uname -a`;
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"uname -a\" command failed.\n\n";
} else {
    print "fnlo-run-nnlojet.pl: Your operating system is (uname -a):\n$osvers\n";
}
my $procvers = `cat /proc/version`;
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"cat /proc/version\" command failed.\n\n";
} else {
    print "fnlo-run-nnlojet.pl: Your linux version is (/proc/version):\n$procvers\n";
}
my $freemem = `free`;
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"free\" command failed.\n\n";
} else {
    print "fnlo-run-nnlojet.pl: The available memory is:\n$freemem\n";
}
my $cpumod = `cat /proc/cpuinfo | grep \"model name\"`;
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"cat /proc/cpuinfo\" command failed.\n\n";
} else {
    my $cpufrq = `cat /proc/cpuinfo | grep \"cpu MHz\"`;
    print "fnlo-run-nnlojet.pl: The processor type is:\n${cpumod}at\n${cpufrq}\n";
}
if ( $verb ) {
    my $freedisk = `df -h`;
    if ( $? ) {
        print "fnlo-run-nnlojet.pl: Info: \"df -h\" command failed.\n\n";
    } else {
        print "fnlo-run-nnlojet.pl: The available disk space is:\n$freedisk\n";
    }
    my $freenode = `df -hi`;
    if ( $? ) {
        print "fnlo-run-nnlojet.pl: Info: \"df -hi\" command failed.\n\n";
    } else {
        print "fnlo-run-nnlojet.pl: The available inode space is:\n$freenode\n";
    }
}
my $cwd = getcwd();
if ( $? ) {
    print "fnlo-run-nnlojet.pl: Info: \"getcwd()\" command failed.\n\n";
} else {
    print "fnlo-run-nnlojet.pl: The current working directory is:\n$cwd\n\n";
    print "fnlo-run-nnlojet.pl: The current working directory's content is:\n";
    my $ret = system("ls -laR");
    if ( $ret ) {print "fnlo-run-nnlojet.pl: Couldn't list current directory: $ret, skipped!\n";}
}
print "######################################################\n\n";

#
# Set system paths environment
#
my $mybase="/cvmfs/etp.kit.edu/nnlo";
my $gccbase="/cvmfs/cms.cern.ch/slc6_amd64_gcc481";
if ( $ENV{PATH} ) {
    $ENV{PATH} = "${mybase}/bin:$ENV{PATH}";
} else {
    $ENV{PATH} = "${mybase}/bin";
}
if ( $ENV{LD_LIBRARY_PATH} ) {
    $ENV{LD_LIBRARY_PATH} = "${mybase}/lib:${mybase}/lib/root:$ENV{LD_LIBRARY_PATH}";
} else {
    $ENV{LD_LIBRARY_PATH} = "${mybase}/lib:${mybase}/lib/root";
}

#
# Print some location info
#
print "\n######################################################\n";
print "fnlo-run-nnlojet.pl: List NNLOJET process directory for debugging purposes:\n";
print "######################################################\n";
my $ret = system("ls -la ${mybase}/src/NNLOJET-rev3278/driver/process");
if ( $ret ) {die "fnlo-run-nnlojet.pl: ERROR! Couldn't list NNLOJET process directory: $ret, aborted!\n";}

#
# Run fastNLO
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
my $scenname = "${scentype}.${scenord}-${scenmode}";
print "\nfnlo-run-nnlojet.pl: Running fastNLO scenario name ${scenname}: $date\n";

# Run NLO calculation
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nnlojet.pl: Starting calculation: FASTCAL0_$date\n";
my $cmd = "${mybase}/src/NNLOJET-rev3278/driver/NNLOJET -run ${scenname}.run |& tee ${scenname}.log";
print "\nfnlo-run-nnlojet.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "fnlo-run-nnlojet.pl: ERROR! Error $ret in fastNLO run step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nnlojet.pl: Calculation finished: FASTCAL1_$date\n";

# Copy/rename results for storage
print "\nfnlo-run-nnlojet.pl: The final working directory's content is:\n";
$ret = system("ls -laR");
if ( $ret ) {print "fnlo-run-nnlojet.pl: Couldn't list current directory: $ret, skipped!\n";}
my $defnam = "ZJ-LO-";
#my @files = glob("${scentype}*.dat ${scentype}*.log ${scentype}*.run *.root *.tab *.wrm");
my @files = glob("${scentype}*.log *.root *.tab *.wrm");
foreach my $file ( @files ) {
    (my $name, my $dir, my $ext) = fileparse($file,'\.[^\.]*$');
    $ext =~ s/^\.//;
    my $newname = "${name}_${gjobnr}.${ext}";
    $newname =~ s/^zj-//;
# APPLgrid special
    if ( ($ext =~ m/root/) ) {
#&& ($name =~ m/^$defnam/) && !($scenname =~ m/$defnam/) ) { 
	$newname =~ s/^$defnam/${scentype}\.${scenord}\./;
    }
    rename $file, $newname;
    print "fnlo-run-nnlojet.pl: $file has been renamed to ${newname}\n";
}

$cmd = "tar cfz ${scenname}_${gjobnr}.tar.gz ZJ*";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-run-nnlojet.pl: ERROR! Result files could not be archived to ${scenname}_${gjobnr}.tar.gz, aborted!\n";}
# Copy to standard name as gc sandbox output file
$ret = system("cp -p ${scenname}_${gjobnr}.tar.gz ${scenname}.tar.gz");
if ( $ret ) {print "fnlo-run-nnlojet.pl: Couldn't copy output archive, export in gc output files will fail: $ret!\n";}
$ret = system("ls -laR");
if ( $ret ) {print "fnlo-run-nnlojet.pl: Couldn't list current directory: $ret, skipped!\n";}

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n###############################\n";
print "# fnlo-run-nnlojet.pl: fastNLO finished: FASTRUN1_$date\n";
print "###############################\n\n";
exit 0;
