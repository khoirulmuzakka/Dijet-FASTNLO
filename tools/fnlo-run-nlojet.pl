#!/usr/bin/env perl
#
# fastNLO run script for NLOJet++
# Version:
#
# created by K. Rabbertz: 29.01.2006
# adapted by K. Rabbertz from fastrun.pl: 28.06.2014
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
    open STDOUT, "| tee fnlo-run-nlojet_${gjobnr}.log" or die
        "fnlo-run-nlojet.pl: ERROR! Can't tee STDOUT.\n";
    open STDERR, "| tee fnlo-run-nlojet_${gjobnr}.err" or die
        "fnlo-run-nlojet.pl: ERROR! Can't tee STDERR.\n";
}

my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-run-nlojet.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n######################################\n";
print "# fnlo-run-nlojet.pl: Starting run of fastNLO: FASTRUN0_$date\n";
print "######################################\n\n";

#
# Parse options
#
our ( $opt_b, $opt_d, $opt_e, $opt_g, $opt_h, $opt_j,
      $opt_n, $opt_o, $opt_r, $opt_t, $opt_v, $opt_w, $opt_x ) =
    ( "GC", "", "0", "guc", "", "0001",
      "2jet", "LO", "", "", "2.3", "", "" );
getopts('b:de:g:hj:n:o:rt:v:wx:') or die "fnlo-run-nlojet.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-run-nlojet.pl\n";
    print "Usage: fnlo-run-nlojet.pl [switches/options] ([ScenarioType_ScenarioName])\n";
    print "  -b batch        Batch system used: GC (def.),\n";
    print "                  GC (grid-control), LOCAL, GRID, or PBS\n";
    print "  -d debug        Switch debug/verbose mode on\n";
    print "  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)\n";
    print "  -g prot         Grid storage protocol to use: guc (def.), srm\n";
    print "                  Needed only, when the batch system is neither LOCAL nor GC\n";
    print "  -h              Print this text\n";
    print "  -j jobnr        Job number to attach (def.=0001)\n";
    print "  -n njet         NLOJet++ jet mode: 2jet (def.) for 2+, 3jet for 3+ final states\n";
    print "  -o order        LO (def.) or NLO calculation\n";
    print "  -r              Reference calculation incl. pdf access (currently not available!)\n";
    print "  -t dir          Output target directory: ".
        "(def.= {scen}{ref}_{jobnr} with\n                  ".
        "ref. to working directory in fastNLO installation)\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n";
    print "  -w              Warm-up run to determine x limits (def.=F; version 2.3+ only)\n";
    print "  -x seed         Seed for random number generator (def.=F; version 2.3+ only)\n";
    print "\n";
    print "Example:\n";
    print "Run NLOJet++ scenario (to run scenario in reference mode use option -r):\n";
    print "   ./fnlo-run-nlojet.pl [-b GC|LOCAL|GRID|PBS] [-e max-events] [-t ./{scen}{ref}_{jobnr}|tdir] [-v 2.3] InclusiveJets_Example_v23_fix\n\n";
    exit;
}

unless ( $opt_b eq "LOCAL" || $opt_b eq "GC" || $opt_b eq "GRID" || $opt_b eq "PBS" ) {
    die "fnlo-run-nlojet.pl: Error! Illegal batch system: $opt_b, aborted.\n";
}
unless ( $opt_b eq "GC" ) {
    die "fnlo-run-nlojet.pl: Error! Batch system other than GC (grid-control) not updated: $opt_b, aborted.\n";
}
unless ( $opt_e =~ m/\d+/ && $opt_e !~ m/\D+/ ) {
    die "fnlo-run-nlojet.pl: Error! Illegal maximal event number: $opt_e, aborted.\n";
}
unless ( $opt_g eq "guc" || $opt_g eq "srm" ) {
    die "fnlo-run-nlojet.pl: Error! No such grid storage protocol: $opt_g, aborted.\n";
}
unless ( $opt_j =~ m/\d{4}/ && $opt_j !~ m/\D+/ ) {
    die "fnlo-run-nlojet.pl: Error! Illegal job number (nnnn): $opt_j, aborted.\n";
}
unless ( $opt_n eq "2jet" || $opt_n eq "3jet" ) {
    die "fnlo-run-nlojet.pl: Error! Illegal njet mode: $opt_n, aborted.\n";
}
unless ( $opt_o eq "LO" || $opt_o eq "NLO" ) {
    die "fnlo-run-nlojet.pl: Error! Illegal option -o $opt_o, aborted.\n";
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-run-nlojet.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}

my $batch = $opt_b;
my $nmax  = $opt_e;
my $prot  = $opt_g;
my $jobnr = $opt_j;
my $njet  = $opt_n;
my $order = $opt_o;
my $ref   = "";
if ( $opt_r ) { $ref = "ref";}
my $verb  = "";
my $vers  = $opt_v;
my $wrm   = "";
if ( $opt_w ) { $wrm = "wrm";}
if ( $opt_r && $opt_w ) {
    die "fnlo-run-nlojet.pl: Error! Reference and warm-up mode are mutually exclusive!\n";
}
my $seed  = "";
if ( $opt_x ) { $seed = $opt_x ;}
print "fnlo-run-nlojet.pl: Running on batch system $batch\n";
print "fnlo-run-nlojet.pl: Using grid storage protocol $prot\n";
print "fnlo-run-nlojet.pl: Maximal event number: $nmax\n";
print "fnlo-run-nlojet.pl: Attaching job number $jobnr\n";
print "fnlo-run-nlojet.pl: Running in order $order\n";
if ( $wrm eq "wrm" ) {
    print "fnlo-run-nlojet.pl: Running in warm-up mode\n";
}
if ( $ref ) {
    print "fnlo-run-nlojet.pl: Running in reference mode\n";
}
if ( $opt_d ) {
    $verb = 1;
    print "fnlo-run-nlojet.pl: Debug/verbose mode is active.\n";
}

#
# Parse arguments
#
my $scentype = "def";
my $scenname = "def";
unless ( @ARGV == 1 ) {
    die "fnlo-run-nlojet.pl: Error! Need one scenario description in the form of ScenarioType_ScenarioName!\n";
}
if ( @ARGV > 0 ) {
    my $scenstr = shift;
    my @parts = split("_",$scenstr,2);
    $scentype = $parts[0];
    $scenname = $parts[1];
}
print "fnlo-run-nlojet.pl: Running scenario ${scenname} of type ${scentype}.\n";
my $tdir = "${scenname}${ref}_${jobnr}";
if ( $opt_t ) {
    $tdir  = $opt_t;
}
print "fnlo-run-nlojet.pl: Output target directory $tdir\n";

#
# Initialization
#
# Run mode settings: order, ordername, # events
my %runmode;
$runmode{LO}[0]  = "born";
$runmode{LO}[1]  = "100000000";
if ( $nmax > 0 && $nmax < $runmode{LO}[1] ) {
    $runmode{LO}[1]  = "$nmax";
}
$runmode{NLO}[0] = "nlo";
$runmode{NLO}[1] = "10000000";
if ( $nmax > 0 && $nmax < $runmode{NLO}[1] ) {
    $runmode{NLO}[1]  = "$nmax";
}

# NLOJet++ product name
my $prdext = "tab";
my $prdnam = "${scenname}${ref}${wrm}_${jobnr}-hhc-$runmode{$order}[0]-${njet}.${prdext}";
if ( $wrm eq "wrm" ) {
    $prdext = "txt";
    $prdnam = "${scentype}_${scenname}_warmup.${prdext}";
}
print "fnlo-run-nlojet.pl: NLOJet++ product name $prdnam\n";

# Directories
my $rundir = getcwd();
chomp $rundir;

#
# Copy scenario steering file (${scenname}.str) to generic name for scenario type (${scentype}.str)  
#
print "fnlo-run-nlojet.pl: Copy scenario steering ${scenname}.str to ".
    "generic name for scenario type ${scentype}.str\n";
my $ret = system("cp -p ${scenname}.str ${scentype}.str");
if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy scenario steering ${scenname}.str to ".
		 "generic scenario steering name ${scentype}.str: $ret, aborted!\n";}

#
# Print system info (some parts in debug mode only; df commands can get stuck ...)
#
print "\n######################################################\n";
print "fnlo-run-nlojet.pl: System information for debugging purposes:\n";
print "######################################################\n";
my $host = `hostname`;
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"hostname\" command failed.\n\n";
} else {
    chomp $host;
    print "fnlo-run-nlojet.pl: The system's hostname is (hostname):\n$host\n\n";
}
my $osvers = `uname -a`;
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"uname -a\" command failed.\n\n";
} else {
    print "fnlo-run-nlojet.pl: Your operating system is (uname -a):\n$osvers\n";
}
my $procvers = `cat /proc/version`;
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"cat /proc/version\" command failed.\n\n";
} else {
    print "fnlo-run-nlojet.pl: Your linux version is (/proc/version):\n$procvers\n";
}
my $freemem = `free`;
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"free\" command failed.\n\n";
} else {
    print "fnlo-run-nlojet.pl: The available memory is:\n$freemem\n";
}
my $cpumod = `cat /proc/cpuinfo | grep \"model name\"`;
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"cat /proc/cpuinfo\" command failed.\n\n";
} else {
    my $cpufrq = `cat /proc/cpuinfo | grep \"cpu MHz\"`;
    print "fnlo-run-nlojet.pl: The processor type is:\n${cpumod}at\n${cpufrq}\n";
}
if ( $verb ) {
    my $freedisk = `df -h`;
    if ( $? ) {
        print "fnlo-run-nlojet.pl: Info: \"df -h\" command failed.\n\n";
    } else {
        print "fnlo-run-nlojet.pl: The available disk space is:\n$freedisk\n";
    }
    my $freenode = `df -hi`;
    if ( $? ) {
        print "fnlo-run-nlojet.pl: Info: \"df -hi\" command failed.\n\n";
    } else {
        print "fnlo-run-nlojet.pl: The available inode space is:\n$freenode\n";
    }
}
my $cwd = getcwd();
if ( $? ) {
    print "fnlo-run-nlojet.pl: Info: \"getcwd()\" command failed.\n\n";
} else {
    print "fnlo-run-nlojet.pl: The current working directory is:\n$cwd\n\n";
    print "fnlo-run-nlojet.pl: The current working directory's content is:\n";
    my $ret = system("ls -laR");
    if ( $ret ) {print "fnlo-run-nlojet.pl: Couldn't list current directory: $ret, skipped!\n";}
}
print "######################################################\n\n";

#
# Set system paths environment
#
$ENV{FASTJET} = "$rundir";
$ENV{LHAPDF}  = "$rundir";
$ENV{NLOJET}  = "$rundir";
$ENV{FASTNLO} = "$rundir";
if ( $vers eq "2.3" ) {
    $ENV{FASTNLOLIBADDDIR} = "$rundir/lib/fastnlo_interface_nlojet";
} else {
    die "fnlo-run-nlojet.pl: ERROR! Unsupported fastNLO version $vers requested, aborted!\n";
}
if ( $ENV{PATH} ) {
    $ENV{PATH} = "$ENV{NLOJET}/bin:$ENV{PATH}";
} else {
    $ENV{PATH} = "$ENV{NLOJET}/bin";
}
if ( $ENV{LD_LIBRARY_PATH} ) {
    $ENV{LD_LIBRARY_PATH} ="$ENV{NLOJET}/lib:$ENV{NLOJET}/lib64:$ENV{FASTNLOLIBADDDIR}:".
        "$ENV{LD_LIBRARY_PATH}";
} else {
    $ENV{LD_LIBRARY_PATH} ="$ENV{NLOJET}/lib:$ENV{NLOJET}/lib64:$ENV{FASTNLOLIBADDDIR}";
}

#
# Run fastNLO
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nlojet.pl: Running fastNLO version $vers scenario on batch system $batch: $date\n";

# In case of non-local running need to unpack binary fastNLO distribution
if ( $batch ne "LOCAL" ) {

# Fetching and unpacking of fastNLO binary archive
    my $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-nlojet.pl: Get/unpack fastNLO binary package: BINGET0_$date\n";

    my $file = "fastNLO-bin-v${vers}.tgz";
    if ( ! -f $file ) {
        if ( $batch eq "GC" ) {
            die "fnlo-run-nlojet.pl: ERROR! Could not find binary tgz of fastNLO $file working dir, aborted!\nCheck your grid storage options for grid-control!\n";
        } elsif ( $batch eq "GRID" ) {
            grid_storage("FETCH","fastNLO_archives",$file,".",$file,$prot);
        } elsif ( $batch eq "PBS" ) {
            die "fnlo-run-nlojet.pl: ERROR! Could not find binary tgz of fastNLO $file in PBS working dir, aborted!\n";
        }
    }
    if ( -f $file ) {
        system ("tar xfz $file");
    } else {
        die "fnlo-run-nlojet.pl: ERROR! Could not find binary tgz of fastNLO $file, aborted!\n";
    }

    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-nlojet.pl: Unpacked fastNLO binary package: BINGET1_$date\n";
}

my $scendir = "$rundir";
chdir "$scendir" or die
    "fnlo-run-nlojet.pl: ERROR! Could not cd to dir $scendir!\n";

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nlojet.pl: Running fastNLO: $date\n";
my $cmd;
if ( $vers eq "2.3" ) {
    $cmd = "$ENV{NLOJET}/bin/nlojet++ --calculate ".
        "--save-after $runmode{$order}[1] ".
        "-c$runmode{$order}[0] ".
        "-d $tdir ".
        "-n ${scenname}${ref}${wrm}_${jobnr} ".
        "-u $ENV{NLOJET}/lib/fastnlo_interface_nlojet/lib${scentype}.la ";
    if ( $nmax ) { $cmd .= "--max-event $nmax "; }
    if ( $seed ) { $cmd .= "-s $seed "; }
} else {
    die "fnlo-run-nlojet.pl: ERROR! Unsupported fastNLO version $vers requested, aborted!\n";
}

# Run NLO calculation
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nlojet.pl: Starting calculation: FASTCAL0_$date\n";
print "\nfnlo-run-nlojet.pl: Running command (time $cmd) 2>&1 in foreground\n";
$ret = system("(time $cmd) 2>&1");
if ( $ret ) {die "fnlo-run-nlojet.pl: ERROR! Error $ret in fastNLO run step, aborted!\n";}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-run-nlojet.pl: Calculation finished: FASTCAL1_$date\n";

# Copy/rename results for grid storage
print "\nfnlo-run-nlojet.pl: The final working directory's content is:\n";
$ret = system("ls -laR");
if ( $ret ) {print "fnlo-run-nlojet.pl: Couldn't list current directory: $ret, skipped!\n";}
my $spath = "${scendir}/${tdir}";
print "\nfnlo-run-nlojet.pl: Source path: $spath\n";
my $sfile = "${prdnam}";
print "fnlo-run-nlojet.pl: Source file: $sfile\n";
my $tpath = "fastNLO_tables_v${vers}/${scentype}_${scenname}";
$tpath =~ s/\/\.\//\//g;
print "fnlo-run-nlojet.pl: Target path: $tpath\n";
my $tfile = "${sfile}";
# Eliminate default job number and attach the GC one
$tfile =~ s/_0001//;
$tfile =~ s/\.${prdext}/_${gjobnr}\.${prdext}/;
print "fnlo-run-nlojet.pl: Target file: $tfile\n";

if ( $batch ne "LOCAL" ) {
# Copy/rename for grid storage via GC
    if ( $batch eq "GC" ) {
	print "fnlo-run-nlojet.pl: Info: Batch mode $batch: Grid storage done by grid-control.\n";
	my $spathdir = `dirname $spath`;
	if ( ! -f "$rundir/$sfile" ) {
	    print "fnlo-run-nlojet.pl: Copy product to current directory for storage.\n";
	    my $ret = system("cp -p $spath/$sfile $rundir/$sfile");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy product into ".
			     "current directory $rundir: $ret, aborted!\n";}
	}
	chdir $rundir or die "fnlo-run-nlojet.pl: ERROR! Couldn't cd to $rundir, aborted!\n";

	print "fnlo-run-nlojet.pl: Rename product to expected file name for storage.\n";
	my $ret = system("mv -f $sfile $tfile");
	if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't rename product ${sfile} into ".
			 "${tfile}: $ret, aborted!\n";}

# Also rename related log and err files in cwd for grid storage, if file size larger than zero
# Remove product extension to rename log/err files
	$tfile =~ s/\.${prdext}//;
	print "fnlo-run-nlojet.pl: Copy log files to correct file names for storage.\n\n";
	if ( -f "job.stderr" && ! -z "job.stderr" ) {
	    my $ret = system("cp -p job.stderr ${tfile}.err");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: ERROR! Couldn't copy job.stderr ".
			     "to ${tfile}.err: $ret, aborted!\n";}
	}
	if ( -f "job.stdout" && ! -z "job.stdout" ) {
	    $ret = system("cp -p job.stdout ${tfile}.log");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: ERROR! Couldn't copy job.stdout ".
			     "to ${tfile}.log: $ret, aborted!\n";}
	}
	$ret = system("pwd");
	if ( $ret ) {print "fnlo-run-nlojet.pl: WARNING! Couldn't print cwd!\n";}
	$ret = system("ls -la");
	if ( $ret ) {print "fnlo-run-nlojet.pl: WARNING! Couldn't list cwd!\n";}
# Do grid storage by ourselves (unused since quite some time, beware)
    } elsif ( $batch eq "GRID" ) {
	grid_storage("TABSAV","$spath","$sfile","$tpath","$tfile",$prot);
	if ( -f "job.stderr" ) {
	    my $ret = system("cp -p job.stderr ${tfile}.err");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy job.stderr ".
			     "fastrun_${gjobnr}.err: $ret, aborted!\n";}
	}
	if ( -f "job.stdout" ) {
	    $ret = system("cp -p job.stdout ${tfile}.log");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy job.stdout ".
			     "fastrun_${gjobnr}.log: $ret, aborted!\n";}
	}
	grid_storage("LOGSAV","$rundir","fastrun_${gjobnr}.err","$tpath","${tfile}.err",$prot);
	grid_storage("LOGSAV","$rundir","fastrun_${gjobnr}.log","$tpath","${tfile}.log",$prot);
# PBS (unused since quite some time, beware)
    } else {
	my $ret = system("cp -p $spath/$sfile $rundir/$tfile");
	if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy table into ".
			 "current directory $rundir: $ret, aborted!\n";}
	if ( -f "job.stderr" ) {
	    my $ret = system("cp -p job.stderr ${tfile}.err");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy job.stderr ".
			     "fastrun_${gjobnr}.err: $ret, aborted!\n";}
	}
	if ( -f "job.stdout" ) {
	    $ret = system("cp -p job.stdout ${tfile}.log");
	    if ( $ret ) {die "fnlo-run-nlojet.pl: Couldn't copy job.stdout ".
			     "fastrun_${gjobnr}.log: $ret, aborted!\n";}
	}
    }
    my $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    print "\nfnlo-run-nlojet.pl: Results stored: TABSAV1_$date\n";
}

$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\n###############################\n";
print "# fnlo-run-nlojet.pl: fastNLO finished: FASTRUN1_$date\n";
print "###############################\n\n";
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
    print "fnlo-run-nlojet.pl: spath $spath\n";
    print "fnlo-run-nlojet.pl: sfile $sfile\n";
    print "fnlo-run-nlojet.pl: tpath $tpath\n";
    print "fnlo-run-nlojet.pl: tfile $tfile\n";
    print "fnlo-run-nlojet.pl: protocol $protocol\n";

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
            print "fnlo-run-nlojet.pl: ERROR! ".
                "globus-url-copy not found: $ret!\n";}
        print "fnlo-run-nlojet.pl: Found command $guccmd version $gucver\n";
        $gucopt = "-cd";
        if ( $gucver <= 2.9 ) {
            print "fnlo-run-nlojet.pl: WARNING! Version of globus-url-copy ".
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
        print "fnlo-run-nlojet.pl: Found srmcp command version $srmver\n";
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

    print "fnlo-run-nlojet.pl: Source: $source\n";
    print "fnlo-run-nlojet.pl: Target: $target\n";

    my $gjobnr = "";
    if ( defined $ENV{MY_JOBID} ) {
        $gjobnr = $ENV{MY_JOBID};
        $gjobnr = substr("0000$gjobnr",-4);
    }

    $date = `date +%d%m%Y_%H%M%S`;
    chomp $date;
    if ( $signam eq "TABSAV" || $signam eq "LOGSAV" ) {
        print "\nfnlo-run-nlojet.pl: Received signal $signam: ${signam}0_$date\n";
        print "fnlo-run-nlojet.pl: Saving source $source to\n";
        print "fnlo-run-nlojet.pl: target $target for\n";
        print "fnlo-run-nlojet.pl: fastNLO job no. $jobnr and\n";
        print "fnlo-run-nlojet.pl: grid-control job no. $gjobnr\n";
        print "fnlo-run-nlojet.pl: Trying to find source $source ...\n";
        system("pwd");
        system("ls -la $spath");
        unless ( -f "$source" ) {
            die "fnlo-run-nlojet.pl: ERROR! Could not find source $source\n";
        }
        if ( -z "$source" ) {
            print "fnlo-run-nlojet.pl: WARNING! Skipping storage of empty file $source\n";
        } else {
            if ($protocol eq "guc") {
                my $gcmd = "globus-url-copy $gucopt ".
                    "file://${source} ".
                    "gsiftp://${sename}/${target}";
                print "Command $gcmd\n";
                my $ret = system("$gcmd");
                if ( $ret ) {die "fnlo-run-nlojet.pl: ERROR! Grid storage of ".
                                 "source $source to target $target failed: $ret!\n";}
            } else {
                $ENV{SRM_PATH} = "";
                my $gcmd = "srmcp ".
                    "file:////${source} ".
                    "srm://${sename}/${target}";
                print "Command $gcmd\n";
                my $ret = system("$gcmd");
                if ( $ret ) {die "fnlo-run-nlojet.pl: ERROR! Grid storage of ".
                                 "source $source to target $target failed: $ret!\n";}
            }
        }
    } elsif ( $signam eq "FETCH" ) {
        print "fnlo-run-nlojet.pl: Fetching binary fastNLO archive $sfile from\n";
        print "fnlo-run-nlojet.pl: SE $sename in path $sepath\n";
        print "fnlo-run-nlojet.pl: Trying to fetch binary fastNLO archive $source on\n";
        print "fnlo-run-nlojet.pl: SE $sename in path $sepath\n";
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
        die "fnlo-run-nlojet.pl: ERROR! Caught unsupported signal $signam, aborted!\n";
    }

    return 0;
}
