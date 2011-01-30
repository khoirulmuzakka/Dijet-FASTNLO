#!/usr/bin/env perl 
#
# fastNLO grid submission archiver
# Version:
# 
# created by K. Rabbertz: 10.12.2007
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
print "\n#####################################################\n";
print "# fastprep.pl: Starting archive creation for fastNLO: FASTPREP_$date\n";
print "#####################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_p, $opt_v ) = ( "", "", "1" );
getopts('hp:v:') or die "fastprep.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastprep.pl\n";
#    print "Usage: fastprep.pl [switches/options] scenario\n";
    print "Usage: fastprep.pl [switches/options]\n";
    print "       WARNING! fastprep.pl works only if all required software is installed \n";
    print "       under a common prefix path given either via \$FASTNLO or the current \n";
    print "       working directory!\n\n";
    print "  -h              Print this text\n";
    print "  -p pdf          Add CTEQ parton densities or LHAPDF\n";
    print "  -v #            Choose between fastNLO version 1 or 2 (def.=1)\n\n";
    exit;
}
unless ( $opt_p eq "" || $opt_p eq "CTEQ" || $opt_p eq "LHAPDF" ) {
    die "fastrun.pl: Error! Illegal option -p $opt_p, aborted.\n";
}
my $pdf   = $opt_p;
my $vers  = $opt_v;
my $arcname = "fastNLO-bin.tgz";
if ( $vers == 1 && $pdf ) {
    $arcname = "fastNLO-bin-${pdf}.tgz";
} elsif ( $vers == 2 ) {
    $arcname = "fastNLO-bin-v${vers}.tgz";
}

#
# Starting archive creation
#
if ( ! $ENV{FASTNLO} ) {
    print "fastprep.pl: WARNING! Environment variable FASTNLO not set!\n";
    print "             Assume current directory to contain installation!\n";
    $ENV{FASTNLO} = getcwd();
}

print "fastprep.pl: Preparing fastNLO archive for submission in directory $ENV{FASTNLO}/..\n";
print "fastprep.pl: Only tested for hh collisions so far ...\n";

if ( $vers == 1 ) {
    chdir "$ENV{FASTNLO}/.." or die "fastprep.pl: ERROR! Could not cd to $ENV{FASTNLO}/..!\n";

    my @gcclibs = `find lib -follow -name \*.so\*`;
    chomp @gcclibs;
    my @njlibs  = `find nlojet*/lib -follow -name \*.so\*`;
    chomp @njlibs;
    my @fjlibs  = `find fastjet*/lib -follow -name \*.so\*`;
    chomp @fjlibs;
    my @fjplugs = `find fastjet*/plugins -follow -name \*.so\*`;
    chomp @fjplugs;
    my @fnlibs1 = `find fastNLO*/trunk -follow -name \*.so\*`;
    chomp @fnlibs1;
    my @fnlibs2 = `find fastNLO*/trunk -follow -name \*.la`;
    chomp @fnlibs2;
    my $nloexe  = `find nlojet*/bin -follow -name nlojet++`;
    chomp $nloexe;

    my $cmd = "tar cfz $arcname --exclude .svn ".
	"@gcclibs @njlibs @fjlibs @fjplugs @fnlibs1 @fnlibs2 $nloexe";
    if ( -d "lib64" ) {
	my @libs64 = `find lib64 -follow -name \*.so\*`;
	chomp @libs64;
	$cmd .= " @libs64";
    }
    if ( $pdf eq "CTEQ" ) {
	my @ctqlnk = `find fastNLO*/trunk/v1.4/author1c/hadron -name ctq61.00.tbl`;
	chomp @ctqlnk;
	my @ctqtbl = `find fastNLO*/trunk/v1.4/author1c/common -name ctq61.00.tbl`;
	chomp @ctqtbl;
	my @comlnk = `find fastNLO*/trunk/v1.4/author1c/hadron -name common`;
	chomp @comlnk;
	$cmd .= " @ctqlnk @ctqtbl @comlnk";
#	$cmd .= " fastNLO/trunk/v1.4/author1c/common/".
#	    " fastNLO/trunk/v1.4/author1c/hadron/common".
#	    " fastNLO/trunk/v1.4/author1c/hadron/ctq61.00.tbl";
    } elsif ( $pdf eq "LHAPDF" ) {
	$cmd .= " lhapdf/lib";
    }
    my $ret = system("$cmd");
    if ( $ret ) {die "fastprep.pl: ERROR! Could not create archive $arcname: $ret\n";}
} else {
    chdir "$ENV{FASTNLO}" or die "fastprep.pl: ERROR! Could not cd to $ENV{FASTNLO}!\n";
# For version 2:
# Unfortunately need to copy some system libs in case of insufficiently equipped target systems
    print "fastprep.pl: Copying some systems lib(s) into directory $ENV{FASTNLO}/lib ...\n";
    print "             This avoids crashes on poor target systems without these.\n";
    my @libc = `ldd bin/nlojet++ | grep libc | cut -d " " -f3`;
    chomp @libc;
    my @libm = `ldd bin/nlojet++ | grep libm | cut -d " " -f3`;
    chomp @libm;
    my @libstdc = `ldd bin/nlojet++ | grep libstdc | cut -d " " -f3`;
    chomp @libstdc;
    my @libg2c = `ldd bin/nlojet++ | grep libg2c | cut -d" " -f3`;
    chomp @libg2c;
    my @libgcc = `ldd bin/nlojet++ | grep libgcc | cut -d " " -f3`;
    chomp @libgcc;
    my @libdl = `ldd bin/nlojet++ | grep libdl | cut -d " " -f3`;
    chomp @libdl;
    my @libltdl = `ldd bin/nlojet++ | grep libltdl | cut -d " " -f3`;
    chomp @libltdl;
    my @sodeps = (@libc,@libm,@libstdc,@libg2c,@libgcc,@libltdl,@libdl);

    print "fastprep.pl: Adding system lib(s) @sodeps if not yet done.\n";
# Find all libs and links corresponding to the detected sodeps ...
    foreach my $sodep ( @sodeps ) {
	my $dir = `dirname $sodep`;
	chomp $dir;
#	print "dir $dir\n";
	my $lib = `basename $sodep`;
	chomp $lib;
#	print "lib $lib\n";
	my @parts = split(/\./,$lib);
	print "parts @parts\n";
	my @addlibs = `find $dir -name $parts[0].so\*`;
	print "addlibs @addlibs\n";
	chomp @addlibs;
	foreach my $copy ( @addlibs ) {
	    my $lib = `basename $copy`;
	    chomp $lib;
	    if (! -e "lib/$lib" ) {
		my $cmd = "cp -P $copy lib";
		my $ret = system("$cmd");
		if ( $ret ) {print "fastprep.pl: Warning! ".
				 "Copying system lib $copy failed.\n"}
	    }
	}
    }

    my @solibs  = `find lib -follow -name \*.so\*`;
    chomp @solibs;
    my @lalibs  = `find lib -follow -name \*.la\*`;
    chomp @lalibs;

    my $cmd = "tar cfz $arcname ".
	"@solibs @lalibs bin/nlojet++";
    my $ret = system("$cmd");
    if ( $ret ) {die "fastprep.pl: ERROR! Could not create archive $arcname: $ret\n";}
}

exit 0;
