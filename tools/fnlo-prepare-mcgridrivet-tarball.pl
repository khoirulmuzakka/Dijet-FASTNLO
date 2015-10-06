#!/usr/bin/env perl
#
# fastNLO MCgrid Rivet analysis archiver
# Version:
#
# created by K. Rabbertz: 10.12.2007
# adapted by K. Rabbertz from fnlo-prepare-tarball.pl: 06.10.2015
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
my $date = `date +%d%m%Y_%H%M%S`;
if ( $? ) {die "fnlo-prepare-mcgridrivet-tarball.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n#####################################################\n";
print "# fnlo-prepare-mcgridrivet-tarball.pl: Starting archive creation for MCgrid Rivet analysis: FASTPREP_$date\n";
print "#####################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_v ) = ( "", "2.3" );
getopts('hv:') or die "fnlo-prepare-mcgridrivet-tarball.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-prepare-mcgridrivet-tarball.pl\n";
    print "Usage: fnlo-prepare-mcgridrivet-tarball.pl [switches/options] scenario\n";
    print "       Here, scenario is equivalent to the Rivet analysis name resp. the directory.\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-prepare-mcgridrivet-tarball.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-prepare-mcgridrivet-tarball.pl: Error! Need one scenario name!\n";
}
my $scen    = shift;
my $arcname = "$scen.tgz";

#
# Starting archive creation
#
if ( ! -d $scen ) {
    print "fnlo-prepare-mcgridrivet-tarball.pl: ERROR! Directory $scen not found!\n";
    print "    Please cd to the directory with the $scen MCgrid Rivet analysis.\n";
    exit(1);
}

print "fnlo-prepare-mcgridrivet-tarball.pl: Preparing MCgrid/Rivet archive $arcname for running on worker nodes.\n";
my @basfiles = ( "$scen/Makefile", "$scen/Results.db", "$scen/Run_prod.dat", "$scen/setup.csh" );
my @strfiles = `find $scen -follow -name \*.str`;
chomp @strfiles;
my @subdirs  = ( "$scen/Process", "$scen/analysis" );
my @prcfiles = `find $scen -follow -name \*.evtcount`;
chomp @prcfiles;
my @wrmfiles = `find $scen -follow -name \*.txt`;
chomp @wrmfiles;

my $cmd = "tar cfzh $arcname @basfiles @strfiles @subdirs @prcfiles @wrmfiles";
print "fnlo-prepare-mcgridrivet-tarball.pl: Executing: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-prepare-mcgridrivet-tarball.pl: ERROR! Could not create archive $arcname: $ret\n";}

exit 0;
