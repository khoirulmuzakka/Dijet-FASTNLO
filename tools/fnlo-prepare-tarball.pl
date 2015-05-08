#!/usr/bin/env perl
#
# fastNLO grid submission archiver
# Version:
#
# created by K. Rabbertz: 10.12.2007
# adapted by K. Rabbertz from fastprep.pl: 27.04.2015
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
if ( $? ) {die "fnlo-prepare-tarball.pl: Error! date command failed, timing not possible,".
               " aborted.\n";}
chomp $date;
print "\n#####################################################\n";
print "# fnlo-prepare-tarball.pl: Starting archive creation for fastNLO: FASTPREP_$date\n";
print "#####################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_v ) = ( "", "2.3" );
getopts('hv:') or die "fnlo-prepare-tarball.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfnlo-prepare-tarball.pl\n";
    print "Usage: fnlo-prepare-tarball.pl [switches/options]\n";
    print "       WARNING! fnlo-prepare-tarball.pl works only if all required software is installed \n";
    print "       under a common prefix path given either via \$FASTNLO or the current \n";
    print "       working directory!\n\n";
    print "  -h              Print this text\n";
    print "  -v #            Choose between fastNLO versions 2.3 and ? (def.=2.3)\n\n";
    exit;
}
unless ( $opt_v eq "2.3" ) {
    die "fnlo-clean-tables.pl: Error! Unsupported fastNLO $opt_v, aborted.\n";
}
my $vers  = $opt_v;
my $arcname = "fastNLO-bin-v${vers}.tgz";

#
# Starting archive creation
#
if ( ! $ENV{FASTNLO} ) {
    print "fnlo-prepare-tarball.pl: WARNING! Environment variable FASTNLO not set!\n";
    print "             Assume current directory to contain installation!\n";
    $ENV{FASTNLO} = getcwd();
}

print "fnlo-prepare-tarball.pl: Preparing fastNLO archive for submission in directory $ENV{FASTNLO}/..\n";
chdir "$ENV{FASTNLO}" or die "fnlo-prepare-tarball.pl: ERROR! Could not cd to $ENV{FASTNLO}!\n";
my @solibs  = `find lib -follow -name \*.so\*`;
chomp @solibs;
my @lalibs  = `find lib -follow -name \*.la\*`;
chomp @lalibs;

# Exclude root libs to save some space
my @libs;
foreach my $lib ( @lalibs, @solibs ) {
    if (! ($lib =~ m/root/ )) {
        push @libs, $lib;
    }
}

my $cmd = "tar cfz $arcname @libs bin/nlojet++";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-prepare-tarball.pl: ERROR! Could not create archive $arcname: $ret\n";}

exit 0;
