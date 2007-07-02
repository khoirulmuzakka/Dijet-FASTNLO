#!/usr/bin/env perl 
#
# fastNLO table addition script
# Version:
# 
# created by K. Rabbertz: 05.03.2006
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
print "\n##################################################\n";
print "# fastadd.pl: Starting table addition for fastNLO: $date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_f, $opt_h ) = ( ".", "187", "" );
getopts('d:f:h') or die "fastadd.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\nfastadd.pl\n";
    print "Usage: fastadd.pl [switches/options] scenario\n";
    print "  -d dir          Installation directory (def.=.)\n";
    print "  -f rev          fastNLO revision to use (def.=187)\n";
    print "  -h              Print this text\n\n";
    exit;
}

unless ( -d $opt_d ) {
    die "fastadd.pl: Error! No such directory: $opt_d, aborted.\n";
}
unless ( $opt_f =~ m/\d+/ && $opt_f !~ m/\D+/ ) {
    die "fastadd.pl: Error! Illegal fastNLO revision number: $opt_f, aborted.\n";
}

my $idir  = $opt_d;
my $frev  = $opt_f;
print "fastadd.pl: Installation directory $idir\n";
print "fastadd.pl: Using fastNLO revision $frev\n";

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fastadd.pl: Error! Need one scenario name!\n";
}
my $scen   = shift;

#
# Initialization
#
my $lodir   = "${scen}_LO_tables"; 
my $nlodir  = "${scen}_NLO_tables"; 
my $loglob  = "${scen}*born*.raw*";
my $nloglob = "${scen}*nlo*.raw*";

# Directories
my $sdir = getcwd();

# Define install hash
my %install;
# First entry (index 0!): Subdirecory name into which the archive is unpacked!
$install{cernlib}[0]    = "cernlib-2003";
$install{lhapdf}[0]     = "lhapdf-5.2.3";
$install{nlojet}[0]     = "nlojet++-2.0.1";
$install{nlojetfix}[0]  = "nlojet++-2.0.1";
$install{fastNLO}[0]    = "fastNLO-rev${frev}";

#
# Find LO tables
#
chdir "$lodir" or die "fastadd.pl: ERROR! Couldn't cd to $lodir\n";
my @lotabs = glob $loglob;
#print "fastadd.pl: DEBUG! lotabs @lotabs\n";

#
# Find NLO tables
#
chdir "../${nlodir}" or die "fastadd.pl: ERROR! Couldn't cd to ../${nlodir}\n";
my @nlotabs = glob $nloglob;
#print "fastadd.pl: DEBUG! nlotabs @nlotabs\n";
chdir "..";

#
# Make fastNLO add
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: Setting environment for fastNLO: $date\n";
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
chdir "$install{fastNLO}[0]/author1c/hadron";
    
my $deb = `pwd`;
chomp $deb;
print "fastadd.pl: DEBUG! Current directory to run NLOJET++: $deb\n";
print "fastadd.pl: DEBUG! ls -la:\n\n";
system("ls -la");
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: Making nlofast-add: $date\n";
system("make nlofast-add");

#
# nlofast-add
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: nlofast-add: $date\n";
chdir "$sdir";
my $cmd = "$ENV{NLOJET}/bin/nlofast-add ";
foreach my $lotab ( @lotabs ) {
    $cmd .= "$lodir/$lotab ";
}
foreach my $nlotab ( @nlotabs ) {
    $cmd .= "$nlodir/$nlotab ";
}
$cmd .= "$scen.tab";
print "fastrun.pl: Running command $cmd\n";
system("$cmd >& $scen.log");
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfastadd.pl: fastNLO finished: $date\n";
exit 0;
