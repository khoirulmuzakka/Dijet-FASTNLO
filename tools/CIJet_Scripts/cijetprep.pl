#!/usr/bin/env perl
#
# CIJET grid submission archiver
# Version:
#
# created by K. Rabbertz: 20.02.2014
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
print "# cijetprep.pl: Starting archive creation for CIJET: CIJETPREP_$date\n";
print "#####################################################\n\n";

#
# Parse options
#
our ( $opt_h, $opt_l, $opt_n ) = ( "", "cteq6ll.LHpdf", "cteq66.LHgrid" );
getopts('hl:n:') or die "cijetprep.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "\ncijetprep.pl\n";
    print "Usage: cijetprep.pl [switches/options]\n";
    print "       WARNING! cijetprep.pl works only if all required software is installed \n";
    print "       under a common prefix path given either via \$CIJET or the current \n";
    print "       working directory!\n\n";
    print "  -h              Print this text\n";
    print "  -l pdf          Filename of LO LHAPDF set to add, def. = cteq6ll.LHpdf\n";
    print "  -n pdf          Filename of NLO LHAPDF set to add, def. = cteq66.LHgrid\n\n";
    exit;
}
my $pdfl    = $opt_l;
my $pdfn    = $opt_n;
my $arcname = "cijet-bin";

#
# Starting archive creation
#
if ( ! $ENV{CIJET} ) {
    print "cijetprep.pl: WARNING! Environment variable CIJET not set!\n";
    print "             Assume current directory to contain installation!\n";
    $ENV{CIJET} = getcwd();
}

print "cijetprep.pl: Preparing CIJET archive for submission in directory $ENV{CIJET}/..\n";
chdir "$ENV{CIJET}" or die "fastprep.pl: ERROR! Could not cd to $ENV{CIJET}!\n";
my @cijetexe = `find . -name dijets4ci_sin`;
my $cijetexe = `basename $cijetexe[0]`;
chomp $cijetexe;
my $cijetdir = `dirname $cijetexe[0]`;
chomp $cijetdir;
print "cijetprep.pl: cijet executable $cijetexe found in $cijetdir ...\n";
my @pdflfile = `find -L share/lhapdf -name $pdfl`;
my $pdflfile = `basename $pdflfile[0]`;
chomp $pdflfile;
my $pdfldir = `dirname $pdflfile[0]`;
chomp $pdfldir;
print "cijetprep.pl: LO PDF $pdflfile found in $pdfldir ...\n";
my @pdfnfile = `find -L share/lhapdf -name $pdfn`;
my $pdfnfile = `basename $pdfnfile[0]`;
chomp $pdfnfile;
my $pdfndir = `dirname $pdfnfile[0]`;
chomp $pdfndir;
print "cijetprep.pl: NLO PDF $pdfnfile found in $pdfndir ...\n";

my $cmd = "tar cf ${arcname}.tar -C $cijetdir $cijetexe";
my $ret = system("$cmd");
if ( $ret ) {die "cijetprep.pl: ERROR! Could not create archive ${arcname}.tar: $ret\n";}
$cmd = "tar rf ${arcname}.tar -h -C $pdfldir $pdflfile";
$ret = system("$cmd");
if ( $ret ) {die "cijetprep.pl: ERROR! Could not append $pdflfile to archive ${arcname}.tar: $ret\n";}
$cmd = "tar rf ${arcname}.tar -h -C $pdfndir $pdfnfile";
$ret = system("$cmd");
if ( $ret ) {die "cijetprep.pl: ERROR! Could not append $pdfnfile to archive ${arcname}.tar: $ret\n";}
$cmd = "gzip ${arcname}.tar";
$ret = system("$cmd");
if ( $ret ) {die "cijetprep.pl: ERROR! Could not gzip archive ${arcname}.tar: $ret\n";}
$cmd = "mv ${arcname}.tar.gz ${arcname}.tgz";
$ret = system("$cmd");
if ( $ret ) {die "cijetprep.pl: ERROR! Could not rename archive ${arcname}.tar: $ret\n";}

exit 0;
