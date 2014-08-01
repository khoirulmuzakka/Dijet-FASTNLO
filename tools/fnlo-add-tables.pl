#!/usr/bin/env perl
#
# fastNLO table addition script:
#    This perl script creates either a total sum table (scen.tab) or
#    sum tables for statistical uncertainty evaluation:
#    - For each single LO table with *one* 'dummy' NLO table
#    - For each single NLO table with *all* LO tables or
#
# Version:
#
# created by K. Rabbertz: 05.03.2006
# last modified:
# adapted by K. Rabbertz from fastadd.pl: 01.08.2014
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
chomp $date;
print "\n##################################################\n";
print "# fnlo-add-tables.pl: Starting table addition for fastNLO: TABADD_$date\n";
print "##################################################\n\n";

#
# Parse options
#
our ( $opt_d, $opt_h, $opt_l, $opt_n, $opt_s, $opt_v, $opt_w ) = ( "", "", "", "", "", "2.3", "" );
getopts('dhl:n:sv:') or die "fnlo-add-tables.pl: Malformed option syntax!\n";
if ( $opt_h ) {
    print "fnlo-add-tables.pl\n";
    print "Usage: fnlo-add-tables.pl [switches/options] scenario\n";
    print "  -d              Verbose output\n";
    print "  -h              Print this text\n";
    print "  -l dir          Directory for LO tables, (def.=scenario)\n";
    print "  -n dir          Directory for NLO/NNLO tables, (def.=scenario)\n";
    print "  -s              Produce tables for statistical evaluation,\n".
        "                  i.e. combinations of each LO with 1 NLO table and\n".
        "                  all LO with each NLO table\n";
    print "  -v #            Choose between fastNLO version 2.3 or ? (def.=2.3)\n\n";
    exit;
}

#
# Parse arguments
#
unless ( @ARGV == 1 ) {
    die "fnlo-add-tables.pl: Error! Need one scenario name!\n";
}
my $scen   = shift;
my $debug  = $opt_d;
my $vers   = $opt_v;

#
# Initialization
#
my $lodir    = "${scen}";
if ( $opt_l ) {$lodir = $opt_l;}
my $nlodir   = "${scen}";
if ( $opt_n ) {$nlodir = $opt_n;}
my $nnlodir  = $nlodir;
my $tabext = "tab";
my $loglob   = "${scen}*born*.${tabext}*";
my $nloglob  = "${scen}*nlo*.${tabext}*";
my $nnloglob = "${scen}*thrcor*.${tabext}*";

# Directory
my $sdir = getcwd();

#
# Check on fnlo-tk-merge
#
my $merger = "fnlo-tk-merge";
my $cmd;
$cmd = `which ${merger}`;
chomp $cmd;
unless ( $cmd ) {
    if ( $ENV{NLOJET} ) {
        if ( -f "$ENV{NLOJET}/bin/${merger}" ) {
            $cmd = "$ENV{NLOJET}/bin/${merger}";
        } else {
            die "fnlo-add-tables.pl: ERROR! ${merger} command not found, ".
                "neither via \`which\` nor in $ENV{NLOJET}/bin, aborted!\n";
        }
    } else {
        die "fnlo-add-tables.pl: ERROR! ${merger} command not found ".
            "via \`which\` and \$NLOJET is not set, aborted!\n";
    }
}

#
# Find LO tables
#
my @lotabs;
if ( -d "$lodir" ) {
    chdir $lodir;
    @lotabs = glob $loglob;
    chdir $sdir;
}
unless ( @lotabs ) {
    print "fnlo-add-tables.pl: No LO table found in $lodir, now looking in $sdir ...\n";
    $lodir  = ".";
    @lotabs = glob $loglob;
}
unless ( @lotabs ) {
    die "fnlo-add-tables.pl: ERROR! No LO table found, aborted!\n";
}
my $ntab = scalar @lotabs;
print "fnlo-add-tables.pl: $ntab LO tables found.\n";
if ( $debug ) {print "fnlo-add-tables.pl: DEBUG! lotabs @lotabs\n";}

#
# Find NLO tables
#
my @nlotabs;
if ( -d "$nlodir" ) {
    chdir $nlodir;
    @nlotabs = glob $nloglob;
    chdir $sdir;
}
unless ( @nlotabs ) {
    print "fnlo-add-tables.pl: No NLO table found in $nlodir, now looking in $sdir ...\n";
    $nlodir  = ".";
    @nlotabs = glob $nloglob;
}
unless ( @nlotabs ) {
    die "fnlo-add-tables.pl: ERROR! No NLO table found, aborted!\n";
}
$ntab = scalar @nlotabs;
print "fnlo-add-tables.pl: $ntab NLO tables found.\n";
if ( $debug ) {print "fnlo-add-tables.pl: DEBUG! nlotabs @nlotabs\n";}

#
# Find NNLO tables
#
my @nnlotabs;
if ( -d "$nnlodir" ) {
    chdir $nnlodir;
    @nnlotabs = glob $nnloglob;
    chdir $sdir;
}
unless ( @nnlotabs ) {
    print "fnlo-add-tables.pl: No NNLO table found in $nnlodir, now looking in $sdir ...\n";
    $nnlodir  = ".";
    @nnlotabs = glob $nnloglob;
}
unless ( @nnlotabs ) {
    print "fnlo-add-tables.pl: WARNING! No NNLO table found!\n";
}
$ntab = scalar @nnlotabs;
print "fnlo-add-tables.pl: $ntab NNLO tables found.\n";
if ( $debug ) {print "fnlo-add-tables.pl: DEBUG! nnlotabs @nnlotabs\n";}

#
# fnlo-tk-merge
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-tables.pl: ${merger}: $date\n";
# Statistics mode: Each LO with first NLO table & all LO with each NLO table
if ( $opt_s ) {
    print "\nfnlo-add-tables.pl: Statistics mode\n";
    my $scmd2 = $cmd;
    foreach my $lotab ( @lotabs ) {
        my $scmd1 = $cmd;
        $scmd1 .= " $lodir/$lotab $nlodir/$nlotabs[0]";
        $scmd2 .= " $lodir/$lotab";
        $scmd1 .= " $lotab";
        $scmd1 =~ s/\.${tabext}$/\.tab/;
        print "fnlo-add-tables.pl: Creating sum table for $lotab ...\n";
        if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
        system("$scmd1 >> ${scen}_addst.log");
    }
    foreach my $nlotab ( @nlotabs ) {
        my $scmd1 = $scmd2;
        $scmd1 .= " $nlodir/$nlotab";
        $scmd1 .= " $nlotab";
        $scmd1 =~ s/\.${tabext}$/\.tab/;
        print "fnlo-add-tables.pl: Creating sum table for $nlotab ...\n";
        if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
        system("$scmd1 >> ${scen}_addst.log");
    }
# Normal mode: All LO with all NLO/NNLO tables
} else {
    my $scmd1 = $cmd;
    foreach my $lotab ( @lotabs ) {
        $scmd1 .= " $lodir/$lotab";
    }
    foreach my $nlotab ( @nlotabs ) {
        $scmd1 .= " $nlodir/$nlotab";
    }
    foreach my $nnlotab ( @nnlotabs ) {
        $scmd1 .= " $nnlodir/$nnlotab";
    }
    $scmd1 .= " $scen.tab";
    print "fnlo-add-tables.pl: Creating total sum table for $scen ...\n";
    if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
    system("$scmd1 > ${scen}_add.log");
}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-tables.pl: fastNLO finished: $date\n";
exit 0;
