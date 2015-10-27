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
    print "  -n dir          Directory for NLO/NNLO or combined tables, (def.=scenario)\n";
    print "  -s              Produce tables for statistical evaluation,\n".
        "                  i.e. combinations of each LO with 1 NLO table and\n".
        "                  all LO with each NLO table\n".
        "                  Do not use with combined tables.\n";
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
my $alldir   = $nlodir;
my $tabext = "tab";
my $loglob   = "${scen}*born*.${tabext}*";
my $nloglob  = "${scen}*nlo*.${tabext}*";
my $nnloglob = "${scen}*thrcor*.${tabext}*";
my $allglob  = "${scen}*_d??-x??-y??_*.${tabext}*";

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
    print "fnlo-add-tables.pl: WARNING! No LO table found!\n";
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
    print "fnlo-add-tables.pl: WARNING! No NLO table found!\n";
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
# Find combined tables
#
my @alltabs;
if ( -d "$alldir" ) {
    chdir $alldir;
    @alltabs = glob $allglob;
    chdir $sdir;
}
unless ( @alltabs ) {
    print "fnlo-add-tables.pl: No combined table found in $alldir, now looking in $sdir ...\n";
    $alldir  = ".";
    @alltabs = glob $allglob;
}
unless ( @alltabs ) {
    print "fnlo-add-tables.pl: WARNING! No combined table found!\n";
}
$ntab = scalar @alltabs;
print "fnlo-add-tables.pl: $ntab combined tables found.\n";
if ( $debug ) {print "fnlo-add-tables.pl: DEBUG! alltabs @alltabs\n";}

unless ( @lotabs || @alltabs ) {
    die "fnlo-add-tables.pl: ERROR! Neither LO nor combined table found!\n";
}

#
# fnlo-tk-merge
#
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-tables.pl: ${merger}: $date\n";
# Statistics mode: Each LO with first NLO table & all LO with each NLO table
#                  Not necessary for combined tables.
if ( $opt_s && ! @alltabs ) {
    print "\nfnlo-add-tables.pl: Statistics mode\n";
    my $scmd2 = $cmd;
    print "\nfnlo-add-tables.pl: Creating LO statistics tables ...\n";
    foreach my $lotab ( @lotabs ) {
        if ( -f $lotab && ! -z $lotab ) {
            print "fnlo-add-tables.pl: Non-empty sum table for $lotab exists already, skipped!\n";
        } else {
            my $scmd1 = $cmd;
            $scmd1 .= " $lodir/$lotab $nlodir/$nlotabs[0]";
            $scmd1 .= " $lotab";
            $scmd1 =~ s/\.${tabext}$/\.tab/;
            print "fnlo-add-tables.pl: Creating sum table for $lotab ...\n";
            if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
            system("$scmd1 >> ${scen}_addst.log");
        }
        $scmd2 .= " $lodir/$lotab";
    }
    my $losum = $lotabs[0];
    $losum =~ s/\d\d\d\d\.tab/sum\.tab/;
    $scmd2 .= " $losum";
    $scmd2 =~ s/\.${tabext}$/\.tab/;
    if ( -f $losum && ! -z $losum ) {
        print "fnlo-add-tables.pl: Non-empty all LO sum table $losum exists already, skipped!\n";
    } else {
        print "\nfnlo-add-tables.pl: Creating temporary LO sum table $losum ...\n";
        if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd2\n";}
        system("$scmd2 >> ${scen}_addst.log");
    }
    print "\nfnlo-add-tables.pl: Creating NLO statistics tables ...\n";
    foreach my $nlotab ( @nlotabs ) {
        if ( -f $nlotab && ! -z $nlotab ) {
            print "fnlo-add-tables.pl: Non-empty sum table for $nlotab exists already, skipped!\n";
        } else {
            my $scmd1 = $cmd;
            $scmd1 .= " $losum";
            $scmd1 .= " $nlodir/$nlotab";
            $scmd1 .= " $nlotab";
            $scmd1 =~ s/\.${tabext}$/\.tab/;
            print "fnlo-add-tables.pl: Creating sum table for $nlotab ...\n";
            if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
            system("$scmd1 >> ${scen}_addst.log");
        }
    }
    $losum =~ s/\.${tabext}$/\.tab/;
    print "\nfnlo-add-tables.pl: Removing temporary LO sum table $losum ...\n";
    unlink $losum;
# Normal mode: All LO with all NLO/NNLO tables
} else {
    my $scentab = ${scen}.".tab";
    if ( -f $scentab && ! -z $scentab ) {
        print "fnlo-add-tables.pl: Non-empty sum table for $scentab exists already, skipped!\n";
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
        $scmd1 .= " $scentab";
        print "fnlo-add-tables.pl: Creating total sum table for $scen ...\n";
        if ( $debug ) {print "fnlo-add-tables.pl: Running command $scmd1\n";}
        system("$scmd1 > ${scen}_add.log");
    }
}
$date = `date +%d%m%Y_%H%M%S`;
chomp $date;
print "\nfnlo-add-tables.pl: fastNLO finished: $date\n";
exit 0;
