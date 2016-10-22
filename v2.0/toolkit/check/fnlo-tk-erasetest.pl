#!/usr/bin/env perl
#
# Unit test to erase bins from a fastNLO table
# Version:
#
# created by K. Rabbertz: 20.10.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

# Remove potentially left-over temporary files
my @tabs = ("InclusiveNJetEvents_fnl5662i_v23_fix_2-2-2", "InclusiveNJetEvents_fnl5662i_v23_fix_2-2");
my @tabfls;
my @tabgzs;
foreach my $tab ( @tabs ) {
    push @tabfls, ${tab}.".tab";
    push @tabgzs, ${tab}.".tab".".gz";
}
foreach my $file ( "erasetest.tab", "erasediff.log", @tabfls, @tabgzs ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
foreach my $tabgz ( @tabgzs ) {
    my $cmd = "cp -f ../data/check/${tabgz} .";
    print "Executing command: $cmd\n";
    my $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-erasetest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
    $cmd = "gunzip ${tabgz}";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-erasetest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}
}

# Erase
my $cmd = "fnlo-tk-modify steerfile=../data/check/SteerModify_EraseBins.str";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-erasetest: Erasing bins from test table $tabfls[0] failed: $ret, aborted!\n";}

# Determine difference to default table with bins removed
$cmd = "diff $tabfls[1] erasetest.tab > erasediff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-erasetest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "erasediff.log" ) {
    print "fnlo-tk-erasetest: Test table with bins erased differs from default:\n";
    $cmd = "cat erasediff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-erasetest: Modification unit test failed, please fix!\n";
}

print "fnlo-tk-erasetest: Modification unit test passed.\n";

exit 0;
