#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 15.05.2014
#
# Create distribution and installation package for
# NLOJet++ interface to fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.59
AC_PREREQ([2.59])
# Define subproject fastNLO_interface_nlojet
AC_INIT([fastNLO_interface_nlojet], [2.3.1pre], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,g.sieber@cern.ch,stober@cern.ch,wobisch@fnal.gov])
# Properly include subprojects
AC_CONFIG_AUX_DIR([.])
# Require minimal automake version of 1.9, switch on warnings
AM_INIT_AUTOMAKE([1.9 -Wall])
# Require minimal libtool version of 1.5
#LT_PREREQ([1.5])
#LT_INIT([dlopen])
#
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set default installation directory
AC_PREFIX_DEFAULT(${HOME}/local)

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_F77([gfortran g77])
AC_PROG_FC([gfortran g77])
#case $FC in
#  gfortran) AC_SUBST([HAVE_GFORTRAN], [yes])
#  g77) AC_SUBST([HAVE_G77], [yes])
#esac
AC_PROG_INSTALL
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET
# AC_PROG_RANLIB

# Checks for sources.

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS
AC_CHECK_LIB([g2c], [getenv_])
AC_CHECK_LIB([LHAPDF], [alphas])
#AC_CHECK_LIB([fastnlotoolkit], [libfastnlotoolkit_la-FastNLOBlockB])
#AC_CHECK_LIB([gfortranbegin], [fmain])
#AC_CHECK_LIB([lhpdf_cteq6], [cteq6])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([memset pow sqrt strstr])

# Announce checks particular for fastNLO_interface_nlojet
AC_MSG_NOTICE([*************************************************])
AC_MSG_NOTICE([Starting checks particular for $PACKAGE_NAME ...])
AC_MSG_NOTICE([*************************************************])

# NLOJet++
AC_SUBST([TMPNAM],["NLOJet++"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for $TMPNAM installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([nlojet],
        [AS_HELP_STRING([--with-nlojet],
        [Specify location of $TMPNAM installation. This takes precedence over
        an NLOJET environment variable, the prefix path or other methods.
        Please note: $TMPNAM is mandatory in order to use $PACKAGE_NAME!])],
        [if test "x$withval" = "xno"; then
            AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
         fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([NLOJET],[$withval])
               AC_MSG_NOTICE([Using given path $NLOJET as location of $TMPNAM installation.])
         fi],
        [])

if test "x$with_nlojet" != "xno"; then
   if test "x$NLOJET" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([NLOJETBIN],[nlojet++],[],[$prefix/bin])
      if test "x$NLOJETBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([NLOJETBIN],[nlojet++])
         if test "x$NLOJETBIN" = "x"; then
            AC_MSG_ERROR([No $TMPNAM installation found. Please install $TMPNAM first. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([NLOJETBIN],[nlojet++],[],[$NLOJET/bin])
      if test "x$NLOJETBIN" = "x"; then
         AC_MSG_ERROR([nlojet++ not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$NLOJETBIN" != "x"; then
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([NLOJET], [`dirname $NLOJETBIN`])
      AC_SUBST([NLOJETPREFIX], [`dirname $NLOJET`])
      AC_SUBST([NLOJETBINPATH], ["$NLOJETPREFIX/bin"])
      AC_SUBST([NLOJETLIBPATH], ["$NLOJETPREFIX/lib"])
      AC_SUBST([NLOJETINCLUDEPATH],["$NLOJETPREFIX/include/nlo++"])
      AC_MSG_RESULT([Using $NLOJETPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $NLOJETBINPATH for the binaries])
      AC_MSG_RESULT([Using $NLOJETLIBPATH for the library])
      AC_MSG_RESULT([Using $NLOJETINCLUDEPATH for the includes])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_NLOJET, test "x$NLOJETPREFIX" != "x")

# fastNLO toolkit
AC_SUBST([TMPNAM],["fastNLO toolkit"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for $TMPNAM installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([fnlotoolkit],
        [AS_HELP_STRING([--with-fnlotoolkit],
        [Specify location of $TMPNAM installation. This takes precedence over
         an FNLOTOOLKIT environment variable, the prefix path or other methods.
         Please note: $TMPNAM is mandatory in order to use $PACKAGE_NAME!])],
        [if test "x$withval" = "xno"; then
            AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
         fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([FNLOTOOLKIT],[$withval])
               AC_MSG_NOTICE([Using given path $FNLOTOOLKIT as location of $TMPNAM installation.])
         fi],
        [])

if test "x$with_fnlotoolkit" != "xno"; then
   if test "x$FNLOTOOLKIT" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([FNLOTOOLKITBIN],[fnlo-tk-config],[],[$prefix/bin])
      if test "x$FNLOTOOLKITBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([FNLOTOOLKITBIN],[fnlo-tk-config])
         if test "x$FNLOTOOLKITBIN" = "x"; then
            AC_MSG_ERROR([No $TMPNAM installation found. Please install $TMPNAM first. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([FNLOTOOLKITBIN],[fnlo-tk-config],[],[$FNLOTOOLKIT/bin])
      if test "x$FNLOTOOLKITBIN" = "x"; then
         AC_MSG_ERROR([fnlo-tk-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$FNLOTOOLKITBIN" != "x"; then
      AC_MSG_NOTICE([$TMPNAM settings:])
#      AC_SUBST([FNLOTOOLKIT], [`dirname $FNLOTOOLKITBIN`])
#      AC_SUBST([FNLOTOOLKITPREFIX], [`dirname $FNLOTOOLKIT`])
      AC_SUBST([FNLOTOOLKITPREFIX],  [`$FNLOTOOLKITBIN --prefix`])
      AC_SUBST([FNLOTOOLKITBINPATH], ["$FNLOTOOLKITPREFIX/bin"])
      AC_SUBST([FNLOTOOLKITLIBPATH], [`$FNLOTOOLKITBIN --libdir`])
      AC_SUBST([FNLOTOOLKITINCLUDEPATH],[`$FNLOTOOLKITBIN --incdir`])
      AC_MSG_RESULT([Using $FNLOTOOLKITPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $FNLOTOOLKITBINPATH for the binaries])
      AC_MSG_RESULT([Using $FNLOTOOLKITLIBPATH for the library])
      AC_MSG_RESULT([Using $FNLOTOOLKITINCLUDEPATH for the includes])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_FNLOTOOLKIT, test "x$FNLOTOOLKITPREFIX" != "x")

# FastJet
AC_SUBST([TMPNAM],["FastJet"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for $TMPNAM installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([fastjet],
        [AS_HELP_STRING([--with-fastjet],
         [Specify location of $TMPNAM installation. This takes precedence over
          a FASTJET environment variable, the prefix path or other methods.
          Please note: $TMPNAM is not required but warmly recommended in order to apply jet algorithms!])],
        [if test "x$withval" = "xno"; then
            AC_MSG_NOTICE([Building $PACKAGE_NAME without $TMPNAM!])
         else
            if test "x$withval" != "xyes"; then
               AC_SUBST([FASTJET],[$withval])
               AC_MSG_NOTICE([Using given path $FASTJET as location of $TMPNAM installation.])
            fi
         fi],
        [])

if test "x$withval" != "xno"; then
   if test "x$FASTJET" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([FASTJETBIN],[fastjet-config],[],[$prefix/bin])
      if test "x$FASTJETBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([FASTJETBIN],[fastjet-config])
         if test "x$FASTJETBIN" = "x"; then
            AC_MSG_NOTICE([Building $PACKAGE_NAME without $TMPNAM!])
         fi
      fi
   else
      AC_PATH_PROG([FASTJETBIN],[fastjet-config],[],[$FASTJET/bin])
      if test "x$FASTJETBIN" = "x"; then
         AC_MSG_ERROR([fastjet-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$FASTJETBIN" != "x"; then
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([FASTJETPREFIX],  [`$FASTJETBIN --prefix`])
      AC_SUBST([FASTJETBINPATH], ["$FASTJETPREFIX/bin"])
      AC_SUBST([FASTJETLIBPATH], ["$FASTJETPREFIX/lib"])
      AC_SUBST([FASTJETINCLUDEPATH],["$FASTJETPREFIX/include"])
      AC_SUBST([FASTJETLIBS],[`$FASTJETBINPATH/fastjet-config --libs`])
      AC_MSG_RESULT([Using $FASTJETPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $FASTJETBINPATH for the binaries])
      AC_MSG_RESULT([Using $FASTJETLIBPATH for the library])
      AC_MSG_RESULT([Using $FASTJETINCLUDEPATH for the includes])
      AC_MSG_RESULT([Using $FASTJETLIBS for the linker flags])
   else
      AC_MSG_NOTICE([No $TMPNAM library found. Building $PACKAGE_NAME without $TMPNAM!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([$TMPNAM support switched off by user. Building $PACKAGE_NAME without $TMPNAM!])
fi
AM_CONDITIONAL(HAVE_FASTJET, test "x$FASTJETPREFIX" != "x")
AM_CONDITIONAL(HAVE_FASTJETCONFIG, test "x$FASTJETBIN" != "x")

# Add subversion revision number, recipe from Autotools by J. Calcote
# Modified on 18.12.13 to fix changed behaviour of svnversion to give "Unversioned directory" instead of "exported"
SVNREV=`( svnversion $srcdir | sed 's/:.*//' ) 2>/dev/null`
if { ! ( svnversion ) >/dev/null 2>&1 || test "x$SVNREV" = "xexported" || test "x$SVNREV" = "xUnversioned directory"; } ;
  then SVNREV=`cat $srcdir/SVNREV`
  else echo $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Define project substitution variables (to be replaced e.g. in config files below)
AC_SUBST([AUTHORS],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, G. Sieber, M. Wobisch"])
AC_SUBST([AUTHORSv14],["T. Kluge, K. Rabbertz, M. Wobisch"])
AC_SUBST([AUTHORSv2],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"])
AC_SUBST([QUOTEv14],["hep-ph/0609285"])
AC_SUBST([QUOTEv2],["arXiv:1109.1310"])
AC_SUBST([WEBPAGE],["http://projects.hepforge.org/fastnlo"])
AC_SUBST([YEARS],["2011-`date +%Y`"])

# Output
AC_CONFIG_FILES([Makefile
                 interface/Makefile
                 interface/check/Makefile
                 interface/data/Makefile
                 interface/data/Examples/Makefile
                 interface/data/InclusiveJets/Makefile
                 interface/data/DijetMass/Makefile
                 interface/data/DijetAngular/Makefile
                 interface/data/AzimuthalDecorrelation/Makefile
                 interface/data/InclusiveNJetEvents/Makefile
                 interface/data/RivetAdditions/Makefile
                 interface/data/check/Makefile
                 interface/hadron/Makefile
                 interface/include/Makefile
                 interface/include/fnlo_int_nlojet/Makefile
                 interface/tools/Makefile
                 interface/jetalgos/Makefile
                ])

AC_OUTPUT
