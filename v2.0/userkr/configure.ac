#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 08.06.2012
#
# Create separate distribution and installation package for my allunc
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.59
AC_PREREQ([2.59])
# Define subproject fastNLO_allunc
AC_INIT([fastNLO_allunc], [2.1.0], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,stober@cern.ch,wobisch@fnal.gov])
# Properly include subprojects
AC_CONFIG_AUX_DIR([.])
# Require minimal automake version of 1.9, switch on warnings
# Don't allow silent compilation/linking for now, incompatible with SLC5
#AM_INIT_AUTOMAKE([1.9 -Wall silent-rules])
AM_INIT_AUTOMAKE([1.9 -Wall])
# Require minimal libtool version of 1.5
#LT_PREREQ([1.5])
#LT_INIT([dlopen])
#
#AC_CONFIG_SRCDIR([common/fnloTable.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set default installation directory
AC_PREFIX_DEFAULT(${HOME}/local)

# Detect host system, required for conditional static CERNLIB linking
AC_CANONICAL_HOST
AM_CONDITIONAL(HAVE_64BIT, test "$host_cpu" = "x86_64")

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_F77([gfortran g77])
AC_PROG_FC([gfortran g77])
#case $FC in
#  gfortran) AC_SUBST([HAVE_GFORTRAN], [yes])
#  g77) AC_SUBST([HAVE_G77], [yes])
#esac
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
#AC_PROG_AWK
#AC_PROG_INSTALL
#AC_PROG_LN_S
#AC_PROG_MAKE_SET

# Checks for sources.

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS
AC_CHECK_LIB([g2c], [getenv_])
# FIXME: Replace `main' with a function in `-lLHAPDF' etc.
#AC_CHECK_LIB([LHAPDF], [LHAPDF])
#AC_CHECK_LIB([gfortranbegin], [fmain])
#AC_CHECK_LIB([lhpdf_cteq6], [cteq6])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt strstr])

# Announce checks particular for fastNLO_allunc
AC_MSG_NOTICE([*************************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO_allunc ...])
AC_MSG_NOTICE([*************************************************])

# LHAPDF
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([lhapdf],
	[AS_HELP_STRING([--with-lhapdf],
	[Specify location of LHAPDF installation. This takes precedence over
	an LHAPDF environment variable, the prefix path or other methods.
	Please note: LHAPDF is mandatory in order to use fastNLO_allunc!])],
	[if test "x$withval" = "xno"; then
	    AC_MSG_ERROR([fastNLO_allunc CANNOT be used without LHAPDF! Configuration aborted!])
	 fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([LHAPDF],[$withval])
	       AC_MSG_NOTICE([Using given path $LHAPDF as location of LHAPDF installation.])
	 fi],
	[])

if test "x$withval" != "xno"; then
   if test "x$LHAPDF" = "x"; then
      AC_MSG_NOTICE([No LHAPDF installation path set. Try to find LHAPDF in prefix path ...])
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$prefix/bin])
      if test "x$LHAPDFBIN" = "x"; then
         AC_MSG_NOTICE([No LHAPDF installation found in prefix path. Try to find LHAPDF in system path ...])
         AC_PATH_PROG([LHAPDFBIN],[lhapdf-config])
         if test "x$LHAPDFBIN" = "x"; then
            AC_MSG_ERROR([No LHAPDF library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$LHAPDF/bin])
      if test "x$LHAPDFBIN" = "x"; then
         AC_MSG_ERROR([lhapdf-config not found in given path to LHAPDF! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$LHAPDFBIN" != "x"; then
      AC_MSG_NOTICE([LHAPDF settings:])
      AC_SUBST([LHAPDFPREFIX],  [`$LHAPDFBIN --prefix`])
      AC_SUBST([LHAPDFSETPATH], [`$LHAPDFBIN --pdfsets-path`])
      AC_SUBST([LHAPDFBINPATH], ["$LHAPDFPREFIX/bin"])
      AC_SUBST([LHAPDFLIBPATH], [`$LHAPDFBIN --libdir`])
      AC_SUBST([LHAPDFINCLUDEPATH],[`$LHAPDFBIN --incdir`])
      AC_MSG_RESULT([Using $LHAPDFPREFIX for the LHAPDF installation path])
      AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])
      AC_MSG_RESULT([Using $LHAPDFBINPATH for the binaries])
      AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
      AC_MSG_RESULT([Using $LHAPDFINCLUDEPATH for the includes])
   else
      AC_MSG_ERROR([No LHAPDF library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_ERROR([fastNLO_allunc CANNOT be used without LHAPDF! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_LHAPDF, test "x$LHAPDFPREFIX" != "x")
AM_CONDITIONAL(HAVE_LHAPDFCONFIG, test "x$LHAPDFBIN" != "x")

# CERNLIB
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for CERNLIB installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([cernlib],
	[AS_HELP_STRING([--with-cernlib],
	[Specify location of CERNLIB installation. This takes precedence over
	a CERNLIB environment variable, the prefix path or other methods.
	Please note: CERNLIB is mandatory in order to use fastNLO_allunc!])],
	[if test "x$withval" = "xno"; then
	    AC_MSG_ERROR([fastNLO_allunc CANNOT be used without CERNLIB! Configuration aborted!])
	 fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([CERNLIB],[$withval])
	       AC_MSG_NOTICE([Using given path $CERNLIB as location of CERNLIB installation.])
	 fi],
	[])

if test "x$withval" != "xno"; then
   if test "x$CERNLIB" = "x"; then
      AC_MSG_NOTICE([No CERNLIB installation path set. Try to find CERNLIB in prefix path ...])
      AC_PATH_PROG([CERNLIBBIN],[cernlib],[],[$prefix/bin])
      if test "x$CERNLIBBIN" = "x"; then
         AC_MSG_NOTICE([No CERNLIB installation found in prefix path. Try to find CERNLIB in system path ...])
         AC_PATH_PROG([CERNLIBBIN],[cernlib])
         if test "x$CERNLIBBIN" = "x"; then
            AC_MSG_ERROR([No CERNLIB library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([CERNLIBBIN],[cernlib],[],[$CERNLIB/bin])
      if test "x$CERNLIBBIN" = "x"; then
         AC_MSG_NOTICE([cernlib command not found in given path to CERNLIB, pretend to have found it anyway!])
         AC_MSG_NOTICE([Keep your fingers crossed ...])
         AC_SUBST([CERNLIBBIN], ["$CERNLIB/bin/cernlib"])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$CERNLIBBIN" != "x"; then
      AC_MSG_NOTICE([CERNLIB settings:])
      AC_SUBST([CERNLIBCMD], ["$CERNLIBBIN"])
      AC_SUBST([CERNLIBBINPATH], [`dirname $CERNLIBBIN`])
      AC_SUBST([CERNLIBPREFIX], [`dirname $CERNLIBBINPATH`])
      AC_SUBST([CERNLIBLIBPATH], ["$CERNLIBPREFIX/lib"])
      AC_MSG_RESULT([Using $CERNLIBCMD for the cernlib command])
      AC_MSG_RESULT([Using $CERNLIBPREFIX for the CERNLIB installation path])
      AC_MSG_RESULT([Using $CERNLIBBINPATH for the binaries])
      AC_MSG_RESULT([Using $CERNLIBLIBPATH for the libraries])
   else
      AC_MSG_ERROR([No CERNLIB library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_ERROR([fastNLO_allunc CANNOT be used without CERNLIB! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_CERNLIB, test "x$CERNLIBLIBPATH" != "x")
AM_CONDITIONAL(HAVE_CERNLIBCONFIG, test "x$CERNLIBCMD" != "x")

# Add subversion revision number, recipe from Autotools by J. Calcote
SVNREV=`svnversion $srcdir | sed 's/:.*//'`
which svnversion > /dev/null; \
if [[ $? -ne 0 ] || [ "x$SVNREV" = "xexported" ]];
  then SVNREV=`cat $srcdir/SVNREV`
  else echo -n $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Define project substitution variables (to be replaced e.g. in config files below)
AC_SUBST([AUTHORSv14],["T. Kluge, K. Rabbertz, M. Wobisch"])
AC_SUBST([AUTHORSv2],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"])
AC_SUBST([QUOTEv14],["hep-ph/0609285"])
AC_SUBST([QUOTEv2],["arXiv:1109.1310"])
AC_SUBST([WEBPAGE],["http://projects.hepforge.org/fastnlo"])

# Output
AC_CONFIG_FILES([Makefile
                 main.f
	        ])

AC_OUTPUT
