#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 10.11.2011
#
# Create separate distribution and installation package for fastNLO_reader 
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.59
AC_PREREQ([2.59])
# Define subproject fastNLO reader
AC_INIT([fastNLO_reader], [2.1.0], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,stober@cern.ch,wobisch@fnal.gov])
AM_CONDITIONAL([HAVE_SUBPROJECT], test "x" = "x" )
# Require minimal automake version of 1.9, switch on warnings
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

# Announce checks particular for fastNLO_reader
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO_reader ...])
AC_MSG_NOTICE([******************************************])

# LHAPDF
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([lhapdf],
	[AS_HELP_STRING([--with-lhapdf],
	[Specify location of LHAPDF installation. This takes precedence over
	an LHAPDF environment variable, the prefix path or other methods.
	Please note: LHAPDF is mandatory in order to use fastNLO_reader!])],
	[if test "x$withval" = "xno"; then
	    AC_MSG_ERROR([fastNLO_reader CANNOT be used without LHAPDF! Configuration aborted!])
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
   AC_MSG_ERROR([fastNLO_reader CANNOT be used without LHAPDF! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_LHAPDF, test "x$LHAPDFPREFIX" != "x")
AM_CONDITIONAL(HAVE_LHAPDFCONFIG, test "x$LHAPDFBIN" != "x")

# Add subversion revision number, recipe from Autotools by J. Calcote 
SVNREV=`svnversion $srcdir | sed 's/:.*//'`
which svnversion > /dev/null; \
if [[ $? -ne 0 ] || [ "x$SVNREV" = "xexported" ]];
  then SVNREV=`cat $srcdir/SVNREV`
  else echo -n $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Output
AC_CONFIG_FILES([Makefile
                 reader_cc/Makefile
                 reader_f/Makefile
	        ])

AC_OUTPUT
