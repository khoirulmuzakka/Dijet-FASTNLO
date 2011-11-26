#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 26.11.2011
#
# Create separate distribution and installation package for fastNLO_converter 
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.59
AC_PREREQ([2.59])
# Define subproject fastNLO converter
AC_INIT([fastNLO_converter], [2.1.0], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,stober@cern.ch,wobisch@fnal.gov])
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

# Announce checks particular for fastNLO_converter
AC_MSG_NOTICE([*******************************************************])
AC_MSG_NOTICE([No particular checks required for fastNLO_converter ...])
AC_MSG_NOTICE([*******************************************************])

# Add subversion revision number, recipe from Autotools by J. Calcote 
SVNREV=`svnversion $srcdir | sed 's/:.*//'`
which svnversion > /dev/null; \
if [[ $? -ne 0 ] || [ "x$SVNREV" = "xexported" ]];
  then SVNREV=`cat $srcdir/SVNREV`
  else echo -n $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Output
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
