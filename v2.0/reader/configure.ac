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
# Define subproject fastNLO_reader
AC_INIT([fastNLO_reader], [2.1.0], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,stober@cern.ch,wobisch@fnal.gov])
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

# Use Doxygen if available (doesn't work with autotools from SLC5)
#DX_HTML_FEATURE(ON)
#DX_CHM_FEATURE(OFF)
#DX_CHI_FEATURE(OFF)
#DX_MAN_FEATURE(OFF)
#DX_RTF_FEATURE(OFF)
#DX_XML_FEATURE(OFF)
#DX_PDF_FEATURE(OFF)
#DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(fastNLO_reader, doxygen.cfg, doc/doxygen)

# Checks for sources.

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS
AC_CHECK_LIB([g2c], [getenv_])
AC_CHECK_LIB([LHAPDF], [alphas])
AC_CHECK_LIB([fastnloreader], [libfastnloreader_la-FastNLOBlockB])
#AC_CHECK_LIB([gfortranbegin], [fmain])
#AC_CHECK_LIB([lhpdf_cteq6], [cteq6])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([memset pow sqrt strstr])

# Announce checks particular for fastNLO_reader
AC_MSG_NOTICE([*************************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO_reader ...])
AC_MSG_NOTICE([*************************************************])

# LHAPDF
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([lhapdf],
        [AS_HELP_STRING([--with-lhapdf],
        [Specify location of LHAPDF installation. This takes precedence over
        an LHAPDF environment variable, the prefix path or other methods.
        Please note: LHAPDF is mandatory in order to use fastNLO_reader!])],
        [if test "x$withval" == "xno"; then
            AC_MSG_ERROR([fastNLO_reader CANNOT be used without LHAPDF! Configuration aborted!])
         fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([LHAPDF],[$withval])
               AC_MSG_NOTICE([Using given path $LHAPDF as location of LHAPDF installation.])
         fi],
        [])

if test "x$with_lhapdf" != "xno"; then
   if test "x$LHAPDF" == "x"; then
      AC_MSG_NOTICE([No LHAPDF installation path set. Try to find LHAPDF in prefix path ...])
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$prefix/bin])
      if test "x$LHAPDFBIN" == "x"; then
         AC_MSG_NOTICE([No LHAPDF installation found in prefix path. Try to find LHAPDF in system path ...])
         AC_PATH_PROG([LHAPDFBIN],[lhapdf-config])
         if test "x$LHAPDFBIN" == "x"; then
            AC_MSG_ERROR([No LHAPDF library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$LHAPDF/bin])
      if test "x$LHAPDFBIN" == "x"; then
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

# QCDNUM
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether QCDNUM support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([qcdnum],
        [AS_HELP_STRING([--with-qcdnum],
        [Specify location of QCDNUM installation. This takes precedence over
        a QCDNUM environment variable, the prefix path or other methods.
        Please note: QCDNUM can optionally be used in fastNLO_reader!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([QCDNUM],[$withval])
            AC_MSG_NOTICE([Using given path $QCDNUM as location of QCDNUM installation.])
        fi],
        [with_qcdnum="no"])

if test "x$with_qcdnum" != "xno"; then
   if test "x$QCDNUM" == "x"; then
      AC_MSG_NOTICE([No QCDNUM installation path set. Try to find QCDNUM in prefix path ...])
      AC_CHECK_FILE([$prefix/lib/libqcdnum.a],
      [preqcdnum="$ac_cv_file_$prefix/lib/libqcdnum.a"],
      [preqcdnum=""])
      if test "x$preqcdnum" == "x"; then
         AC_MSG_NOTICE([No QCDNUM library found in prefix path. Try to link with systems QCDNUM ...])
         AC_CHECK_LIB([qcdnum],[libqcdnum_a-qcdnum],
         [sysqcdnum=$ac_cv_lib_libqcdnum_a-qcdnum],
         [sysqcdnum=""])
         if test "x$sysqcdnum" == "x"; then
            AC_MSG_ERROR([No QCDNUM library found. Please install QCDNUM first. Configuration aborted!])
         else
            AC_SUBST([QCDNUMLIBPATH], [`dirname $sysqcdnum`])
         fi
      else
         AC_SUBST([QCDNUMLIBPATH], [`dirname $preqcdnum`])
      fi
      AC_SUBST([QCDNUM], [`dirname $QCDNUMLIBPATH`])
   else
      AC_CHECK_FILE([$QCDNUM/lib/libqcdnum.a],
      [],
      [AC_MSG_ERROR([QCDNUM library not found at given path $QCDNUM/lib. Please install QCDNUM first. Configuration aborted!])])
   fi
   AC_SUBST([QCDNUMPREFIX], "$QCDNUM")
   AC_SUBST([QCDNUMLIBPATH], "$QCDNUMPREFIX/lib")
   AC_SUBST([QCDNUMLIBS], ["-lzmstf -lhqstf -lqcdnum -lmbutil"])
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([QCDNUM settings:])
   AC_MSG_RESULT([Using $QCDNUMPREFIX for the QCDNUM installation path])
   AC_MSG_RESULT([Using $QCDNUMLIBPATH for the libraries])
else
   AC_MSG_NOTICE([No QCDNUM support requested. Continue without QCDNUM.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_QCDNUM, test "x$QCDNUMPREFIX" != "x")

# Hoppet
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether Hoppet support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([hoppet],
        [AS_HELP_STRING([--with-hoppet],
        [Specify location of Hoppet installation. This takes precedence over
        a HOPPET environment variable, the prefix path or other methods.
        Please note: Hoppet can optionally be used in fastNLO_reader!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([HOPPET],[$withval])
            AC_MSG_NOTICE([Using given path $HOPPET as location of HOPPET installation.])
        fi],
        [with_hoppet="no"])

if test "x$with_hoppet" != "xno"; then
   if test "x$HOPPET" == "x"; then
      AC_MSG_NOTICE([No Hoppet installation path set. Try to find Hoppet in prefix path ...])
      AC_CHECK_FILE([$prefix/lib/libhoppet_v1.a],
      [prehoppet="$ac_cv_file_$prefix/lib/libhoppet_v1.a"],
      [prehoppet=""])
      if test "x$prehoppet" == "x"; then
         AC_MSG_NOTICE([No Hoppet library found in prefix path. Try to link with systems Hoppet ...])
         AC_CHECK_LIB([libhoppet_v1],[hoppetstart],
         [syshoppet=$ac_cv_lib_libhoppet_v1_hoppetstart],
         [syshoppet=""])
         if test "x$syshoppet" == "x"; then
            AC_MSG_ERROR([No Hoppet library found. Please install Hoppet first. Configuration aborted!])
         else
            AC_SUBST([HOPPETLIBPATH], [`dirname $syshoppet`])
         fi
      else
         AC_SUBST([HOPPETLIBPATH], [`dirname $prehoppet`])
      fi
      AC_SUBST([HOPPET], [`dirname $HOPPETLIBPATH`])
   else
      AC_CHECK_FILE([$HOPPET/lib/libhoppet_v1.a],
      [],
      [AC_MSG_ERROR([HOPPET library not found at given path $HOPPET/lib. Please install Hoppet first. Configuration aborted!])])
   fi
   AC_SUBST([HOPPETPREFIX], ["$HOPPET"])
   AC_SUBST([HOPPETBIN], ["$HOPPETPREFIX/bin"])
   AC_SUBST([HOPPETLDFLAGS], ["$($HOPPETBIN/hoppet-config --libs)"])
   AC_SUBST([HOPPETCPPFLAGS], ["$($HOPPETBIN/hoppet-config --cxxflags)"])
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([Hoppet settings:])
   AC_MSG_RESULT([Using $HOPPETPREFIX for the Hoppet installation path])
   AC_MSG_RESULT([Using $HOPPETBIN for the Hoppet bin path])
   AC_MSG_RESULT([Using $HOPPETLDFLAGS for the libraries])
else
   AC_MSG_NOTICE([No Hoppet support requested. Continue without Hoppet.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_HOPPET, test "x$HOPPETPREFIX" != "x")

# Python extension
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether Python extension is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_ENABLE([pyext],
        [AC_HELP_STRING([--enable-pyext],
        [Optionally build Python module to interface C++ library. Default is no.])],
        [],
        [enable_pyext="no"])
# Basic Python checks
if test "x$enable_pyext" == "xyes"; then
   AX_PYTHON_DEVEL([>='2.3.1'])
   AM_PATH_PYTHON
fi
# Check on existence of SWIG to interface to Python
if test "x$enable_pyext" == "xyes" -a ! -e pyext/FastNLOReader_wrap.cpp; then
   AC_PATH_PROG(SWIG, swig)
   if test "x$SWIG" ==  "x"; then
      AC_MSG_ERROR([swig not found, interfacing to Python not possible! Please install swig first or disable Python extension. Configuration aborted!])
   else
      AC_MSG_NOTICE([swig found, interface to Python possible.])
   fi
fi
AC_MSG_NOTICE([******************************************])
if test "x$enable_pyext" == "xyes"; then
   AC_MSG_NOTICE([All Python checks successful. Requested Python extension will be built.])
else
   AC_MSG_NOTICE([No Python extension requested. Continue without pyext.])
fi
AC_MSG_NOTICE([******************************************])
AM_CONDITIONAL(ENABLE_PYEXT, [test "x$enable_pyext" == "xyes"])

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
                 reader_cc/Makefile
                 reader_cc/fastnloreader/Makefile
                 reader_cc/fastnloreader/FastNLOReaderVersion.cc
                 reader_cc/fastnloreader/include/Makefile
                 reader_cc/fastnloreader/include/fastnlo/Makefile
                 reader_cc/pyext/Makefile
                 reader_cc/pyext/setup.py
                 reader_cc/src/Makefile
                 reader_f/Makefile
                 reader_f/fnlo-fread.f
                ])

AC_OUTPUT
