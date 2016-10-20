#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.10.2013
#
# Create distribution and installation package for fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.62
AC_PREREQ([2.62])
# Define subproject fastNLO_toolkit
AC_INIT([fastNLO_toolkit], [2.3.1pre], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,g.sieber@cern.ch,stober@cern.ch,wobisch@fnal.gov])
# Properly include subprojects
AC_CONFIG_AUX_DIR([.])
# Require minimal automake version of 1.9, switch on warnings
AM_INIT_AUTOMAKE([1.9 -Wall])
# Require minimal libtool version of 1.5
#LT_PREREQ([1.5])
#LT_INIT([dlopen])
#
# AC_DEFINE preprocessor definitions are written to config.h
# which must be included first into all respective source files
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

# Announce checks particular for fastNLO_toolkit
AC_MSG_NOTICE([*************************************************])
AC_MSG_NOTICE([Starting checks particular for $PACKAGE_NAME ...])
AC_MSG_NOTICE([*************************************************])

# LHAPDF
AC_SUBST([TMPNAM],["LHAPDF"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for $TMPNAM installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([lhapdf],
        [AS_HELP_STRING([--with-lhapdf],
        [Specify location of $TMPNAM installation. This takes precedence over
        an LHAPDF environment variable, the prefix path or other methods.
        Please note: $TMPNAM is mandatory in order to use $PACKAGE_NAME!])],
        [if test "x$withval" = "xno"; then
            AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
         fi
         if test "x$withval" != "xyes"; then
               AC_SUBST([LHAPDF],[$withval])
               AC_MSG_NOTICE([Using given path $LHAPDF as location of $TMPNAM installation.])
         fi],
        [])

if test "x$with_lhapdf" != "xno"; then
   if test "x$LHAPDF" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path $prefix/bin ...])
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$prefix/bin])
      if test "x$LHAPDFBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([LHAPDFBIN],[lhapdf-config])
         if test "x$LHAPDFBIN" = "x"; then
            AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$LHAPDF/bin])
      if test "x$LHAPDFBIN" = "x"; then
         AC_MSG_ERROR([lhapdf-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$LHAPDFBIN" != "x"; then
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([LHAPDFPREFIX],  [`$LHAPDFBIN --prefix`])
      AC_SUBST([LHAPDFVERSION], [`$LHAPDFBIN --version`])
      LHAPDFVMIN="6"
      LHAPDFVORD=`echo -e "$LHAPDFVERSION\n$LHAPDFVMIN" | sort -t"." -g | head -n1`
      if test "x$LHAPDFVORD" = "x$LHAPDFVMIN"; then
         AC_SUBST([LHAPDFSETPATH], [`$LHAPDFBIN --datadir`])
      else
         AC_SUBST([LHAPDFSETPATH], [`$LHAPDFBIN --pdfsets-path`])
      fi
      AC_SUBST([LHAPDFBINPATH], ["$LHAPDFPREFIX/bin"])
      AC_SUBST([LHAPDFLIBPATH], [`$LHAPDFBIN --libdir`])
      AC_SUBST([LHAPDFINCLUDEPATH],[`$LHAPDFBIN --incdir`])
      AC_MSG_RESULT([Using $TMPNAM version $LHAPDFVERSION])
      AC_MSG_RESULT([Using $LHAPDFPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])
      AC_MSG_RESULT([Using $LHAPDFBINPATH for the binaries])
      AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
      AC_MSG_RESULT([Using $LHAPDFINCLUDEPATH for the includes])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_ERROR([$PACKAGE_NAME CANNOT be used without $TMPNAM! Configuration aborted!])
fi
AM_CONDITIONAL(HAVE_LHAPDF, test "x$LHAPDFPREFIX" != "x")
AM_CONDITIONAL(HAVE_LHAPDFCONFIG, test "x$LHAPDFBIN" != "x")

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
      AC_SUBST([FASTJETVERSION], [`$FASTJETBIN --version`])
      AC_SUBST([FASTJETBINPATH], ["$FASTJETPREFIX/bin"])
      AC_SUBST([FASTJETLIBPATH], ["$FASTJETPREFIX/lib"])
      AC_SUBST([FASTJETINCLUDEPATH],["$FASTJETPREFIX/include"])
      AC_SUBST([FASTJETLIBS],[`$FASTJETBINPATH/fastjet-config --libs`])
      AC_MSG_RESULT([Using $TMPNAM version $FASTJETVERSION])
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

# QCDNUM
AC_SUBST([TMPNAM],["QCDNUM"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether $TMPNAM support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([qcdnum],
        [AS_HELP_STRING([--with-qcdnum],
        [Specify location of $TMPNAM installation. This takes precedence over
        a QCDNUM environment variable, the prefix path or other methods.
        Please note: $TMPNAM can optionally be used in $PACKAGE_NAME!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([QCDNUM],[$withval])
            AC_MSG_NOTICE([Using given path $QCDNUM as location of $TMPNAM installation.])
        fi],
        [with_qcdnum="no"])

if test "x$with_qcdnum" != "xno"; then
   if test "x$QCDNUM" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([QCDNUMBIN],[qcdnum-config],[],[$prefix/bin])
      if test "x$QCDNUMBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([QCDNUMBIN],[qcdnum-config])
         if test "x$QCDNUMBIN" = "x"; then
            AC_MSG_NOTICE([Building $PACKAGE_NAME without $TMPNAM!])
         fi
      fi
   else
      AC_PATH_PROG([QCDNUMBIN],[qcdnum-config],[],[$QCDNUM/bin])
      if test "x$QCDNUMBIN" = "x"; then
         AC_MSG_ERROR([qcdnum-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$QCDNUMBIN" != "x"; then
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([QCDNUM],        [`$QCDNUMBIN --prefix`])
      AC_SUBST([QCDNUMPREFIX],  [`$QCDNUMBIN --prefix`])
      AC_SUBST([QCDNUMVERSION], [`$QCDNUMBIN --version`])
      AC_SUBST([QCDNUMBINPATH], ["$QCDNUMPREFIX/bin"])
      AC_SUBST([QCDNUMLIBPATH], [`$QCDNUMBIN --libdir`])
      AC_SUBST([QCDNUMINCLUDEPATH],["$QCDNUMPREFIX/include"])
      AC_SUBST([QCDNUMLDFLAGS], [`$QCDNUMBIN --ldflags`])
      AC_MSG_RESULT([Using $TMPNAM version $QCDNUMVERSION])
      AC_MSG_RESULT([Using $QCDNUMPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $QCDNUMBINPATH for the binaries])
      AC_MSG_RESULT([Using $QCDNUMLIBPATH for the libraries])
      AC_MSG_RESULT([Using $QCDNUMINCLUDEPATH for the includes])
      AC_MSG_RESULT([Using $QCDNUMLDFLAGS for the linker flags])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([No $TMPNAM support requested. Continue without $TMPNAM.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_QCDNUM, test "x$QCDNUMPREFIX" != "x")

# HOPPET
AC_SUBST([TMPNAM],["HOPPET"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether $TMPNAM support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([hoppet],
        [AS_HELP_STRING([--with-hoppet],
        [Specify location of $TMPNAM installation. This takes precedence over
        a HOPPET environment variable, the prefix path or other methods.
        Please note: $TMPNAM can optionally be used in $PACKAGE_NAME!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([HOPPET],[$withval])
            AC_MSG_NOTICE([Using given path $HOPPET as location of $TMPNAM installation.])
        fi],
        [with_hoppet="no"])

if test "x$with_hoppet" != "xno"; then
   if test "x$HOPPET" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_CHECK_FILE([$prefix/lib/libhoppet_v1.a],
      [prehoppet="$ac_cv_file_$prefix/lib/libhoppet_v1.a"],
      [prehoppet=""])
      if test "x$prehoppet" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM library found in prefix path. Try to link with systems $TMPNAM ...])
         AC_CHECK_LIB([libhoppet_v1],[hoppetstart],
         [syshoppet=$ac_cv_lib_libhoppet_v1_hoppetstart],
         [syshoppet=""])
         if test "x$syshoppet" = "x"; then
            AC_MSG_ERROR([No $TMPNAM library found. Please install $TMPNAM first. Configuration aborted!])
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
      [AC_MSG_ERROR([$TMPNAM library not found at given path $HOPPET/lib. Please install $TMPNAM first. Configuration aborted!])])
   fi
   AC_SUBST([HOPPETPREFIX], ["$HOPPET"])
   AC_SUBST([HOPPETBIN], ["$HOPPETPREFIX/bin"])
   AC_SUBST([HOPPETLDFLAGS], ["$($HOPPETBIN/hoppet-config --libs)"])
   AC_SUBST([HOPPETCXXFLAGS], ["$($HOPPETBIN/hoppet-config --cxxflags)"])
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([$TMPNAM settings:])
   AC_MSG_RESULT([Using $HOPPETPREFIX for the $TMPNAM installation path])
   AC_MSG_RESULT([Using $HOPPETBIN for the $TMPNAM bin path])
   AC_MSG_RESULT([Using $HOPPETCXXFLAGS for the compilation])
   AC_MSG_RESULT([Using $HOPPETLDFLAGS for the libraries])
else
   AC_MSG_NOTICE([No $TMPNAM support requested. Continue without $TMPNAM.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_HOPPET, test "x$HOPPETPREFIX" != "x")

# YODA
AC_SUBST([TMPNAM],["YODA"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether $TMPNAM support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([yoda],
        [AS_HELP_STRING([--with-yoda],
        [Specify location of $TMPNAM installation. This takes precedence over
        a YODA environment variable, the prefix path or other methods.
        Please note: $TMPNAM can optionally be used in $PACKAGE_NAME!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([YODA],[$withval])
            AC_MSG_NOTICE([Using given path $YODA as location of $TMPNAM installation.])
        fi],
        [with_yoda="no"])

if test "x$with_yoda" != "xno"; then
   if test "x$YODA" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([YODABIN],[yoda-config],[],[$prefix/bin])
      if test "x$YODABIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([YODABIN],[yoda-config])
         if test "x$YODABIN" = "x"; then
            AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([YODABIN],[yoda-config],[],[$YODA/bin])
      if test "x$YODABIN" = "x"; then
         AC_MSG_ERROR([yoda-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$YODABIN" != "x"; then
      AC_DEFINE([WITH_YODA],[1],[Enable conditional compilation with YODA])
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([YODAPREFIX],  [`$YODABIN --prefix`])
      AC_SUBST([YODAVERSION], [`$YODABIN --version`])
      AC_SUBST([YODABINPATH], ["$YODAPREFIX/bin"])
      AC_SUBST([YODALIBPATH], [`$YODABIN --libdir`])
      AC_SUBST([YODAINCLUDEPATH],[`$YODABIN --includedir`])
      AC_SUBST([YODACPPFLAGS], [`$YODABIN --cppflags`])
      AC_SUBST([YODALDFLAGS], [`$YODABIN --libs`])
      AC_MSG_RESULT([Using $TMPNAM version $YODAVERSION])
      AC_MSG_RESULT([Using $YODAPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $YODABINPATH for the binaries])
      AC_MSG_RESULT([Using $YODALIBPATH for the library])
      AC_MSG_RESULT([Using $YODAINCLUDEPATH for the includes])
      AC_MSG_RESULT([Using $YODACPPFLAGS for the CPP flags])
      AC_MSG_RESULT([Using $YODALDFLAGS for the libraries])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([No $TMPNAM support requested. Continue without $TMPNAM.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_YODA, test "x$YODAPREFIX" != "x")

# ROOT
AC_SUBST([TMPNAM],["ROOT"])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether $TMPNAM support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([root],
        [AS_HELP_STRING([--with-root],
        [Specify location of $TMPNAM installation. This takes precedence over
        a ROOTSYS environment variable, the prefix path or other methods.
        Please note: $TMPNAM can optionally be used in $PACKAGE_NAME!])],
        [if test "x$withval" != "xyes" -a "x$withval" != "xno"; then
            AC_SUBST([ROOTSYS],[$withval])
            AC_MSG_NOTICE([Using given path $ROOTSYS as location of $TMPNAM installation.])
        fi],
        [with_root="no"])

if test "x$with_root" != "xno"; then
   if test "x$ROOTSYS" = "x"; then
      AC_MSG_NOTICE([No $TMPNAM installation path set. Try to find $TMPNAM in prefix path ...])
      AC_PATH_PROG([ROOTBIN],[root-config],[],[$prefix/bin])
      if test "x$ROOTBIN" = "x"; then
         AC_MSG_NOTICE([No $TMPNAM installation found in prefix path. Try to find $TMPNAM in system path ...])
         AC_PATH_PROG([ROOTBIN],[root-config])
         if test "x$ROOTBIN" = "x"; then
            AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
         fi
      fi
   else
      AC_PATH_PROG([ROOTBIN],[root-config],[],[$ROOTSYS/bin])
      if test "x$ROOTBIN" = "x"; then
         AC_MSG_ERROR([root-config not found in given path to $TMPNAM! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$ROOTBIN" != "x"; then
      AC_DEFINE([WITH_ROOT],[1],[Enable conditional compilation with ROOT])
      AC_MSG_NOTICE([$TMPNAM settings:])
      AC_SUBST([ROOTPREFIX],  [`$ROOTBIN --prefix`])
      AC_SUBST([ROOTVERSION], [`$ROOTBIN --version`])
      AC_SUBST([ROOTBINPATH], [`$ROOTBIN --bindir`])
      AC_SUBST([ROOTLIBPATH], [`$ROOTBIN --libdir`])
      AC_SUBST([ROOTINCLUDEPATH],[`$ROOTBIN --incdir`])
      AC_SUBST([ROOTCPPFLAGS], [`$ROOTBIN --cflags`])
      AC_SUBST([ROOTCXXFLAGS], [`$ROOTBIN --libs`])
      AC_SUBST([ROOTLDFLAGS], [`$ROOTBIN --ldflags`])
      AC_MSG_RESULT([Using $TMPNAM version $ROOTVERSION])
      AC_MSG_RESULT([Using $ROOTPREFIX for the $TMPNAM installation path])
      AC_MSG_RESULT([Using $ROOTBINPATH for the binaries])
      AC_MSG_RESULT([Using $ROOTLIBPATH for the library])
      AC_MSG_RESULT([Using $ROOTINCLUDEPATH for the includes])
      AC_MSG_RESULT([Using $ROOTCPPFLAGS for the precompiler flags])
      AC_MSG_RESULT([Using $ROOTCXXFLAGS for the compiler flags])
      AC_MSG_RESULT([Using $ROOTLDFLAGS for the linker flags])
   else
      AC_MSG_ERROR([No $TMPNAM library found. Configuration aborted!])
   fi
else
   AC_MSG_NOTICE([No $TMPNAM support requested. Continue without $TMPNAM.])
   AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_ROOT, test "x$ROOTPREFIX" != "x")

# Create code documentation
# Check on DOXYGEN
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Check for presence of doxygen ...])
AC_MSG_NOTICE([******************************************])
AC_CHECK_PROGS(DOXYGEN, doxygen, doxygen)
if test "x$DOXYGEN" = "x"; then
   AC_MSG_WARN([Did not find doxygen. Continue without doxygen documentation.])
else
   DOXYVMIN="1.5"
   DOXYVERS=`doxygen --version`
   DOXYVORD=`echo -e "$DOXYVERS\n$DOXYVMIN" | sort -t"." -g | head -n1`
   if test "x$DOXYVORD" = "x$DOXYVMIN"; then
      AC_MSG_NOTICE([Creating code documentation using doxygen version $DOXYVERS.])
      AC_MSG_NOTICE([******************************************])
      DX_DOT_FEATURE(ON)
      AC_CHECK_PROGS(DOT, dot, dot)
      if test "x$DOT" = "x"; then
         AC_MSG_WARN([Did not find dot. Continue doxygen documentation without dot graphics.])
         DX_DOT_FEATURE(OFF)
      fi
      DX_HTML_FEATURE(ON)
      DX_CHM_FEATURE(OFF)
      DX_CHI_FEATURE(OFF)
      DX_MAN_FEATURE(OFF)
      DX_RTF_FEATURE(OFF)
      DX_XML_FEATURE(OFF)
      DX_PDF_FEATURE(ON)
      DX_PS_FEATURE(OFF)
   else
      AC_MSG_WARN([Doxygen version $DOXYVERS too old. Continue without doxygen documentation.])
      DOXYGEN=""
   fi
fi
AM_CONDITIONAL(HAVE_DOXYGEN, test "x$DOXYGEN" != "x")

# Create documentation
# Check first on PDFLATEX ...
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Check for presence of pdflatex ...])
AC_MSG_NOTICE([******************************************])
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test "x$PDFLATEX" = "x"; then
  AC_MSG_WARN([Unable to create PDF documentation with pdflatex, check for latex & dvipdfm.])
# Then LATEX & DVIPDFM
  AC_MSG_NOTICE([******************************************])
  AC_MSG_NOTICE([Check for presence of latex & dvipdfm ...])
  AC_MSG_NOTICE([******************************************])
  AC_CHECK_PROG(LATEX, latex, latex)
  if test "x$LATEX" = "x"; then
    AC_MSG_WARN([Unable to create PDF documentation with latex & dvipdfm, continue without.])
  else
    AC_CHECK_PROG(DVIPDFM, dvipdfm, dvipdfm)
    if test "x$DVIPDFM" = "x"; then
      AC_MSG_WARN([Unable to create PDF documentation with latex & dvipdfm, continue without.])
    else
      AC_MSG_NOTICE([Creating PDF documentation using latex & dvipdfm.])
      AC_MSG_NOTICE([******************************************])
    fi
  fi
else
  AC_MSG_NOTICE([Creating PDF documentation using pdflatex.])
  AC_MSG_NOTICE([******************************************])
fi
AM_CONDITIONAL(HAVE_PDFLATEX, test "x$PDFLATEX" != "x")
AM_CONDITIONAL(HAVE_LATEXDVIPDFM, test "x$LATEX" != "x" && test "x$DVIPDFM" != "x")

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
if test "x$enable_pyext" = "xyes"; then
   AX_PYTHON_DEVEL([>='2.3.1'])
   AM_PATH_PYTHON
fi
# Check on existence of SWIG to interface to Python
if test "x$enable_pyext" = "xyes" -a ! -e pyext/FastNLOReader_wrap.cpp; then
   AC_PATH_PROG(SWIG, swig)
   if test "x$SWIG" =  "x"; then
      AC_MSG_ERROR([swig not found, interfacing to Python not possible! Please install swig first or disable Python extension. Configuration aborted!])
   else
      AC_MSG_NOTICE([swig found, interface to Python possible.])
   fi
fi
AC_MSG_NOTICE([******************************************])
if test "x$enable_pyext" = "xyes"; then
   AC_MSG_NOTICE([All Python checks successful. Requested Python extension will be built.])
else
   AC_MSG_NOTICE([No Python extension requested. Continue without pyext.])
fi
AC_MSG_NOTICE([******************************************])
AM_CONDITIONAL(ENABLE_PYEXT, [test "x$enable_pyext" = "xyes"])

# Add subversion revision number, recipe from Autotools by J. Calcote
# Modified on 18.12.13 to fix changed behaviour of svnversion to give "Unversioned directory" instead of "exported"
SVNREV=`( svnversion $srcdir | sed 's/:.*//' ) 2>/dev/null`
if { ! ( svnversion ) >/dev/null 2>&1 || test "x$SVNREV" = "xexported" || test "x$SVNREV" = "xUnversioned directory"; } ;
  then SVNREV=`cat $srcdir/SVNREV`
  else echo $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Initialize Doxygen doc
if test "x$DOXYGEN" != "x"; then
   DX_INIT_DOXYGEN($PROJECT, doxygen.cfg, doc/doxygen)
# or set missing conditionals to false
else
   DX_COND_doc_TRUE=
   DX_COND_doc_FALSE='#'
   DX_COND_dot_TRUE=
   DX_COND_dot_FALSE='#'
   DX_COND_man_TRUE=
   DX_COND_man_FALSE='#'
   DX_COND_rtf_TRUE=
   DX_COND_rtf_FALSE='#'
   DX_COND_html_TRUE=
   DX_COND_html_FALSE='#'
   DX_COND_chm_TRUE=
   DX_COND_chm_FALSE='#'
   DX_COND_chi_TRUE=
   DX_COND_chi_FALSE='#'
   DX_COND_xml_TRUE=
   DX_COND_xml_FALSE='#'
   DX_COND_latex_TRUE=
   DX_COND_latex_FALSE='#'
   DX_COND_ps_TRUE=
   DX_COND_ps_FALSE='#'
   DX_COND_pdf_TRUE=
   DX_COND_pdf_FALSE='#'
fi

# Define project substitution variables (to be replaced e.g. in config files below)
AC_SUBST([AUTHORS],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, G. Sieber, M. Wobisch"])
AC_SUBST([AUTHORSv14],["T. Kluge, K. Rabbertz, M. Wobisch"])
AC_SUBST([AUTHORSv2],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"])
AC_SUBST([QUOTEv14],["hep-ph/0609285"])
AC_SUBST([QUOTEv2],["arXiv:1109.1310"])
AC_SUBST([SUBPROJECT],["toolkit"])
AC_SUBST([WEBPAGE],["http://projects.hepforge.org/fastnlo"])
AC_SUBST([YEARS],["2005-`date +%Y`"])

# Output
AC_CONFIG_FILES([Makefile
                 data/Makefile
                 data/check/Makefile
                 data/Modify/Makefile
                 doc/Makefile
                 fastnlotoolkit/Makefile
                 fastnlotoolkit/include/Makefile
                 fastnlotoolkit/include/fastnlotk/Makefile
                 fastnlotoolkit/include/fastnlotk/fastNLOConstants.h
                 pyext/Makefile
                 pyext/setup.py
                 scripts/Makefile
                 src/Makefile
                 src/fnlo-tk-config
                ])

AC_OUTPUT

#                 doc/paper/Makefile
