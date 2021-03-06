#                                               -*- Autoconf -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.02.2009
#
# Create distribution and installation package for fastNLO v2
#
# -----------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.59
AC_PREREQ([2.59])
# Define project fastNLO
AC_INIT([fastNLO], [2.1.0], [daniel.britzger@desy.de,klaus.rabbertz@cern.ch,g.sieber@cern.ch,stober@cern.ch,wobisch@fnal.gov])
# Properly include subprojects
AC_CONFIG_AUX_DIR([.])
# Decouple the 'reader' and 'processor/converter' subprojects ?
AC_CONFIG_SUBDIRS([reader processor/converter])
# Require minimal automake version of 1.9, switch on warnings
AM_INIT_AUTOMAKE([1.9 -Wall])
# Require minimal libtool version of 1.5
#LT_PREREQ([1.5])
#LT_INIT([dlopen])
#
AC_CONFIG_SRCDIR([common/fnloTable.h])
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

# Announce checks particular for fastNLO
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO ...])
AC_MSG_NOTICE([******************************************])

# NLOJet++
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for NLOJet++ installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([nlojet],
        [AS_HELP_STRING([--with-nlojet],
        [Specify location of NLOJet++ installation. This takes precedence over
        an NLOJET environment variable, the prefix path or other methods.
        Please note: fastNLO CANNOT be used without NLOJet++!])],
        [if test "x$withval" = "xno"; then
        AC_MSG_ERROR([fastNLO CANNOT be used without NLOJet++! Configuration aborted!])
        fi
        if test "x$withval" != "xyes"; then
        AC_SUBST([NLOJET],[$withval])
        AC_MSG_NOTICE([Using given path $NLOJET as location of NLOJet++ installation.])
        fi],
        [])

if test "x$NLOJET" = "x"; then
   AC_MSG_NOTICE([No NLOJet++ installation path set. Try to find NLOJet++ in prefix path ...])
   AC_PATH_PROG([NLOJETBIN],[nlojet++],[],[$prefix/bin])
   if test "x$NLOJETBIN" = "x"; then
      AC_MSG_NOTICE([No NLOJet++ installation found in prefix path. Try to find NLOJet++ in system path ...])
      AC_PATH_PROG([NLOJETBIN],[nlojet++])
      if test "x$NLOJETBIN" = "x"; then
         AC_MSG_ERROR([No NLOJet++ installation found. Please install NLOJet++ first. Configuration aborted!])
      fi
   fi
   AC_SUBST([NLOJETBINPATH], [`dirname $NLOJETBIN`])
   AC_SUBST([NLOJET], [`dirname $NLOJETBINPATH`])
fi

AC_SUBST([NLOJETPREFIX], "$NLOJET")
AC_SUBST([NLOJETBINPATH], "$NLOJETPREFIX/bin")
AC_SUBST([NLOJETLIBPATH], "$NLOJETPREFIX/lib")
AC_SUBST([NLOJETINCLUDEPATH], "$NLOJETPREFIX/include/nlo++")
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([NLOJet++ settings:])
AC_MSG_RESULT([Using $NLOJETPREFIX for the NLOJET++ installation path])
AC_MSG_RESULT([Using $NLOJETBINPATH for the binaries])
AC_MSG_RESULT([Using $NLOJETLIBPATH for the libraries])
AC_MSG_RESULT([Using $NLOJETINCLUDEPATH for the includes])

# CTEQ-PDF (comes with NLOJET++ 4.1.3)
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for CTEQ-PDF (comes with NLOJET++ 4.1.3) installation ...])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Attention! These are the PDFs coming with NLOJET++ version 4.1.x!])
AC_MSG_NOTICE([Do not confuse with LHAPDF!])

AC_ARG_WITH([cteq-pdf],
        [AS_HELP_STRING([--with-cteq-pdf],
        [Specify location of CTEQ-PDF installation. This takes precedence over
        a ZNPDF environment variable, the prefix path or other methods.
        Please note: fastNLO currently CANNOT be used without CTEQ-PDF
        which is required by NLOJet++!])],
        [if test "x$withval" = "xno"; then
        AC_MSG_ERROR([fastNLO CANNOT be used without CTEQ-PDF! Configuration aborted!])
        fi
        if test "x$withval" != "xyes"; then
        AC_SUBST([ZNPDF],[$withval])
        AC_MSG_NOTICE([Using given path $ZNPDF as location of CTEQ-PDF installation.])
        fi],
        [])

if test "x$ZNPDF" = "x"; then
   AC_MSG_NOTICE([No ZNPDF installation path set. Try to find CTEQ-PDF in prefix path ...])
   AC_CHECK_FILE([$prefix/lib/libcteqpdf.la],
   [precteqpdf="$ac_cv_file_$prefix/lib/libcteqpdf.la"],
   [precteqpdf=""])
   if test "x$precteqpdf" = "x"; then
      AC_MSG_NOTICE([No CTEQ-PDF library found in prefix path. Try to link with system's CTEQ-PDF ...])
      AC_CHECK_LIB([cteqpdf],[libcteqpdf_la-cteqpdf],
      [syscteqpdf=$ac_cv_lib_libcteqpdf_la-cteqpdf],
      [syscteqpdf=""])
      if test "x$syscteqpdf" = "x"; then
         AC_MSG_ERROR([No CTEQ-PDF library found. Please install CTEQ-PDF first. Configuration aborted!])
      else
         AC_SUBST([ZNPDFLIBPATH], [`dirname $syscteqpdf`])
      fi
   else
      AC_SUBST([ZNPDFLIBPATH], [`dirname $precteqpdf`])
   fi
   AC_SUBST([ZNPDF], [`dirname $ZNPDFLIBPATH`])
fi

AC_SUBST([ZNPDFPREFIX], "$ZNPDF")
AC_SUBST([ZNPDFLIBPATH], "$ZNPDFPREFIX/lib")
AC_SUBST([ZNPDFINCLUDEPATH], "$ZNPDFPREFIX/include/cteq")
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([CTEQ-PDF settings:])
AC_MSG_RESULT([Using $ZNPDFPREFIX for the CTEQ-PDF installation path])
AC_MSG_RESULT([Using $ZNPDFLIBPATH for the libraries])
AC_MSG_RESULT([Using $ZNPDFINCLUDEPATH for the includes])

# LHPDF (comes with NLOJET++ 4.0.1)
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHPDF (comes with NLOJET++ 4.0.1) installation ...])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Attention! These are the PDFs coming with NLOJET++ version 4.0.x!])
AC_MSG_NOTICE([Do not confuse with LHAPDF!])

AC_ARG_WITH([lhpdf],
        [AS_HELP_STRING([--with-lhpdf],
        [Specify location of LHPDF installation. This takes precedence over
        an LHPDF environment variable, the prefix path or other methods.
        Please note: fastNLO currently CANNOT be used without LHPDF
        which is required by NLOJet++!])],
        [if test "x$withval" = "xno"; then
        AC_MSG_ERROR([fastNLO CANNOT be used without LHPDF! Configuration aborted!])
        fi
        if test "x$withval" != "xyes"; then
        AC_SUBST([LHPDF],[$withval])
        AC_MSG_NOTICE([Using given path $LHPDF as location of LHPDF installation.])
        fi],
        [])

if test "x$LHPDF" = "x"; then
   AC_MSG_NOTICE([No LHPDF installation path set. Try to find LHPDF in prefix path ...])
   AC_CHECK_FILE([$prefix/lib/liblhpdf.la],
   [prelhpdf="$ac_cv_file_$prefix/lib/liblhpdf.la"],
   [prelhpdf=""])
   if test "x$prelhpdf" = "x"; then
      AC_MSG_NOTICE([No LHPDF library found in prefix path. Try to link with system's LHPDF ...])
      AC_CHECK_LIB([lhpdf],[lhpdfset],
      [syslhpdf=$ac_cv_lib_lhpdfset],
      [syslhpdf=""])
      if test "x$syslhpdf" = "x"; then
         AC_MSG_ERROR([No LHPDF library found. Please install LHPDF first. Configuration aborted!])
      else
         AC_SUBST([LHPDFLIBPATH], [`dirname $syslhpdf`])
      fi
   else
      AC_SUBST([LHPDFLIBPATH], [`dirname $prelhpdf`])
   fi
   AC_SUBST([LHPDF], [`dirname $LHPDFLIBPATH`])
fi

AC_SUBST([LHPDFPREFIX], "$LHPDF")
AC_SUBST([LHPDFLIBPATH], "$LHPDFPREFIX/lib")
AC_SUBST([LHPDFINCLUDEPATH], "$LHPDFPREFIX/include/lhpdf")
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([LHPDF settings:])
AC_MSG_RESULT([Using $LHPDFPREFIX for the LHPDF installation path])
AC_MSG_RESULT([Using $LHPDFLIBPATH for the libraries])
AC_MSG_RESULT([Using $LHPDFINCLUDEPATH for the includes])

# LHAPDF
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([lhapdf],
        [AS_HELP_STRING([--with-lhapdf],
        [Specify location of LHAPDF installation. This takes precedence over
        an LHAPDF environment variable, the prefix path or other methods.
        Please note: LHAPDF is mandatory in order to use fastNLO!])],
        [if test "x$withval" = "xno"; then
            AC_MSG_ERROR([fastNLO CANNOT be used without LHAPDF! Configuration aborted!])
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
            AC_MSG_NOTICE([Building fastNLO without LHAPDF!])
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
      AC_MSG_NOTICE([No LHAPDF library found. Building fastNLO without LHAPDF!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([LHAPDF support switched off by user. Building fastNLO without LHAPDF!])
fi
AM_CONDITIONAL(HAVE_LHAPDF, test "x$LHAPDFPREFIX" != "x")
AM_CONDITIONAL(HAVE_LHAPDFCONFIG, test "x$LHAPDFBIN" != "x")

# FASTJET
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for fastjet installation ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([fastjet],
        [AS_HELP_STRING([--with-fastjet],
         [Specify location of FASTJET installation. This takes precedence over
          a FASTJET environment variable, the prefix path or other methods.])],
        [if test "x$withval" = "xno"; then
            AC_MSG_NOTICE([Building fastNLO without FASTJET!])
         else
            if test "x$withval" != "xyes"; then
               AC_SUBST([FASTJET],[$withval])
               AC_MSG_NOTICE([Using given path $FASTJET as location of FASTJET installation.])
            fi
         fi],
        [])

if test "x$withval" != "xno"; then
   if test "x$FASTJET" = "x"; then
      AC_MSG_NOTICE([No FASTJET installation path set. Try to find FASTJET in prefix path ...])
      AC_PATH_PROG([FASTJETBIN],[fastjet-config],[],[$prefix/bin])
      if test "x$FASTJETBIN" = "x"; then
         AC_MSG_NOTICE([No FASTJET installation found in prefix path. Try to find FASTJET in system path ...])
         AC_PATH_PROG([FASTJETBIN],[fastjet-config])
         if test "x$FASTJETBIN" = "x"; then
            AC_MSG_NOTICE([Building fastNLO without FASTJET!])
         fi
      fi
   else
      AC_PATH_PROG([FASTJETBIN],[fastjet-config],[],[$FASTJET/bin])
      if test "x$FASTJETBIN" = "x"; then
         AC_MSG_ERROR([fastjet-config not found in given path to FASTJET! Configuration aborted!])
      fi
   fi
   AC_MSG_NOTICE([******************************************])
   if test "x$FASTJETBIN" != "x"; then
      AC_MSG_NOTICE([FASTJET settings:])
      AC_SUBST([FASTJETPREFIX],  [`$FASTJETBIN --prefix`])
      AC_SUBST([FASTJETBINPATH], ["$FASTJETPREFIX/bin"])
      AC_SUBST([FASTJETLIBPATH], ["$FASTJETPREFIX/lib"])
      AC_SUBST([FASTJETINCLUDEPATH],["$FASTJETPREFIX/include"])
      AC_SUBST([FASTJETLIBS],[`$FASTJETBINPATH/fastjet-config --libs`])
      AC_MSG_RESULT([Using $FASTJETPREFIX for the FASTJET installation path])
      AC_MSG_RESULT([Using $FASTJETBINPATH for the binaries])
      AC_MSG_RESULT([Using $FASTJETLIBPATH for the library])
      AC_MSG_RESULT([Using $FASTJETINCLUDEPATH for the includes])
      AC_MSG_RESULT([Using $FASTJETLIBS for the linker flags])
   else
      AC_MSG_NOTICE([No FASTJET library found. Building fastNLO without FASTJET!])
   fi
else
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([FASTJET support switched off by user. Building fastNLO without FASTJET!])
fi
AM_CONDITIONAL(HAVE_FASTJET, test "x$FASTJETPREFIX" != "x")
AM_CONDITIONAL(HAVE_FASTJETCONFIG, test "x$FASTJETBIN" != "x")
#AC_MSG_ERROR([STOPP])

# CERNLIB
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for CERN libraries ...])
AC_MSG_NOTICE([******************************************])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_NOTICE([CERN_ROOT environment variable is not set. Continue without CERNlib.])
else
  AC_MSG_RESULT([Using environment variable CERN_ROOT=$CERN_ROOT to locate CERNlib installation])
  AC_SUBST([CERNLIBPATH], ["$CERN_ROOT/lib"])
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib"])
  AC_MSG_NOTICE([******************************************])
  AC_MSG_NOTICE([CERNlib settings:])
  AC_MSG_RESULT([Using $CERN_ROOT for CERNlib installation path])
  AC_MSG_RESULT([Using $CERNLIBPATH for the libraries])
  AC_MSG_RESULT([Using -L$CERNLIBPATH -lmathlib -lkernlib -lpacklib for the linker flag])
fi
AM_CONDITIONAL(HAVE_CERNLIB, test "x$CERNLIBPATH" != "x")

# ROOT
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Checking whether ROOT support is requested ...])
AC_MSG_NOTICE([******************************************])

AC_ARG_WITH([root],
        [AS_HELP_STRING([--with-root],
         [Specify location of ROOT installation. By default fastNLO is built
          without ROOT support.])],
        [if test "x$withval" = "xno"; then
            AC_MSG_NOTICE([Building fastNLO without ROOT support!])
         else
            if test "x$withval" != "xyes"; then
               AC_SUBST([ROOTSYS],[$withval])
               AC_MSG_NOTICE([Using given path $ROOTSYS as location of ROOT installation.])
               AC_PATH_PROG([ROOTBIN],[root-config],[],[$ROOTSYS/bin])
               if test "x$ROOTBIN" = "x"; then
                  AC_MSG_ERROR([ROOT binaries not found in $ROOTSYS/bin! Configuration aborted!])
               fi
               AC_PATH_PROG([ROOTCINT],[rootcint],[],["$ROOTSYS/bin"])
               if test "x$ROOTCINT" = "x"; then
                  AC_MSG_ERROR([ROOTSYS variable set and/or root-config found, but no rootcint!
                                ROOT support not possible, configuration aborted!])
               fi
            else
               AC_MSG_ERROR(["yes" is not a valid ROOTSYS path! Configuration aborted!])
            fi
         fi],
        [])

if test "x$ROOTBIN" != "x"; then
   AC_SUBST([ROOTPREFIX],  [`$ROOTBIN --prefix`])
   AC_SUBST([ROOTBINPATH], [`$ROOTBIN --bindir`])
   AC_SUBST([ROOTLIBPATH], [`$ROOTBIN --libdir`])
   AC_SUBST([ROOTINCLUDEPATH],[`$ROOTBIN --incdir`])
   AC_SUBST([ROOTLIBS],[`$ROOTBIN --libs`])
   AC_MSG_NOTICE([******************************************])
   AC_MSG_NOTICE([ROOT settings:])
   AC_MSG_RESULT([Using $ROOTPREFIX for the ROOT installation path])
   AC_MSG_RESULT([Using $ROOTBINPATH for the binaries])
   AC_MSG_RESULT([Using $ROOTLIBPATH for the libraries])
   AC_MSG_RESULT([Using $ROOTINCLUDEPATH for the includes])
   AC_MSG_RESULT([Using $ROOTLIBS for the linker flags])
else
   AC_MSG_NOTICE([No ROOT support requested. Continue without ROOT.])
fi
AM_CONDITIONAL(HAVE_ROOT, test "x$ROOTCINT" != "x")
AC_MSG_NOTICE([******************************************])

# Create PDF documentation?
# First PDFLATEX ...
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

# Add subversion revision number, recipe from Autotools by J. Calcote
SVNREV=`svnversion $srcdir | sed 's/:.*//'`
which svnversion > /dev/null; \
if [[ $? -ne 0 ] || [ "x$SVNREV" = "xexported" ]];
  then SVNREV=`cat $srcdir/SVNREV`
  else echo -n $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Define project substitution variables (to be replaced e.g. in config files below)
AC_SUBST([AUTHORS],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, G. Sieber, M. Wobisch"])
AC_SUBST([AUTHORSv14],["T. Kluge, K. Rabbertz, M. Wobisch"])
AC_SUBST([AUTHORSv2],["D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"])
AC_SUBST([QUOTEv14],["hep-ph/0609285"])
AC_SUBST([QUOTEv2],["arXiv:1109.1310"])
AC_SUBST([WEBPAGE],["http://projects.hepforge.org/fastnlo"])

# Output (not listed as subproject in AC_CONFIG_SUBDIRS above)
AC_CONFIG_FILES([Makefile
        common/Makefile
        creator/ep/Makefile
        creator/hadron/Makefile
        creator/Makefile
        doc/arXiv_1109.1310/Makefile
        doc/Makefile
        processor/Makefile
        processor/modifier/Makefile
])


AC_OUTPUT
