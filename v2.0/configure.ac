#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Require minimal autoconf version of 2.63
AC_PREREQ([2.59])
# Define project fastNLO
AC_INIT([fastNLO], [2.0.0], [klaus.rabbertz@cern.ch,wobisch@fnal.gov])
# Require minimal automake version of 1.11, switch on warnings
AM_INIT_AUTOMAKE([1.9 -Wall])

#
AC_CONFIG_SRCDIR([common/fnloTable.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set default installation directory
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_F77
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
AC_MSG_RESULT([Using $LHPDFPREFIX for the LHPDF installation path])
AC_MSG_RESULT([Using $LHPDFLIBPATH for the libraries])
AC_MSG_RESULT([Using $LHPDFINCLUDEPATH for the includes])

# LHAPDF
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([First try to find LHAPDF via lhapdf-config executable in system PATH ...])
AC_PATH_PROG([LHAPDFBIN],[lhapdf-config])
if test "x$LHAPDFBIN" = "x"; then
  AC_MSG_WARN([No path to LHAPDF binaries found, try to use LHAPDF environment variable ...])
  if test "x$LHAPDF" = "x"; then
    AC_MSG_WARN([LHAPDF environment variable is not set, please install LHAPDF. fastNLO user code may not work without!])
  else
    AC_MSG_RESULT([Using environment variable LHAPDF=$LHAPDF to locate LHAPDF installation ...])
    AC_MSG_RESULT([Searching for executable lhapdf-config in $LHAPDF/bin ...])
    AC_PATH_PROG([LHAPDFBIN],[lhapdf-config],[],[$LHAPDF/bin])
    if test "x$LHAPDFBIN" = "x"; then
      AC_MSG_WARN([LHAPDF binary not found in $LHAPDF/bin, please install LHAPDF. fastNLO user code may not work without!])
    else
      AC_SUBST([LHAPDFPREFIX], [`$LHAPDFBIN --prefix`])
    fi
  fi
else
  AC_SUBST([LHAPDFPREFIX], `$LHAPDFBIN --prefix`)
fi
AC_MSG_NOTICE([******************************************])
if test "x$LHAPDFPREFIX" = "x"; then
  AC_MSG_WARN([No installation of LHAPDF found. fastNLO user code may not work without!])
else
  AC_SUBST([LHAPDFSETPATH], [`$LHAPDFBIN --pdfsets-path`])
  AC_SUBST([LHAPDFBINPATH], ["$LHAPDFPREFIX/bin"])
  AC_SUBST([LHAPDFLIBPATH], [`$LHAPDFBIN --libdir`])
  AC_SUBST([LHAPDFINCLUDEPATH],[`$LHAPDFBIN --incdir`])
  AC_MSG_RESULT([Using $LHAPDFPREFIX for the LHAPDF installation path])
  AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])
  AC_MSG_RESULT([Using $LHAPDFBINPATH for the binaries])
  AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
  AC_MSG_RESULT([Using $LHAPDFINCLUDEPATH for the includes])
fi

# FASTJET
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for fastjet installation ...])
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([First try to find fastjet via fastjet-config executable in system PATH ...])
AC_PATH_PROG([FASTJETBIN],[fastjet-config])
if test "x$FASTJETBIN" = "x"; then
  AC_MSG_WARN([No path to fastjet binary found, try to use FASTJET environment variable ...])
  if test "x$FASTJET" = "x"; then
    AC_MSG_WARN([FASTJET environment variable is not set, please install fastjet. Parts of fastNLO may not work without!])
  else
    AC_MSG_RESULT([Using environment variable FASTJET=$FASTJET to locate FASTJET installation ...])
    AC_MSG_RESULT([Searching for executable fastjet-config in $FASTJET/bin ...])
    AC_PATH_PROG([FASTJETBIN],[fastjet-config],[],[$FASTJET/bin])
    if test "x$FASTJETBIN" = "x"; then
      AC_MSG_WARN([FASTJET binaries not found in $FASTJET/bin, please install fastjet. Parts of fastNLO may not work without!])
    else
      AC_SUBST([FASTJETPREFIX], [`$FASTJETBIN --prefix`])
    fi
  fi
else
  AC_SUBST([FASTJETPREFIX], [`$FASTJETBIN --prefix`])
fi
AC_MSG_NOTICE([******************************************])
if test "x$FASTJETPREFIX" = "x"; then
  AC_MSG_WARN([No installation of fastjet found. Parts of fastNLO may not work!])
else
  AC_SUBST([FASTJETBINPATH], ["$FASTJETPREFIX/bin"])
  AC_SUBST([FASTJETLIBPATH], ["$FASTJETPREFIX/lib"])
  AC_SUBST([FASTJETINCLUDEPATH],["$FASTJETPREFIX/include"])
  AC_SUBST([FASTJETLIBS],[`$FASTJETBINPATH/fastjet-config --libs`])
  AC_MSG_RESULT([Using $FASTJETPREFIX for the FASTJET installation path])
  AC_MSG_RESULT([Using $FASTJETBINPATH for the binaries])
  AC_MSG_RESULT([Using $FASTJETLIBPATH for the library])
  AC_MSG_RESULT([Using $FASTJETINCLUDEPATH for the includes])
  AC_MSG_RESULT([Using $FASTJETLIBS for the linker flags])
fi

# CERNLIB
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for CERN libraries ...])
AC_MSG_NOTICE([******************************************])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_WARN([CERN_ROOT environment variable is not set. Continue without CERNlib.])
else
  AC_MSG_RESULT([Using environment variable CERN_ROOT=$CERN_ROOT to locate CERNlib installation])
  AC_SUBST([CERNLIBPATH], ["$CERN_ROOT/lib"])
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using $CERN_ROOT for CERNlib installation path])
  AC_MSG_RESULT([Using $CERNLIBPATH for the libraries])
  AC_MSG_RESULT([Using -L$CERNLIBPATH -lmathlib -lkernlib -lpacklib for the linker flag])
fi

# ROOT
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Looking for ROOT installation ...])
AC_MSG_NOTICE([******************************************])
#AC_PATH_PROG([ROOTBIN],[root-config])
#if test "x$ROOTBIN" = "x"; then
#  AC_MSG_WARN([No path to ROOT binaries found, try to use ROOTSYS environment variable ...])
#  if test "x$ROOTSYS" = "x"; then
#    AC_MSG_WARN([ROOTSYS environment variable is not set. Continue without ROOT.])
#  else
#    AC_MSG_RESULT([Using environment variable ROOTSYS=$ROOTSYS to locate ROOT installation ...])
#    AC_MSG_RESULT([Searching for executable root-config in $ROOTSYS/bin ...])
#    AC_PATH_PROG([ROOTBIN],[root-config],[],[$ROOTSYS/bin])
#    if test "x$ROOTBIN" = "x"; then
#      AC_MSG_WARN([ROOT binaries not found in $ROOTSYS/bin. Continue without ROOT.])
#    else
#      AC_SUBST([ROOTPREFIX], [`$ROOTBIN --prefix`])
#    fi
#  fi
#else
#  AC_SUBST([ROOTPREFIX], [`$ROOTBIN --prefix`])
#fi
#AC_MSG_NOTICE([******************************************])
#if test "x$ROOTPREFIX" = "x"; then
#  AC_MSG_WARN([No installation of ROOT found, will work without ROOT dictionaries.])
#else
#  AC_SUBST([ROOTSYS], [$ROOTPREFIX])
#  AC_SUBST([ROOTBINPATH], [`$ROOTBIN --bindir`])
#  AC_SUBST([ROOTLIBPATH], [`$ROOTBIN --libdir`])
#  AC_SUBST([ROOTINCLUDEPATH], [`$ROOTBIN --incdir`])
#fi
AC_MSG_WARN([ROOT will be ignored for the moment!])
#AC_CHECK_PROGS([have_rootcint],[rootcint],[],["$ROOTSYS/bin"])
#if test "x$have_rootcint" = "x"; then
#  AC_MSG_WARN([ROOTSYS environment variable set and/or root-config found, but no rootcint. Not a big problem, but library will be build without ROOT dictionary.])
#else
# No use of ROOT for the moment
#  AC_SUBST([ROOTCINTPATH], ["$ROOTSYS/bin/rootcint"])  
  AC_SUBST([ROOTCINTPATH], [""])  
#  AC_MSG_RESULT([Using $ROOTSYS for the ROOT installation path])
#  AC_MSG_RESULT([Using $ROOTBINPATH for the binaries])
#  AC_MSG_RESULT([Using $ROOTINCLUDEPATH for the includes])
#  AC_MSG_RESULT([Using $ROOTLIBPATH for the libraries])
#  AC_MSG_RESULT([Using $ROOTCINTPATH for rootcint path])
#fi
AM_CONDITIONAL(HAVEROOT, test "x$ROOTCINTPATH" != "x")

# Add subversion revision number, recipe from Autotools by J. Calcote 
SVNREV=`svnversion $srcdir | sed 's/:.*//'`
which svnversion > /dev/null; \
if [ $? -ne 0 ] || [ "x$SVNREV" = "xexported" ]
  then SVNREV=`cat $srcdir/SVNREV`
  else echo -n $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

# Output
#                 user/Makefile
#                 author/hadron/Makefile
#	         tableconvert/Makefile
AC_CONFIG_FILES([Makefile
                 author/ep/Makefile
                 common/Makefile
                 merger/Makefile
	        ])

AC_OUTPUT
