# Process this file with autoconf to produce a configure script.
AC_INIT([fastnlo], [2.0], [thomas.kluge@cern.ch])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_CXX
AC_PROG_F77
AC_PROG_CPP
AM_PROG_LIBTOOL

# Checks for sources.
AC_CONFIG_SRCDIR([common/fnloconstants.h])

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS

# Announce checks particluar for fastNLO
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO ...])
AC_MSG_NOTICE([******************************************])

#NLOJET++
AC_MSG_NOTICE([Looking for NLOJET++ installation ...])
if test "x$NLOJET" = "x"; then
  AC_MSG_WARN([NLOJET environment variable is not set, trying to find create-nlojet-user in PATH ...])
  AC_PATH_PROG([NLOJETBINPATH],[create-nlojet-user],[:])
  if test "$NLOJETBINPATH" = :; then
    AC_MSG_ERROR([No path to NLOJET++ binaries found, configuration aborted!])
  else
    AC_SUBST([NLOJETPREFIX], "$NLOJETBINPATH/..")
    AC_SUBST([NLOJETINCLUDEPATH], "$NLOJETPREFIX/include/nlo++")
  fi
else
  AC_MSG_RESULT([Using environment variable NLOJET=$NLOJET to locate NLOJET++ installation])
  AC_SUBST([NLOJETPREFIX], "$NLOJET")
  AC_SUBST([NLOJETBINPATH], "$NLOJET/bin")
  AC_SUBST([NLOJETINCLUDEPATH], "$NLOJETPREFIX/include/nlo++")
fi
AC_MSG_RESULT([Using $NLOJETPREFIX for the NLOJET++ installation path])
AC_MSG_RESULT([Using $NLOJETBINPATH for the binaries])
AC_MSG_RESULT([Using $NLOJETINCLUDEPATH for the includes])

#LHPDF (Attention! This is the version coming with NLOJET++ version 4!) 
AC_MSG_NOTICE([Looking for LHPDF installation ...])
AC_MSG_NOTICE([Attention! These are the pdf's coming with NLOJET++ version 4!])
AC_MSG_NOTICE([Do not confuse with LHAPDF which will also be required by fastNLO!])
if test "x$LHPDF" = "x"; then
  AC_MSG_ERROR([LHPDF environment variable is not set, configuration aborted!])
else
  AC_MSG_RESULT([Using environment variable LHPDF=$LHPDF to locate LHPDF installation])
  AC_SUBST([LHPDFPREFIX], "$LHPDF")
  AC_SUBST([LHPDFLIBPATH], "$LHPDF/lib")
  AC_SUBST([LHPDFINCLUDEPATH], "$LHPDFPREFIX/include/lhpdf") fi
AC_MSG_RESULT([Using $LHPDFPREFIX for the LHPDF installation path])
AC_MSG_RESULT([Using $LHPDFLIBPATH for the libraries])
AC_MSG_RESULT([Using $LHPDFINCLUDEPATH for the includes])

#LHAPDF
AC_MSG_NOTICE([Looking for LHAPDF installation ...])
if test "x$LHAPDF" = "x"; then
  AC_MSG_WARN([LHAPDF environment variable is not set, trying to find lhapdf-config in PATH ...])
  AC_PATH_PROG([LHAPDFBINPATH],[lhapdf-config],[:])
  if test "$LHAPDFBINPATH" = :; then
    AC_MSG_ERROR([No path to LHAPDF binaries found, configuration aborted!])
  else
    AC_SUBST([LHAPDFPREFIX], [`lhapdf-config --prefix`])
    AC_SUBST([LHAPDFSETPATH], [`lhapdf-config --pdfsets-path`])
  fi
else
  AC_MSG_RESULT([Using environment variable LHAPDF=$LHAPDF to locate LHAPDF installation])
  AC_SUBST([LHAPDFPREFIX], ["$LHAPDF"])
  AC_SUBST([LHAPDFSETPATH], [`$LHAPDF/bin/lhapdf-config --pdfsets-path`])
fi
AC_SUBST([LHAPDFBINPATH], ["$LHAPDFPREFIX/bin"])
AC_SUBST([LHAPDFLIBPATH], ["$LHAPDFPREFIX/lib"])
AC_SUBST([LHAPDFINCLUDEPATH],["$LHAPDFPREFIX/include/LHAPDF"])
AC_MSG_RESULT([Using $LHAPDFPREFIX for the LHAPDF installation path])
AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])
AC_MSG_RESULT([Using $LHAPDFBINPATH for the binaries])
AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
AC_MSG_RESULT([Using $LHAPDFINCLUDEPATH for the includes])

#CERNLIB
AC_MSG_NOTICE([Looking for CERN libraries ...])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_ERROR([CERN_ROOT environment variable is not set!
  CERN_ROOT/lib must be the path to -lkernlib, -lmathlib, -lpacklib])
else
  AC_MSG_RESULT([Using environment variable CERN_ROOT=$CERN_ROOT to locate CERNlib installation])
  AC_SUBST([CERNLIBPATH], ["$CERN_ROOT/lib"])
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using $CERN_ROOT for CERNlib installation path])
  AC_MSG_RESULT([Using $CERNLIBPATH for the libraries])
  AC_MSG_RESULT([Using -L$CERNLIBPATH -lmathlib -lkernlib -lpacklib for the linker flag])
fi

#ROOT
AC_MSG_NOTICE([Looking for ROOT installation ...])
if test "x$ROOTSYS" = "x"; then
  AC_MSG_WARN([ROOTSYS environment variable is not set, trying to find root-config in PATH ...])
  AC_PATH_PROG([ROOTBINPATH],[root-config],[:])
  if test "$ROOTBINPATH" = :; then
    AC_MSG_ERROR([No path to ROOT binaries found, configuration aborted!])
  else
    AC_SUBST([ROOTSYS], ["$ROOTBINPATH/.."])
    AC_SUBST([ROOTINCLUDEPATH], [`root-config --incdir`])
  fi
else
  AC_SUBST([ROOTSYS], ["$ROOTSYS"])	
  AC_SUBST([ROOTINCLUDEPATH], [`$ROOTSYS/bin/root-config --incdir`])
fi
AC_MSG_RESULT([Using $ROOTSYS for the ROOT installation path])
AC_MSG_RESULT([Using $ROOTINCLUDEPATH for the includes])
AC_CHECK_PROGS([have_rootcint],[rootcint],[:],["$ROOTSYS/bin/"])
if test "$have_rootcint" = ":"; then
  AC_MSG_WARN([ROOTSYS environment variable set and/or root-config found, but no rootcint. Not a big problem, but library will be build without ROOT dictionary.])
else
  AC_SUBST([ROOTCINTPATH], ["$ROOTSYS/bin/rootcint"])  
  AC_MSG_RESULT([Using $ROOTCINTPATH])
fi
AM_CONDITIONAL(HAVEROOT, test "x$ROOTCINTPATH" != "x")

# Output
AC_CONFIG_FILES([Makefile
                 author/hadron/Makefile
                 common/Makefile
                 merger/Makefile
                 user/Makefile
	         tableconvert/Makefile
	        ])
AC_OUTPUT
