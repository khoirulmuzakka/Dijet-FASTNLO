# Process this file with autoconf to produce a configure script.
AC_INIT([fastnlo], [2.0], [thomas.kluge@cern.ch,klaus.rabbertz@cern.ch])
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

# Announce checks particular for fastNLO
AC_MSG_NOTICE([******************************************])
AC_MSG_NOTICE([Starting checks particular for fastNLO ...])
AC_MSG_NOTICE([******************************************])

#NLOJET++
AC_MSG_NOTICE([Looking for NLOJET++ installation ...])
AC_MSG_NOTICE([First try to find NLOJET++ via create-nlojet-user executable in system PATH ...])
AC_PATH_PROG([NLOJETBINPATH],[create-nlojet-user],[no])
if test "x$NLOJETBINPATH" = "xno"; then
  AC_MSG_WARN([No path to NLOJET++ binaries found, try to use NLOJET environment variable ...])
  if test "x$NLOJET" = "x"; then
    AC_MSG_ERROR([NLOJET environment variable is not set, please install NLOJET++. Configuration aborted!])
  else
    AC_MSG_RESULT([Using environment variable NLOJET=$NLOJET to locate NLOJET++ installation ...])
    AC_PATH_PROG([NLOJETBINPATH],[create-nlojet-user],[no],[path=$NLOJET/bin])
    if test "x$NLOJETBINPATH" = "xno"; then
      AC_MSG_ERROR([NLOJET binaries not found in $NLOJET/bin, please install NLOJET++ correctly. Configuration aborted!])
    fi
  fi
fi
AC_SUBST([NLOJETBINPATH], [`dirname $NLOJETBINPATH`])
AC_SUBST([NLOJETPREFIX], [`dirname $NLOJETBINPATH`])
AC_SUBST([NLOJETLIBPATH], "$NLOJETPREFIX/lib")
AC_SUBST([NLOJETINCLUDEPATH], "$NLOJETPREFIX/include/nlo++")
AC_MSG_RESULT([Using $NLOJETPREFIX for the NLOJET++ installation path])
AC_MSG_RESULT([Using $NLOJETBINPATH for the binaries])
AC_MSG_RESULT([Using $NLOJETLIBPATH for the libraries])
AC_MSG_RESULT([Using $NLOJETINCLUDEPATH for the includes])

#LHPDF (Attention! This is the version coming with NLOJET++ version 4!) 
AC_MSG_NOTICE([Looking for LHPDF installation ...])
AC_MSG_NOTICE([Attention! These are the pdf's coming with NLOJET++ version 4!])
AC_MSG_NOTICE([Do not confuse with LHAPDF which will also be required by fastNLO!])
# Does not work yet
AC_MSG_NOTICE([First try to find lhpdf in system libs or LD_LIBRARY_PATH ...])
saved_LIBS="$LIBS"
saved_LDFLAGS="$LDFLAGS"
found=no
for i in "yes" $(echo "$LD_LIBRARY_PATH" | tr : ' '); do
  if test "x$i" = "xyes"; then
    LIBS="$saved_LIBS"
    LDFLAGS="$saved_LDFLAGS"
  else
    LIBS="$saved_LIBS -l$i"
    LDFLAGS="$saved_LDFLAGS -L$i"
  fi
  echo "LIBS $LIBS"
  echo "LDFLAGS $LDFLAGS"
#  AC_CHECK_LIB([lhpdf], [lhpdf], [found="$i"; break])
  AC_SEARCH_LIBS([lhpdf],[lhpdf])
done

AC_SUBST([LHPDFPREFIX], "")
if test "x$found" = "xyes"; then
  AC_MSG_RESULT([Found lhpdf lib amongst system libs, assuming /usr as system install path])
  AC_SUBST([LHPDFPREFIX], "/usr")
else
  if test "x$found" = "xno"; then
    AC_MSG_WARN([lhpdf lib not found, try to use LHPDF environment variable ...])
    if test "x$LHPDF" = "x"; then
      AC_MSG_ERROR([LHPDF environment variable is not set, please install LHPDF. Configuration aborted!])
    else
      AC_MSG_RESULT([Using environment variable LHPDF=$LHPDF to locate LHPDF installation ...])
#      LIBS="$saved_LIBS -L$LHPDF/lib"
#      LDFLAGS="$saved_LDFLAGS -L$LHPDF/lib"
#      AC_CHECK_LIB([lhpdf], [lhpdf_pdf_evolveas], [found="$LHPDF/lib"])
      found=$LHPDF/lib
      if test "x$found" = "xno"; then
        AC_MSG_ERROR([lhpdf lib not found in $LHPDF/lib, configuration aborted])
      else
        AC_SUBST([LHPDFPREFIX], "$LHPDF")
      fi
    fi
  else
    AC_MSG_RESULT([Found lhpdf lib in $found from LD_LIBRARY_PATH])
    AC_SUBST([LHPDFPREFIX], `dirname $found`)
  fi
fi
LIBS="$saved_LIBS"
LDFLAGS="$saved_LDFLAGS"
AC_SUBST([LHPDFLIBPATH], "$LHPDFPREFIX/lib")
AC_SUBST([LHPDFINCLUDEPATH], "$LHPDFPREFIX/include/lhpdf")
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

#FASTJET
AC_MSG_NOTICE([Looking for fastjet installation ...])
if test "x$FASTJET" = "x"; then
  AC_MSG_WARN([FASTJET environment variable is not set, trying to find fastjet-config in PATH ...])
  AC_PATH_PROG([FASTJETBINPATH],[fastjet-config],[:])
  if test "$FASTJETBINPATH" = :; then
    AC_MSG_ERROR([No path to FASTJET binaries found, configuration aborted!])
  else
    AC_SUBST([FASTJETPREFIX], [`fastjet-config --prefix`])
  fi
else
  AC_MSG_RESULT([Using environment variable FASTJET=$FASTJET to locate FASTJET installation])
  AC_SUBST([FASTJETPREFIX], ["$FASTJET"])
fi
AC_SUBST([FASTJETBINPATH], ["$FASTJETPREFIX/bin"])
AC_SUBST([FASTJETLIBPATH], ["$FASTJETPREFIX/lib"])
AC_SUBST([FASTJETINCLUDEPATH],["$FASTJETPREFIX/include"])
AC_MSG_RESULT([Using $FASTJETPREFIX for the FASTJET installation path])
AC_MSG_RESULT([Using $FASTJETBINPATH for the binaries])
AC_MSG_RESULT([Using $FASTJETLIBPATH for the library])
AC_MSG_RESULT([Using $FASTJETINCLUDEPATH for the includes])

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
    AC_SUBST([ROOTBINPATH], [`dirname $ROOTBINPATH`])
    AC_SUBST([ROOTSYS], [`dirname $ROOTBINPATH`])
    AC_SUBST([ROOTLIBPATH], [`root-config --libdir`])
    AC_SUBST([ROOTINCLUDEPATH], [`root-config --incdir`])
  fi
else
  AC_SUBST([ROOTSYS], ["$ROOTSYS"])	
  AC_SUBST([ROOTBINPATH], [`$ROOTSYS/bin/root-config --bindir`])
  AC_SUBST([ROOTLIBPATH], [`$ROOTSYS/bin/root-config --libdir`])
  AC_SUBST([ROOTINCLUDEPATH], [`$ROOTSYS/bin/root-config --incdir`])
fi
AC_MSG_RESULT([Using $ROOTSYS for the ROOT installation path])
AC_MSG_RESULT([Using $ROOTBINPATH for the binaries])
AC_MSG_RESULT([Using $ROOTINCLUDEPATH for the includes])
AC_MSG_RESULT([Using $ROOTLIBPATH for the libraries])
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
