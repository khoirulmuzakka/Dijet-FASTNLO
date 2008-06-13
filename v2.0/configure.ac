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

#CERNLIB
AC_MSG_CHECKING([for CERN libraries])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_WARN([CERN_ROOT environment variable is not set!
  CERN_ROOT/lib must be the path to -lkernlib, -lmathlib, -lpacklib])
else
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib ])
fi

#LHAPDF
AC_CHECK_PROGS([have_lhapdf_config],[lhapdf-config],[:])
if test "$have_lhapdf_config" = :; then
   AC_MSG_WARN([lhapdf-config is not in path! This means LHAPDF is not properly installed.])
fi
AC_MSG_CHECKING([for LHAPDF installation])
AC_SUBST([LHAPDFPREFIX], [`lhapdf-config --prefix`])
AC_SUBST([LHAPDFLIBPATH],$LHAPDFPREFIX/lib)
AC_SUBST([LHAPDFSETPATH], [`lhapdf-config --pdfsets-path`])
AC_SUBST([LHAPDFINCLUDEPATH],"$LHAPDFPREFIX/include/lhapdf")
AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
AC_MSG_RESULT([Using $LHAPDFINCLUDEPATH for the includes])
AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])

#ROOT
AC_MSG_CHECKING([for ROOT installation])
if test "x$ROOTSYS" = "x"; then
  AC_MSG_WARN([ROOTSYS environment variable is not set.
  Not a big problem, but library will be build without ROOT dictionary.])
else
  AC_SUBST([ROOTSYS], ["$ROOTSYS"])	
  AC_CHECK_PROGS([have_rootcint],[rootcint],[:],["$ROOTSYS/bin/"])
  if test "$have_rootcint" = ":"; then
    AC_MSG_WARN([ROOTSYS environment variable set, but no rootcint found.
    Not a big problem, but library will be build without ROOT dictionary.])
  else
    AC_SUBST([ROOTCINTPATH], ["$ROOTSYS/bin/rootcint"])  
    AC_MSG_RESULT([Using $ROOTCINTPATH])
  fi
fi
AM_CONDITIONAL(HAVEROOT, test "x$ROOTCINTPATH" != "x")

# Output
AC_CONFIG_FILES([Makefile
                 common/Makefile
                 merger/Makefile
                 user/Makefile
	         tableconvert/Makefile
	        ])
AC_OUTPUT
