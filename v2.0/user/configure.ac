# Process this file with autoconf to produce a configure script.
AC_INIT([user], [1.0], [wobisch@fnal.gov])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_F77
AC_PROG_CPP

# Checks for sources.
AC_CONFIG_SRCDIR([fnx9999.f])

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS

#CERNLIB
AC_MSG_CHECKING([for CERN libraries])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_ERROR([CERN_ROOT environment variable is not set!
  CERN_ROOT/lib must be the path to -lkernlib, -lmathlib, -lpacklib])
else
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib -lkernlib -lpacklib ])
fi

#LHAPDF
AC_CHECK_PROGS([have_lhapdf_config],[lhapdf-config],[:])
if test "$have_lhapdf_config" = :; then
   AC_MSG_ERROR([lhapdf-config is not in path! This means LHAPDF is not properly installed.])
fi
AC_MSG_CHECKING([for LHAPDF installation])
AC_SUBST([LHAPDFPREFIX], [`lhapdf-config --prefix`])
AC_SUBST([LHAPDFLIBPATH],$LHAPDFPREFIX/lib)
AC_SUBST([LHAPDFSETPATH], [`lhapdf-config --pdfsets-path`])
AC_MSG_RESULT([Using $LHAPDFLIBPATH for the library])
AC_MSG_RESULT([Using $LHAPDFSETPATH for the PDF sets])



# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
