#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 13.12.2011
#
# Define common libraries
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main libraries with the fastNLO classes inside
pkglib_LTLIBRARIES = libfastnlo.la libfastnloJetAlgos.la
# Needed for direct use in ROOT
# libfastnloUser.la

# The associated header files
pkginclude_HEADERS = fnloTable.h fnloTableUser.h LinkDef.h fnloBlockBNlojet.h

libfastnlo_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc fnloBlockBNlojet.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h fnloBlockA2.h fnloBlockB.h fnloBlockBNlojet.h \
  fnloCmdLine.h fnloCmdLine.cc fnloUtils.h fnloUtils.hxx fnloUtils.cc

#libfastnloUser_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc \
#  fnloTableUser.cc \
#  fnloconstants.h fnloTable.h fnloBlockA1.h fnloBlockA2.h fnloBlockB.h fnloTableUser.h

# If ROOT is installed, add a dictionary to the library
#if HAVE_ROOT
# Note: fnloTableDict.cc has to be pre-created and therefore is a BUILT_SOURCE !
#  BUILT_SOURCES = fnloTableDict.cc
# But don't distribute it ...
#  nodist_libfastnlo_la_SOURCES = fnloTableDict.cc
#  nodist_libfastnloUser_la_SOURCES = fnloTableDict.cc
# Only these go to the distribution ...
#  libfastnlo_la_SOURCES += fnloTable.h fnloTableUser.h LinkDef.h
#  libfastnloUser_la_SOURCES += fnloTable.h fnloTableUser.h LinkDef.h
# With my ROOT installation the "pre-compile" command with ROOTCINT fails
# with seg fault :-(
# Needs to be checked/fixed
#fnloTableDict.cc: fnloTable.h fnloTableUser.h LinkDef.h
#	${ROOTCINT} -f fnloTableDict.cc -c -I${ROOTINCLUDEPATH} -I${FASTNLOINCLUDEPATH} fnloTable.h fnloTableUser.h LinkDef.h
#endif
# Fake creation to make "make dist" happy
#fnloTableDict.cc: fnloTable.h fnloTableUser.h LinkDef.h
#	touch fnloTableDict.cc

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Preprocessor flags
AM_CPPFLAGS = -I${NLOJETINCLUDEPATH} -I${LHAPDFINCLUDEPATH}
if HAVE_FASTJET
  AM_CPPFLAGS += -I${FASTJETINCLUDEPATH}

# Source files to be distributed for each library
  libfastnloJetAlgos_la_SOURCES = cone-e.cc cone-e.h cone-et.cc cone-et.h \
    fj-ak.cc fj-ak.h fj-d0r2.cc fj-d0r2.h fj-kt.cc fj-kt.h fj-sc.cc fj-sc.h \
    kt-e.cc kt-e.h

endif
if HAVE_ROOT
  AM_CPPFLAGS += -I${ROOTINCLUDEPATH}
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = fnloTableDict.h fnloTableDict.cc \
	configure Makefile.in aclocal.m4 depcomp missing install-sh
