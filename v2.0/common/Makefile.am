AUTOMAKE_OPTIONS = foreign

# the main library with the fastNLO classes inside
pkglib_LTLIBRARIES = libfastnlo.la libfastnloUser.la libfastnloNlojet.la

libfastnlo_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc  \
  fnloconstants.h fnloTable.h fnloBlockA1.h  fnloBlockA2.h  fnloBlockB.h

libfastnloUser_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc  fnloTableUser.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h  fnloBlockA2.h  fnloBlockB.h fnloTableUser.h LHAPDF.h

libfastnloNlojet_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc fnloBlockBNlojet.cc  \
  fnloconstants.h fnloTable.h fnloBlockA1.h  fnloBlockA2.h  fnloBlockB.h fnloBlockBNlojet.h

# if ROOT is installed, add a dictionary to the library
if HAVEROOT
  libfastnlo_la_SOURCES += fnloTableDict.cc
  libfastnloUser_la_SOURCES += fnloTableDict.cc
  fnloTableDict.cc: fnloTable.h fnloTableUser.h LinkDef.h
	-${ROOTCINTPATH} -f fnloTableDict.cc -c -I${LHAPDFINCLUDEPATH} fnloTable.h fnloTableUser.h LinkDef.h
endif

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Preprocessor flags
AM_CPPFLAGS = -I${NLOJETINCLUDEPATH} -I${LHAPDFINCLUDEPATH}
if HAVEROOT
  AM_CPPFLAGS += -I${ROOTINCLUDEPATH}
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = fnloTableDict.h fnloTableDict.cc \
	configure Makefile.in aclocal.m4 depcomp missing install-sh
