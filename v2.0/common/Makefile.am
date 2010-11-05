# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main libraries with the fastNLO classes inside
pkglib_LTLIBRARIES = libfastnlo.la libfastnloUser.la libfastnloNlojet.la

# The associated header files
# Add header files for jet algos ?
pkginclude_HEADERS = fnloTable.h fnloBlockBNlojet.h cone-e-07.h

# Source files to be distributed for each library
libfastnlo_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h fnloBlockA2.h fnloBlockB.h

libfastnloUser_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc \
  fnloTableUser.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h fnloBlockA2.h fnloBlockB.h fnloTableUser.h \
  LHAPDF.h LHAPDFConfig.h
# Copied from LHAPDF 5.8.4. This is far from optimal ...

libfastnloNlojet_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc \
  fnloBlockBNlojet.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h fnloBlockA2.h fnloBlockB.h fnloBlockBNlojet.h
# If ROOT is installed, add a dictionary to the library
#if HAVEROOT
#  libfastnlo_la_SOURCES += fnloTableDict.cc
#  libfastnloUser_la_SOURCES += fnloTableDict.cc
#  fnloTableDict.cc: fnloTable.h fnloTableUser.h LinkDef.h
#	-${ROOTCINTPATH} -f fnloTableDict.cc -c -I${LHAPDFINCLUDEPATH} fnloTable.h fnloTableUser.h LinkDef.h
#endif

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
