#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 22.11.2012
#
# fastNLO_reader: C++ version
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Include files to be distributed
#
fastnloincludedir = $(includedir)/fastnlo
fastnloinclude_HEADERS = FastNLOReader.h FastNLODiffReader.h FastNLODiffUser.h FastNLOBlockB.h Alphas.h speaker.h \
	FastNLOUser.h FastNLOAlphas.h FastNLOLHAPDF.h FastNLOCRunDec.h CRunDec.h FastNLOQCDNUMAS.h FastNLOHoppet.h

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
