#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.10.2013
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Include files to be distributed
#
pkginclude_HEADERS = fastNLOBase.h fastNLOCoeffAddBase.h fastNLOCoeffAddFix.h fastNLOCoeffAddFlex.h fastNLOEvent.h \
	fastNLOCoeffBase.h fastNLOCoeffData.h fastNLOCoeffMult.h fastNLOCoefficients.h fastNLOConstants.h fastNLOGeneratorConstants.h \
	fastNLOCreate.h fastNLOInterpolBase.h fastNLOInterpolCatmullRom.h fastNLOInterpolLagrange.h fastNLOInterpolLinear.h fastNLOInterpolOneNode.h \
	fastNLOAlphas.h fastNLOLHAPDF.h fastNLOPDFLinearCombinations.h fastNLOTools.h \
	fastNLODiffReader.h fastNLODiffUser.h \
	fastNLOReader.h fastNLOTable.h read_steer.h speaker.h Alphas.h CRunDec.h fastNLOCRunDec.h fastNLOQCDNUMAS.h fastNLOHoppet.h

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
