# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main targets, here the libtool plugin libraries to build
lib_LTLIBRARIES = libfnh1001.la libfnh5001.la libfnh5002.la libfnh5002aemHS2012.la libfnh5003ak.la libfnh5003.la libfnhd5102.la libfnh4002.la libfnh5721.la
#lib_LTLIBRARIES += libfnhd1002.la libfnhd1012.la libfnhd1022.la libfnhd0001.la
lib_LTLIBRARIES += libfnh5201.la libfnh5202.la libfnh4401.la
lib_LTLIBRARIES += libfnh5001PhiCut.la
lib_LTLIBRARIES += libfnh5001SucheMuRMean.la
lib_LTLIBRARIES += libfnh50031SucheMuRMean.la
lib_LTLIBRARIES += libfnh5301kt.la libfnh5302kt.la libfnh5303kt.la libfnh5304kt.la 
lib_LTLIBRARIES += libfnh5821_Incl.la libfnh5821_Dij.la
lib_LTLIBRARIES += libfnh5001v22.la
# libfnh2101d.la

#
# All source files which go into the main targets (and have to be distributed)
#
# All targets
AM_SOURCES = pdf-cteq6.h pdf-dis-dummy.h
# Individual sources
libfnh1001_la_SOURCES = fnh1001.cc alew.f $(AM_SOURCES)
libfnh5001_la_SOURCES = fnh5001.cc alew.f $(AM_SOURCES)
libfnh5002_la_SOURCES = fnh5002.cc alew.f $(AM_SOURCES)
libfnh5002aemHS2012_la_SOURCES = fnh5002aemHS2012.cc wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)
libfnh5003_la_SOURCES = fnh5003.cc alew.f $(AM_SOURCES)
libfnh5003ak_la_SOURCES = fnh5003ak.cc alew.f $(AM_SOURCES)
libfnh5721_la_SOURCES = fnh5721.cc alew.f $(AM_SOURCES)
#libfnhd5102_la_SOURCES = fnhd5102.cc alew.f $(AM_SOURCES)
#libfnhd0001_la_SOURCES = fnhd0001.cc alew.f $(AM_SOURCES)
#libfnhd1002_la_SOURCES = fnhd1002.cc alew.f $(AM_SOURCES)
#libfnhd1012_la_SOURCES = fnhd1012.cc h1jet.cc alew.f qcd_2006.f i_2006_fitb.f i_2006_fita.f h12006flux.f h12007jetsdpdf_pdflib_interface.F $(AM_SOURCES)
#libfnhd1022_la_SOURCES = fnhd1022.cc alew.f $(AM_SOURCES)
libfnh4002_la_SOURCES = fnh4002.cc alew.f $(AM_SOURCES)
libfnh4401_la_SOURCES = fnh4401.cc alew.f $(AM_SOURCES)
libfnh5201_la_SOURCES = fnh5201.cc alew.f $(AM_SOURCES)
libfnh5202_la_SOURCES = fnh5202.cc alew.f $(AM_SOURCES)
libfnh5001PhiCut_la_SOURCES = fnh5001PhiCut.cc alew.f $(AM_SOURCES)
libfnh5001SucheMuRMean_la_SOURCES = fnh5001SucheMuRMean.cc alew.f $(AM_SOURCES)
libfnh50031SucheMuRMean_la_SOURCES = fnh50031SucheMuRMean.cc alew.f $(AM_SOURCES)
libfnh5821_Incl_la_SOURCES = fnh5821_Incl.cc alew.f $(AM_SOURCES)
libfnh5821_Dij_la_SOURCES = fnh5821_Dij.cc alew.f $(AM_SOURCES)

libfnh5301kt_la_SOURCES = fnh5301kt.cc alew.f wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)
libfnh5302kt_la_SOURCES = fnh5302kt.cc alew.f wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)
libfnh5304kt_la_SOURCES = fnh5304kt.cc alew.f wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)
libfnh5303kt_la_SOURCES = fnh5303kt.cc alew.f wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)

libfnh5001v22_la_SOURCES = fnh5001v22.cc alew.f wrap_ew.f hadr5n12.f formff.f $(AM_SOURCES)

#
# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
AM_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${FASTJETINCLUDEPATH} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
# Example of individual specification of CXXFLAGS
#libfnh1001_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
#libfnh2101d_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

#
# Fortran compiler flags to be used for the main target
#
# All targets
AM_FFLAGS = -O0 -Wall
# Example of individual specification of FFLAGS
#libfnh1001_la_FFLAGS = -O0 -Wall
#libfnh2101d_la_FFLAGS = -O0 -Wall

#
# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
#AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPAT#H} -llhpdf_cteq6
AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
# Example of individual specification of LDFLAGS
libfnh1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnh5001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lfastnlotoolkit
libfnh5001v22_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lfastnlotoolkit
libfnh5002_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnh5003ak_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet
#libfnh2101d_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution.
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
