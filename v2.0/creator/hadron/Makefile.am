# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main targets, here the libtool plugin libraries to build
# TBD: Better select the "main" targets ...
# TBD: Uncomment here or below to install additional scenarios
lib_LTLIBRARIES = libfnt1001midp.la libfnt2009midp.la libfnt2012test.la libfnt2013norm.la libfnt2014norm.la libfnl1001.la
#
# Require fastjet
#
if HAVE_FASTJET
#   lib_LTLIBRARIES += libfnl2222b.la
   lib_LTLIBRARIES += libfnl2342b.la libfnl2352.la
   lib_LTLIBRARIES += libfnl2412c.la libfnl2412e.la libfnl2422.la
#   lib_LTLIBRARIES += libfnl2452.la libfnl2822num.la libfnl2822den.la libfnl1001ak07.la libfnl1001sc07.la libfnl1001d0r207.la libfnt2012.la
   lib_LTLIBRARIES += libfnl2622f.la libfnl3622g.la
   lib_LTLIBRARIES += libfnl2332y0.la libfnl2332y1.la libfnl2332y2.la libfnl2332y3.la libfnl2332y4.la libfnl2332y5.la libfnl2332y6.la libfnl2332d.la
   lib_LTLIBRARIES += libfnl3332y0.la libfnl3332y1.la libfnl3332y2.la libfnl3332y3.la libfnl3332y4.la libfnl3332y5.la libfnl3332y6.la
   lib_LTLIBRARIES += libfnl4332y0.la libfnl4332y1.la libfnl4332y2.la libfnl4332y3.la libfnl4332y4.la libfnl4332y5.la libfnl4332y6.la libfnl4350.la
   lib_LTLIBRARIES += libfnl5002ptavystar.la libfnl5002pt12.la
   lib_LTLIBRARIES += libfnl2732dena.la libfnl2732denb.la libfnl2732denk.la libfnl2732denl.la
   lib_LTLIBRARIES += libfnl2732numa.la libfnl2732numb.la libfnl2732numk.la libfnl2732numl.la
#   lib_LTLIBRARIES += libfnl2912b.la libfnl2912h.la libfnl2912bptavm3.la
   lib_LTLIBRARIES += libfnl2912by0ptav.la libfnl2912by0m3.la libfnl2912by0ptavm3.la libfnl2912hy0ptav.la libfnl2912hy0m3.la libfnl2912hy0ptavm3.la libfnl2912by1ptav.la libfnl2912by1m3.la libfnl2912by1ptavm3.la libfnl2912hy1ptav.la libfnl2912hy1m3.la libfnl2912hy1ptavm3.la libfnl2912by2ptav.la
   lib_LTLIBRARIES += libfnl2932by0ptav.la libfnl2932by0m3.la libfnl2932by1ptav.la libfnl2932by1m3.la
   lib_LTLIBRARIES += libfnl2370.la
endif

#
# All source files which go into the main targets (and have to be distributed)
# Note: All standard jet algo interfaces like fj-ak.cc and fj-ak.h reside in the
#       "common" directory and are collected into the common library libfastnloJetAlgos.
#       An example how to add one locally here can be found below.
#
# All targets
AM_SOURCES = pdf-cteq6.h pdf-hhc-dummy.h
# cteq6.h is included from lhpdf installation which is different from original
# Individual sources
libfnt1001midp_la_SOURCES      = fnt1001midp.cc      $(AM_SOURCES)
libfnt2009midp_la_SOURCES      = fnt2009midp.cc      $(AM_SOURCES)
libfnt2012test_la_SOURCES      = fnt2012test.cc      $(AM_SOURCES)
libfnt2013norm_la_SOURCES      = fnt2013norm.cc      $(AM_SOURCES)
libfnt2014norm_la_SOURCES      = fnt2014norm.cc      $(AM_SOURCES)
#libfnt2015norm_la_SOURCES      = fnt2015norm.cc      $(AM_SOURCES)
libfnl1001_la_SOURCES          = fnl1001.cc          $(AM_SOURCES)
# Require fastjet
libfnl1001ak07_la_SOURCES      = fnl1001ak07.cc      $(AM_SOURCES)
libfnl1001sc07_la_SOURCES      = fnl1001sc07.cc      $(AM_SOURCES)
libfnl1001d0r207_la_SOURCES    = fnl1001d0r207.cc    $(AM_SOURCES)
libfnt2012_la_SOURCES          = fnt2012.cc          $(AM_SOURCES)
libfnl2222b_la_SOURCES         = fnl2222b.cc         $(AM_SOURCES)
libfnl2342b_la_SOURCES         = fnl2342b.cc         $(AM_SOURCES)
libfnl2352_la_SOURCES          = fnl2352.cc          $(AM_SOURCES)
libfnl2452_la_SOURCES          = fnl2452.cc          $(AM_SOURCES)
libfnl2412c_la_SOURCES         = fnl2412c.cc         $(AM_SOURCES)
libfnl2412e_la_SOURCES         = fnl2412e.cc         $(AM_SOURCES)
libfnl2422_la_SOURCES          = fnl2422.cc          $(AM_SOURCES)
libfnl2622f_la_SOURCES         = fnl2622f.cc         $(AM_SOURCES)
libfnl2732dena_la_SOURCES      = fnl2732dena.cc      $(AM_SOURCES)
libfnl2732denb_la_SOURCES      = fnl2732denb.cc      $(AM_SOURCES)
libfnl2732denk_la_SOURCES      = fnl2732denk.cc      $(AM_SOURCES)
libfnl2732denl_la_SOURCES      = fnl2732denl.cc      $(AM_SOURCES)
libfnl2732numa_la_SOURCES      = fnl2732numa.cc      $(AM_SOURCES)
libfnl2732numb_la_SOURCES      = fnl2732numb.cc      $(AM_SOURCES)
libfnl2732numk_la_SOURCES      = fnl2732numk.cc      $(AM_SOURCES)
libfnl2732numl_la_SOURCES      = fnl2732numl.cc      $(AM_SOURCES)
libfnl2822den_la_SOURCES       = fnl2822den.cc       $(AM_SOURCES)
libfnl2822num_la_SOURCES       = fnl2822num.cc       $(AM_SOURCES)
libfnl2912b_la_SOURCES         = fnl2912b.cc         $(AM_SOURCES)
libfnl2912by0ptav_la_SOURCES   = fnl2912by0ptav.cc   $(AM_SOURCES)
libfnl2912by0m3_la_SOURCES     = fnl2912by0m3.cc     $(AM_SOURCES)
libfnl2912by1ptav_la_SOURCES   = fnl2912by1ptav.cc   $(AM_SOURCES)
libfnl2912by1m3_la_SOURCES     = fnl2912by1m3.cc     $(AM_SOURCES)
libfnl2912by2ptav_la_SOURCES   = fnl2912by2ptav.cc   $(AM_SOURCES)
libfnl2912hy0ptav_la_SOURCES   = fnl2912hy0ptav.cc   $(AM_SOURCES)
libfnl2912hy0m3_la_SOURCES     = fnl2912hy0m3.cc     $(AM_SOURCES)
libfnl2912hy1ptav_la_SOURCES   = fnl2912hy1ptav.cc   $(AM_SOURCES)
libfnl2912hy1m3_la_SOURCES     = fnl2912hy1m3.cc     $(AM_SOURCES)
libfnl2932by0ptav_la_SOURCES   = fnl2932by0ptav.cc   $(AM_SOURCES)
libfnl2932by0m3_la_SOURCES     = fnl2932by0m3.cc     $(AM_SOURCES)
libfnl2932by1ptav_la_SOURCES   = fnl2932by1ptav.cc   $(AM_SOURCES)
libfnl2932by1m3_la_SOURCES     = fnl2932by1m3.cc     $(AM_SOURCES)
libfnl2912h_la_SOURCES         = fnl2912h.cc         $(AM_SOURCES)
libfnl2332y0_la_SOURCES        = fnl2332y0.cc        $(AM_SOURCES)
libfnl2332y1_la_SOURCES        = fnl2332y1.cc        $(AM_SOURCES)
libfnl2332y2_la_SOURCES        = fnl2332y2.cc        $(AM_SOURCES)
libfnl2332y3_la_SOURCES        = fnl2332y3.cc        $(AM_SOURCES)
libfnl2332y4_la_SOURCES        = fnl2332y4.cc        $(AM_SOURCES)
libfnl2332y5_la_SOURCES        = fnl2332y5.cc        $(AM_SOURCES)
libfnl2332y6_la_SOURCES        = fnl2332y6.cc        $(AM_SOURCES)
libfnl2332d_la_SOURCES         = fnl2332d.cc         $(AM_SOURCES)
libfnl3332y0_la_SOURCES        = fnl3332y0.cc        $(AM_SOURCES)
libfnl3332y1_la_SOURCES        = fnl3332y1.cc        $(AM_SOURCES)
libfnl3332y2_la_SOURCES        = fnl3332y2.cc        $(AM_SOURCES)
libfnl3332y3_la_SOURCES        = fnl3332y3.cc        $(AM_SOURCES)
libfnl3332y4_la_SOURCES        = fnl3332y4.cc        $(AM_SOURCES)
libfnl3332y5_la_SOURCES        = fnl3332y5.cc        $(AM_SOURCES)
libfnl3332y6_la_SOURCES        = fnl3332y6.cc        $(AM_SOURCES)
libfnl3622g_la_SOURCES         = fnl3622g.cc         $(AM_SOURCES)
libfnl4332y0_la_SOURCES        = fnl4332y0.cc        $(AM_SOURCES)
libfnl4332y1_la_SOURCES        = fnl4332y1.cc        $(AM_SOURCES)
libfnl4332y2_la_SOURCES        = fnl4332y2.cc        $(AM_SOURCES)
libfnl4332y3_la_SOURCES        = fnl4332y3.cc        $(AM_SOURCES)
libfnl4332y4_la_SOURCES        = fnl4332y4.cc        $(AM_SOURCES)
libfnl4332y5_la_SOURCES        = fnl4332y5.cc        $(AM_SOURCES)
libfnl4332y6_la_SOURCES        = fnl4332y6.cc        $(AM_SOURCES)
libfnl4350_la_SOURCES          = fnl4350.cc          $(AM_SOURCES)
libfnl2370_la_SOURCES          = fnl2370.cc          $(AM_SOURCES)
# v21
libfnl5002ptavystar_la_SOURCES = fnl5002ptavystar.cc $(AM_SOURCES)
libfnl5002pt12_la_SOURCES      = fnl5002pt12.cc      $(AM_SOURCES)
libfnl2912bptavm3_la_SOURCES   = fnl2912bptavm3.cc   $(AM_SOURCES)
libfnl2912by0ptavm3_la_SOURCES = fnl2912by0ptavm3.cc $(AM_SOURCES)
libfnl2912hy0ptavm3_la_SOURCES = fnl2912hy0ptavm3.cc $(AM_SOURCES)
libfnl2912by1ptavm3_la_SOURCES = fnl2912by1ptavm3.cc $(AM_SOURCES)
libfnl2912hy1ptavm3_la_SOURCES = fnl2912hy1ptavm3.cc $(AM_SOURCES)

# Example how to add a jet algo interface locally in this directory
# Should only be done temporarily
#libfnl1001ak07_la_SOURCES   = fnl1001ak07.cc   fj-ak.h      fj-ak.cc      $(AM_SOURCES)

#
# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
#AM_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
if HAVE_FASTJETCONFIG
   AM_CXXFLAGS += `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
else
if HAVE_FASTJET
   AM_CXXFLAGS += -I${FASTJETINCLUDEPATH}
endif
endif

# Example of individual specification of CXXFLAGS
#libfnl1001_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
# Require fastjet, use fastjet-config if possible
#if HAVE_FASTJETCONFIG
#   libfnl1001sc07_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
#else
#   libfnl1001sc07_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}
#endif

#
# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
#AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -lfastnloJetAlgo
AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -lfastnloJetAlgos -L${LHPDFLIBPATH} -llhpdf_cteq6
if HAVE_FASTJETCONFIG
   AM_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs --plugins`
else
if HAVE_FASTJET
   AM_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
endif
endif

# Example of individual specification of LDFLAGS
#libfnl1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
# Require fastjet, use fastjet-config if possible
#if HAVE_FASTJETCONFIG
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 `${FASTJETBINPATH}/fastjet-config --libs --plugins`
#else
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
#endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
