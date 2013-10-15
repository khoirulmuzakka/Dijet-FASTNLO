#Rabbertz, 10.11.2011
#
## fastNLO_reader: C++ version~
#
## -----------------------------------------------------------------------------

if ENABLE_PYEXT

EXTRA_DIST = fastnloreader.i fastnloreader_wrap.cpp fastnloreader.py
fastnloreader_wrap.cpp fastnloreader.py: fastnloreader.i
	swig -c++ -python -I$(top_srcdir)/include -o fastnloreader_wrap.cpp $<

AM_CXXFLAGS = -Wno-long-long

MAINTAINERCLEANFILES = fastnloreader.py fastnloreader_wrap.cpp

all-local: fastnloreader_wrap.cpp
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/lhapdf*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc

endif

