# Determine Rivet analysis name
MAKEFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
ANALYSIS := MCgrid_$(notdir $(patsubst %/,%,$(dir $(MAKEFILEPATH))))

# Rivet flags
RFLAGS=$(shell rivet-config --cppflags)
RLDFLAGS=$(shell rivet-config --ldflags)
RDATA=$(shell rivet-config --datadir)
RLIBS=$(shell rivet-config --libdir)

# APPLgrid flags
# The ":=" must be there to have a single-expanded make variable!
APPLVERSION:=$(shell applgrid-config --version 2>/dev/null)
APPLFLAGS=$(shell applgrid-config --cxxflags)
APPLLDFLAGS=$(shell applgrid-config --ldflags)
APPLSHARE=$(shell applgrid-config --share)

# fastNLO flags
# The ":=" must be there to have a single-expanded make variable!
FNLOVERSION:=$(shell fnlo-tk-config --version 2>/dev/null)
FNLOFLAGS=$(shell fnlo-tk-config --cxxflags)
FNLOLDFLAGS=$(shell fnlo-tk-config --ldflags)
FNLOSHARE=$(shell fnlo-tk-config --share)

# MCgrid flags
MCFLAGS=$(shell pkg-config mcgrid --cflags)
MCLDFLAGS=$(shell pkg-config mcgrid --libs)

# combine flags
PLGFLAGS=$(RFLAGS) $(AFLAGS) $(MCFLAGS)
PLGLDFLAGS=$(RLDFLAGS) $(ALDFLAGS) $(MCLDFLAGS)

# plugin input/output and data
RTARGET = analysis/Rivet$(ANALYSIS)Plugins.so
RSOURCES = $(wildcard ./analysis/*.cc)
RYODAFILE = $(wildcard ./analysis/*.yoda)
RDATAFILES = $(wildcard ./analysis/*.info) $(wildcard ./analysis/*.plot) $(RYODAFILE)

all: plugin-applgrid plugin-fastnlo

.PHONY: all plugin-applgrid plugin-fastnlo clean install install-applgrid

plugin-applgrid:
ifdef APPLVERSION
	@echo "Building APPLgrid plugin with version $(APPLVERSION)"
	rivet-buildplugin $(RTARGET) $(RSOURCES) $(PLGFLAGS) $(PLGLDFLAGS) -DUSE_APPL=1
else
	@echo "APPLgrid not found, no plugin built!"
endif

plugin-fastnlo:
ifdef FNLOVERSION
	@echo "Building fastNLO plugin with version $(FNLOVERSION)"
	rivet-buildplugin $(RTARGET) $(RSOURCES) $(PLGFLAGS) $(PLGLDFLAGS) -DUSE_FNLO=1
else
	@echo "fastNLO not found, no plugin built!"
endif

clean:
	rm -f $(RTARGET)

install:
	cp $(RDATAFILES) $(RDATA)
	cp $(RTARGET) $(RLIBS)

install-applgrid:
	cp ./subprocesses/*.config $(ASHARE)

# Temporarily hide the data file before plotting as rivet-mkhtml will otherwise normalise against data
# which is not very useful if we want to compare the tiny deviations between grid and generator
# (on the other hand, data must be available for the Rivet init phase, so it has to be installed)
comparison-plot:
	rm -f $(RDATA)/$(ANALYSIS).yoda
	rivet-mkhtml *yoda
	cp $(RYODAFILE) $(RDATA)
