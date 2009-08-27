AUTOMAKE_OPTIONS = foreign
SUBDIRS = common merger user tableconvert author/hadron author/ep

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = config.guess config.status config.sub configure \
	Makefile.in Makefile aclocal.m4 auto4mte.cache \
	depcomp install-sh libtool ltmain.sh missing
