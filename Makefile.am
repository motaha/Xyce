# Process this file with "automake --foreign" to produce a Makefile.in
# NEVER edit the Makefile.in produced by automake or the Makefile produced
# by configure, as it will get overwritten when automake or configure
# are run

ACLOCAL_AMFLAGS=-I config

# define which subdirectories to do recursive makes in.
SUBDIRS=src $(MAYBE_PLUGIN)

if BUILD_XYCE_SHAREABLE

plugin:
	echo 'Making plugin in user_plugin'
	cd user_plugin && $(MAKE) $(AM_MAKEFLAGS) plugin

install-plugin:
	echo 'Making install-plugin in user_plugin'
	cd user_plugin && $(MAKE) $(AM_MAKEFLAGS) install-plugin
endif

# an ugly, ugly kludge
if SGI
clean-generic:
	rm -rf $(top_builddir)/ii_files $(top_builddir)/src/ii_files $(top_builddir)/src/*/ii_files $(top_builddir)/src/*/*/ii_files
endif

# another ugly, ugly kludge
if DEC
clean-generic:
	rm -rf $(top_builddir)/cxx_repository $(top_builddir)/src/cxx_repository $(top_builddir)/src/*/cxx_repository $(top_builddir)/src/*/*/cxx_repository
endif


# Yet another ugly, ugly kludge
if SUN
clean-generic:
	rm -rf $(top_builddir)/SunWS_cache $(top_builddir)/src/SunWS_cache $(top_builddir)/src/*/SunWS_cache $(top_builddir)/src/*/*/SunWS_cache
endif

dist-hook::
	rm -rf `find $(top_distdir) -type d -name .svn`
	rm -rf `find $(top_distdir) -type d -name .git`
	rm -rf `find $(top_distdir) -type d -name CVS`
