DASCH_PREFIX = /dasch
PREFIX = $(DASCH_PREFIX)
CC = gcc
CFLAGS = -Wall -g -O2 -Isrc -I$(DASCH_PREFIX)/include -DHAVE_CONFIG_H

ANET_KDTREE_LIBS = $(DASCH_PREFIX)/lib/libkd.a $(DASCH_PREFIX)/lib/libanutils.a $(DASCH_PREFIX)/lib/libqfits.a $(DASCH_PREFIX)/lib/libanbase.a

CFITSIO_LIBS = $(DASCH_PREFIX)/lib/libcfitsio.a -lz

# Expectation is that this will be provided by Giza:
CPGPLOT_LIBS = -L$(DASCH_PREFIX)/lib -lcpgplot

GD_CFLAGS = $(shell pkg-config --cflags gdlib)
GD_LIBS = $(shell pkg-config --libs gdlib)

MYSQL_CFLAGS = $(shell pkg-config --cflags mariadb)
MYSQL_LIBS = $(shell pkg-config --libs mariadb)

OPENMPI_CFLAGS = $(shell pkg-config --cflags ompi-c)
OPENMPI_LIBS = $(shell pkg-config --libs ompi-c)

PLPLOT_CFLAGS = $(shell pkg-config --cflags plplot)
PLPLOT_LIBS = $(shell pkg-config --libs plplot)

STARBASE_LIBS = $(DASCH_PREFIX)/lib/libtable.a $(DASCH_PREFIX)/lib/libutil.a

WCS_LIBS = $(DASCH_PREFIX)/lib/libwcs.a -lm

PROGRAMS = \
  bin/checksize \
  bin/comparefits \
  bin/contour_plot \
  bin/dasch-imwcs \
  bin/dasch_scat \
  bin/dasch-querycat \
  bin/dasch-querylc \
  bin/dasch-queryplates \
  bin/dasch-queryunmatched \
  bin/db2xyls \
  bin/filter_defect \
  bin/filter_multiple \
  bin/filter_wedge \
  bin/filterblended \
  bin/getdirectory \
  bin/getlocation \
  bin/getseries \
  bin/juliandate \
  bin/magdepcalibrate \
  bin/matchstars \
  bin/prepare_octave \
  bin/recover_points \
  bin/resort_magfiles \
  bin/restoretan \
  bin/search_close \
  bin/setFitWCS \
  bin/setScampRMS \
  bin/showflags \
  bin/update_limiting \
  bin/update_photometry \
  bin/update_sextractor \
  bin/update_summary2 \
  bin/votable \
  bin/web_authorize \
  bin/web_data \
  bin/web_edit \
  bin/web_plot \
  bin/web_point \
  bin/web_query

OCTAVES = \
  src/lowess.oct

all: $(PROGRAMS) $(OCTAVES)

install:
	mkdir -p $(PREFIX)/bin $(PREFIX)/pipeline
	cp -a bin/* $(PREFIX)/bin
	cp -a bin_scripts/* $(PREFIX)/bin
	cp -a scripts/* $(PREFIX)/pipeline
	cp -a src/*.oct $(PREFIX)/pipeline

clean:
	rm -f bin/* src/*.a src/*.o src/dasch-imwcs/*.o src/*.oct

src/libphotometry.a: src/galaxyutils.o src/kdtree.o src/parse_coo.o src/photdb.o src/photometryutils.o src/pipelineutils.o src/plottransient.o src/RA2fpix.o src/scandb.o
	rm -f $@ && ar cru $@ $^

bin/checksize: src/checksize.c
	$(CC) $(CFLAGS) -o $@ $^

bin/comparefits: src/comparefits.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS)

bin/contour_plot: src/contour_plot.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(CPGPLOT_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS) $(STARBASE_LIBS)

bin/dasch-imwcs: src/dasch-imwcs/imwcs.c src/dasch-imwcs/imsetwcs.o src/dasch-imwcs/matchstar.o src/dasch-imwcs/dasch_imwcs.h
	$(CC) $(CFLAGS) -o $@ $^ $(WCS_LIBS) $(CFITSIO_LIBS)

bin/dasch_scat: src/dasch_scat.c src/pipelineutils.o src/scandb.o src/scanread.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS) $(STARBASE_LIBS)

bin/dasch-querycat: src/dasch-querycat.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/dasch-querylc: src/dasch-querylc.c src/pipelineutils.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/dasch-queryplates: src/dasch-queryplates.c src/pipelineutils.o src/scandb.o src/scanread.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS) $(STARBASE_LIBS)

bin/dasch-queryunmatched: src/dasch-queryunmatched.c src/pipelineutils.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/db2xyls: src/db2xyls.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS) $(STARBASE_LIBS)

bin/filter_defect: src/filter_defect.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(CPGPLOT_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/filter_multiple: src/filter_multiple.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/filter_wedge: src/filter_wedge.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/filterblended: src/filterblended.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/getdirectory: src/getdirectory.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/getlocation: src/getlocation.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/getseries: src/getseries.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/juliandate: src/juliandate.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/magdepcalibrate: src/magdepcalibrate.c src/magdeputil.o src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(CPGPLOT_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/matchstars: src/matchstars.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/prepare_octave: src/prepare_octave.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/recover_points: src/recover_points.c src/magdeputil.o src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(CPGPLOT_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/resort_magfiles: src/resort_magfiles.c src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/restoretan: src/restoretan.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS)

bin/search_close: src/search_close.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/setFitWCS: src/setFitWCS.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/setScampRMS: src/setScampRMS.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(CFITSIO_LIBS)

bin/showflags: src/showflags.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/update_limiting: src/update_limiting.c src/libphotometry.a src/magdeputil.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/update_photometry: src/update_photometry.c src/libphotometry.a src/magdeputil.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/update_sextractor: src/update_sextractor.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(ANET_KDTREE_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/update_summary2: src/update_summary2.c src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/votable: src/votable.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS)

bin/web_authorize: src/web_authorize.c src/pipelineutils.o src/scandb.o
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS)

bin/web_data: src/web_data.c src/pipelineutils.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/web_edit: src/web_edit.c src/pipelineutils.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/web_plot: src/web_plot.c src/pipelineutils.o src/urlencode.o src/libphotometry.a
	$(CC) $(CFLAGS) $(GD_CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(GD_LIBS) $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/web_point: src/web_point.c src/pipelineutils.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

bin/web_query: src/web_query.c src/pipelineutils.o src/urlencode.o src/libphotometry.a
	$(CC) $(CFLAGS) $(MYSQL_CFLAGS) -o $@ $^ $(MYSQL_LIBS) $(WCS_LIBS) $(STARBASE_LIBS) $(PLPLOT_LIBS)

.c.o:
	$(CC) -c $(CFLAGS) $(MYSQL_CFLAGS) $(PLPLOT_CFLAGS) -o $@ $<

%.oct: %.cc
	mkoctfile -o $@ -D_GLIBCXX_ASSERTIONS $<
