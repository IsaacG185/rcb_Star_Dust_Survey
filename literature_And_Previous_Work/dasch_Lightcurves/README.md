# DASCH Pipeline

This repository contains the tools and scripts underlying the legacy DASCH
scientific processing pipeline.


## Direct Linux build (not recommended)

Build with a simple `make` ... once you install all of the dependencies:

- The [Astrometry.Net](http://astrometry.net/doc/readme.html) local tools
- [Astropy](https://www.astropy.org/)
- [Funtools](https://github.com/ericmandel/funtools/) (this is barely used)
- [libgd](https://libgd.github.io/) (only for `web_plot`)
- [Giza](https://danieljprice.github.io/giza/) (PGPLOT-compatible library)
- ImageMagick (required by SCAMP)
- A **highly customized version of** [wcstools](http://tdc-www.harvard.edu/wcstools/)
  - This is provided by the `daschapps` Docker image described in the next
    section.
- MariaDB/MySQL
- Matplotlib
- [Octave](https://octave.org/)
- [PLplot](https://plplot.sourceforge.net/)
- [SCAMP](https://www.astromatic.net/software/scamp/)
- [SExtractor](https://www.astromatic.net/software/sextractor/)
- [Starbase](http://hopper.si.edu/wiki/mmti/Starbase)


## Containerized build

You can build a Dockerized version of the DASCH pipeline like so:

```
docker build . -t daschpipe:latest
```

This process derives from the `daschapps:latest` image, built from the DASCH
"apps" repository.

Then, generate an Apptainer/Singularity image from the Docker image:

```
docker run \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):/output:rw,z \
  --privileged \
  -t --rm singularityware/docker2singularity -n daschpipe daschpipe:latest
```

This will generate a file named `daschpipe.simg`.


## Brief notes

The pipeline steps are now orchestrated in `daschdata/pipeline/` and the
`daschops` repository. It might be helpful to merge all of these trees into a
monorepo to help coordinate updates.

See `scripts/ww_prep_stage.sh` and `scripts/photometry_stage.sh` for some chunks
of steps, but these do not capture the full breadth of the pipeline processing.

The original forensics about core pipeline operations were based on
`runpipeline.c` and the `full*.csh` files.


## Environment variables

The pipeline code cares about the following environment variables:

- Database:
  - `DASCH_MYSQLHOST` - hostname of the database server (use 127.0.0.1 instead
    of `localhost` for Reasons)
  - `DASCH_USERNAME` - username for accessing the database server
  - `DASCH_PASSWORD` - password for accessing the database server
  - `DASCH_MYSQLUSERHOST` - set this to `${DASCH_USERNAME}@${DASCH_MYSQLHOST}`
  - `MYSQL_TCP_PORT` - override MySQL connection port from the default of 3306
  - `DASCH_PHOT_MYSQLHOST` - host for photometry database server
  - `DASCH_PHOT_USERNAME` - username for photometry database server
  - `DASCH_PHOT_PASSWORD` - password for photometry database server
- `DASCH_ASTROMETRY` - directory for astrometry solution data
- `DASCH_BINOUTPUT` - directory for annular-bin analysis data
- `DASCH_CATALOG` - internal catalog file path, e.g. `/n/dasch12/catalogs/gsc232bin.dat`
  - This will be mutated in a bunch of ways by various scripts
  - No `.` may appear in the path except at the very end (the `.dat`)
- `DASCH_CATALOGALL` - directory for Octave files
- `DASCH_CATALOGBIN` - a different directory for Octave files
- `DASCH_GID` - the name of the DASCH Unix group (`dasch_project`)
- `DASCH_HEADERS` - directory for saving headers
- `DASCH_INGEST` - directory for ingested Octave data
- `DASCH_MATCH` - source-matching data directory; historical locations:
  - `boslfs:dasch14/Pipeline/match`
  - tiny backup in `boslfs:dasch4_backup/dasch/raid014/Pipeline/match`
  - and on `dasch4` corresponding to the above
- `DASCH_NUMBINS` - set to `9` (used by Octave scripts)
- `DASCH_PHOT_ROOT` - root directory of photometry database
  - Internally, `DASCH_PHOT_MAGNITUDES@KEY@` specify individual refcat data directories
- `DASCH_PLOT` - set to `YES` or `NO` (used by Octave scripts)
- `DASCH_RAID_OVERRIDE` - if present, will cause the `getdirectory` program to
  use its value as a prefix for file locations, rather than `/dasch/raidNNN` based on the
  database. This is useful for compartmentalized local testing.
- `DASCH_SCAMPIMAGES` - a directory to store images generated during SCAMP processing
- `DASCH_SCRATCH` - a scratch directory
- `DASCH_SCRIPTS` - pipeline scripts location; should be `/dasch/pipeline` in
  the containerized install
- `DASCH_WCSFIT_DATA` - directory for storing data for the initial WCS fitting


## History

This codebase was originally imported into a Git repository by a reconstruction
of various sets of files strewn throughout the DASCH storage systems.
Modification history prior to 2023 is hinted at by changelog messages embedded
in some of the files, and captured coarsely in a set of Git commits
reconstructing historical checkpoints. Use version control, folks.

The `attic` subdirectory contains old code that was not needed to complete
DASCH scanning or Data Release 7 in the 2023-2024 timeframe. Everything else
played some kind of role in these activities.


## Availability

An export of this Git repository was archived as part of DASCH Data Release 7.
Starting in 2023 the codebase was managed in a private GitLab repository at this
URL:

https://gitlab.com/HarvardRC/rse/cfa-dasch/legacy_code/legacy_pipeline

The President and Fellows of Harvard College are the copyright owners of all
original DASCH work. Everything in this repository is licensed under the
permissive MIT license unless specified otherwise.


## Authorship

Prior to 2023, this codebase was primarily the work of Edward J. Los, with
contributions from Sumin Tang and Mathieu Servillat as well as third-party
sources. Work done in the 2023-2024 period was done by Peter K. G. Williams.
