#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Do all of the steps to prepare `ww` mosaics for photometry.
#
# Combined inputs:
#
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww.fit
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_16${rot}ww.fit
# - if solnum > 0:
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${prev_sol_tag}_tnx.db
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}.hdr
#   - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01${rot}ww_s${solnum}_none.db
#
# Combined outputs:
#
# - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.hdr
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.db
# - if solnum = 0:
#   - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww_tnx.fit
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_fgroups_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_referr1d_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_referr2d_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_distort_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_drad_map.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_drad_offset.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww.refcat.ucac5.fits.reg
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_scamp.sexcat_filtered.fits.reg
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_scamp.sexcat_filtered.fits.xy.reg
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_16${rot}ww.hdr
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww.db
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww_background.db
#
# Database updates:
#
# - None, really?
#
# Variables used:
#
# - $DASCH_ASTROMETRY
# - $DASCH_HEADERS
# - $DASCH_MATCH
# - $DASCH_SCAMPIMAGES
# - $DASCH_SCRATCH
# - $DASCH_SCRIPTS

set -xeuo pipefail

$DASCH_SCRIPTS/run_sextractor.sh "$@"
$DASCH_SCRIPTS/search_close.sh "$@"
$DASCH_SCRIPTS/run_scamp.sh "$@"
$DASCH_SCRIPTS/run_sextractor_second.sh "$@"
$DASCH_SCRIPTS/find_astrometry01.sh "$@"
$DASCH_SCRIPTS/run_scamp_01.sh "$@"
$DASCH_SCRIPTS/run_scamp.sh "$@" 2 # <== note "2" mode
$DASCH_SCRIPTS/run_sextractor_second.sh "$@" 2 # <== note "2" mode
$DASCH_SCRIPTS/save_header.sh "$@"
$DASCH_SCRIPTS/filter_wedge.sh "$@"
