#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Prepare a `tnx` mosaic for photometry. This does some of the same steps as
# `ww_prep_stage.sh`, but under the assumption that we have a good TNX that we
# do not want to alter. This is different than `recreate_ww`, which also starts
# with a TNX but sets up to rerun the distortion (TNX) calibration and generate
# a new TNX file.
#
# Combined inputs:
#
# - if solnum = 0:
#   - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww_tnx.fit
# - else:
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${prev_sol_tag}_tnx.db
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}_tnx.hdr
#
# Combined outputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.db
# - if solnum = 0:
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww_background.db
#
# Database updates:
#
# - None
#
# Variables used:
#
# - $DASCH_ASTROMETRY
# - $DASCH_HEADERS
# - $DASCH_MATCH
# - $DASCH_SCRATCH
# - $DASCH_SCRIPTS

set -xeuo pipefail

$DASCH_SCRIPTS/run_sextractor_second.sh "$@"
$DASCH_SCRIPTS/search_close.sh "$@"
$DASCH_SCRIPTS/run_sextractor_second.sh "$@" 2 # <== note "2" mode
