#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Do photometric calibration and extract the results.
#
# Combined inputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww_background.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.db
#   ** for solution numbers $s$, 0 < s <= solnum **
# - if solnum > 0:
#   - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_allobjects.db
#     ** for solution numbers $s$, 0 <= s < solnum **
#
# Combined outputs. Files annotated with `(**)` may be missing
#
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01${rot}ww_s${sol_num + 1}${qual}_none.db (**)
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.grid (**)
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.out (**)
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.para (**)
# - $DASCH_CATALOGALL/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.db
# - $DASCH_CATALOGALL/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_extinction.db
# - $DASCH_CATALOGALL/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a0.db
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.db
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.colorterm.txt
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.db.ra.b
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.db.REF.i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.db.dec.i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.spatial_bins.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.REF-i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.dec-i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.ra-b
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_magdepcalibrate.db (**)
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_allobjects.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor_copy.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_u.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_b.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim1.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim2.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim3.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_drad.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_defect.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_summary_defect.db
#
# Database updates:
#
# - Potentially updates `FitWCS` in `scanner.mosaics` to set or clear `ColortermCrash` field
#
# Variables used:
#
# - Directories:
#   - $DASCH_ASTROMETRY
#   - $DASCH_BINOUTPUT
#   - $DASCH_CATALOGALL
#   - $DASCH_CATALOGBIN
#   - $DASCH_HEADERS
#   - $DASCH_INGEST
#   - $DASCH_MATCH
#   - $DASCH_SCRATCH
#   - $DASCH_SCRIPTS
# - Settings:
#   - $DASCH_CATALOG
#   - $DASCH_NUMBINS
#   - $DASCH_PLOT

set -xeuo pipefail

$DASCH_SCRIPTS/run_match_second.sh "$@"
$DASCH_SCRIPTS/filter_defect.sh "$@"
$DASCH_SCRIPTS/prepare_octave.sh "$@"
$DASCH_SCRIPTS/run_divide_annul9.sh "$@"
$DASCH_SCRIPTS/run_colorterm.sh "$@"
$DASCH_SCRIPTS/check_colorterm.sh "$@"
$DASCH_SCRIPTS/run_annular9.sh "$@"
$DASCH_SCRIPTS/ingest_matlab2.sh "$@"
$DASCH_SCRIPTS/run_local_calibration.sh "$@"
$DASCH_SCRIPTS/magdepcalibrate.sh "$@"
$DASCH_SCRIPTS/recover_points.sh "$@"
$DASCH_SCRIPTS/filter_multiple.sh "$@"

