#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Remove perl scripts and replace with prepare_mosaic and save_results
#  Mar 10, 2008 Edward J. Los - Remove save_results
#  Jul 16, 2008 Edward J. Los - Use catrms=0.02 for webda, sdss, and KeplerCam
#  Jan 13, 2009 Edward J. Los - Add colorterm.m  
#  Mar 17, 2009 Edward J. Los - Add update_photometry
#  Sep 28, 2009 Edward J. Los - Add save_header
#  Sep 29, 2009 Edward J. Los - Add multiple exposure support
#  Nov  3, 2009 Edward J. Los - Make find_astrometry2 a separate script
#  Feb 27, 2009 Edward J. Los - Add filter_multiple for pass2 scripts
#  Jul 23, 2012 Edward J. Los - Add solutionNumber to search_close.csh
#  Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
#  
#  Note: for apass solution 0 reprocessing, run search_close, run_sextractor_second, filter_wedge, run_match_second, and filter_defect
#                  solution n reprocessing, run_sextractor_second only is necessary
#
#  Note: use catrms=0.22 for gsc catalog
limit maxproc unlimited
set solutionNumber = $1
set imagelist = $2
@ nextSolution = $solutionNumber + 1
cd $DASCH_SCRIPTS
echo "directory: ${DASCH_SCRIPTS} imagelist: ${imagelist} solution: ${solutionNumber}"
date
#
#  Retrieve backup
#
foreach plate ( `cat $imagelist` ) 
    set matchsource = "/n/dasch14/Pipeline/match/${plate}*"
    echo "matchsource $matchsource"
    echo "cp $matchsource $DASCH_MATCH"
    cp $matchsource $DASCH_MATCH
end
#
# Copy to local node
#
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/copypipeline $imagelist -f
$DASCH_SCRIPTS/copypipeline $imagelist -f -q apass
#$DASCH_SCRIPTS/copypipeline $imagelist -f -q gaia
$DASCH_SCRIPTS/copypipeline $imagelist -f -q atlas
#$DASCH_SCRIPTS/copypipeline $imagelist -f -q experimental
#
#
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/compress_sextractor.csh $solutionNumber $imagelist
$DASCH_SCRIPTS/decompress_sextractor.csh $solutionNumber $imagelist
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_scratch.csh $solutionNumber $imagelist
#
#
# Main Pipeline
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_mosaic.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_sextractor.csh $solutionNumber $imagelist 
$DASCH_SCRIPTS/run_sextractor_second.csh $solutionNumber $imagelist 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/search_close.csh $solutionNumber $imagelist 
#$DASCH_SCRIPTS/search_close3.csh $solutionNumber $imagelist 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_scamp.csh $solutionNumber $imagelist 
#cd $DASCH_SCRIPTS 
#$DASCH_SCRIPTS/run_sextractor_second.csh $solutionNumber $imagelist 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/find_astrometry01.csh $solutionNumber $imagelist 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_scamp_01.csh $solutionNumber $imagelist 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_scamp.csh $solutionNumber $imagelist 2 
#cd $DASCH_SCRIPTS 
#$DASCH_SCRIPTS/run_sextractor_second.csh $solutionNumber $imagelist  2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/save_header.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_wedge.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist
#
##################### END OF VERSION INDEPENDENT CODE ###############
#
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/check_colorterm.csh  $solutionNumber $imagelist
cd $DASCH_SCRIPTS
run_annular9.csh $solutionNumber $imagelist 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.22 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/find_astrometry2.csh $nextSolution $imagelist
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist
#
#####################################################################
#
#  PASS TWO
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  pass2  id_all_M44_total.db
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist pass2 
#cd $DASCH_SCRIPTS
#run_annular9.csh $solutionNumber $imagelist pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.22 pass2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist pass2
#
#  webda2
#
#setenv DASCH_NUMBINS     "1"
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  webda2  webda.db
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist webda2
#cd $DASCH_SCRIPTS
#run_annular9.csh $solutionNumber $imagelist webda2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist webda2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist webda2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist webda2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 webda2
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist webda2
#
#  kepler
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  kepler 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#run_annular9.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist kepler
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist kepler
#
#  experimental
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  experimental 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#run_annular9.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist experimental
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist experimental
#
#  apass  NOTE: prerequisites: search_close, run_sextractor_second, and filter_defect
#
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  apass 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
run_annular9.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist apass
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist apass
#
#  gaia  NOTE: prerequisites: search_close, run_sextractor_second, and filter_defect
#
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  gaia 
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#run_annular9.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist gaia
#cd $DASCH_SCRIPTS
#$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist gaia
#
#  atlas  NOTE: prerequisites: search_close, run_sextractor_second, and filter_defect
#
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_match_second.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_defect.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/prepare_octave.csh $solutionNumber $imagelist  atlas 
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_divide_annul9.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_colorterm.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
run_annular9.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/ingest_matlab2.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/run_local_calibration.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/magdepcalibrate.csh $solutionNumber $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/recover_points.csh $solutionNumber $imagelist catrms=0.02 atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/filter_multiple.csh $nextSolution $imagelist atlas
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/clean_intermediate.csh $solutionNumber $imagelist atlas
#
# Copy back operation
#
cd $DASCH_SCRIPTS
$DASCH_SCRIPTS/copypipeline $imagelist  -r
$DASCH_SCRIPTS/copypipeline $imagelist -r -q apass
#$DASCH_SCRIPTS/copypipeline $imagelist -r -q gaia
$DASCH_SCRIPTS/copypipeline $imagelist -r -q atlas
#$DASCH_SCRIPTS/copypipeline $imagelist -r -q kepler
#$DASCH_SCRIPTS/copypipeline $imagelist -r -q experimental
#
#
