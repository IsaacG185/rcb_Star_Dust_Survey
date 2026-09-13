#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Insert pipleline data in the photometry database
#
#  echo "source run_phot.csh >& ~/backup/2012_XX_XX/phot_XXX.log" | at now
#
#
#
echo "apass"
date
update_photometry -q apass -a  -l /dasch/Pipeline/photometry3.tmp -o /dasch/Pipeline/ingest/file_apass_update_photometry.log
date
update_limiting -q apass -l /dasch/Pipeline/ingest/file_apass_update_limiting.log -o /dasch/Pipeline/catalogs/tmpapasslimiting.dat
mv /dasch/Pipeline/catalogs/apasslimiting.dat /dasch/Pipeline/catalogs/oldapasslimiting.dat
mv /dasch/Pipeline/catalogs/apasslimiting.idx /dasch/Pipeline/catalogs/oldapasslimiting.idx
mv /dasch/Pipeline/catalogs/tmpapasslimiting.dat /dasch/Pipeline/catalogs/apasslimiting.dat
mv /dasch/Pipeline/catalogs/tmpapasslimiting.idx /dasch/Pipeline/catalogs/apasslimiting.idx
date
resort_magfiles -q apass  -o /dasch/Pipeline/ingest/resort_magfiles_apass.log 
date
search_none  -q apass -l 40.0 -c 1 -i   -o /dasch/Pipeline/ingest/search_none_apass_2013_03_30.log -d /dasch/Pipeline/ingest/search_none_apass_2013_03_30.db  -p /dasch/Pipeline/ingest/plate_outlier_apass_2013_03_30.db -u /dasch/Pipeline/ingest/id_unmatched_apass_2013_03_30.db -m /dasch/Pipeline/ingest/propermotion_2013_03_30.log
date
update_summary  -q apass  -n 2  -i -o /dasch/Pipeline/ingest/file_apass_update_summary.log
#
#
#
#
date
