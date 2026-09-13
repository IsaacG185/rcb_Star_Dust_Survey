#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Prepare the tables in /home/scanner/scanner/docs/office/aavso8/aavso8rev7.pdf
#
#  echo "source run_phot.csh >& ~/backup/2012_XX_XX/phot_XXX.log" | at now
#
#  NOTE: The RELEASE_LEVEL in pipelineutils.h should be set to the desired public release level
#
#  After: run idstats.c on all id tables.
#
echo "apass"
date
resort_magfiles -r -q apass  -o /dasch/Pipeline/ingest/resort_magfiles_apass_DR6A.log  -p /dasch/Pipeline/ingest/plates_apass_DR6A.txt -y /dasch/Pipeline/ingest/countvsyear_apass_DR6A.db -m /dasch/Pipeline/ingest/countvsmag_apass_DR6A.db -l /dasch/Pipeline/ingest/limvsyear_apass_DR6A.db
date
update_summary2  -q apass -v -l -90 -n 2 -i   -o /dasch/Pipeline/ingest/file_apass_update_summary2_2019_07_16.log  -u /dasch/Pipeline/ingest/id_apass_2019_07_16.db
date
#
#
#
echo "gsc2.3.2"
date
resort_magfiles -r  -o /dasch/Pipeline/ingest/resort_magfiles_gsc_DR6A.log -p /dasch/Pipeline/ingest/plates_gsc_DR6A.txt -y /dasch/Pipeline/ingest/countvsyear_gsc_DR6A.db -m /dasch/Pipeline/ingest/countvsmag_gsc_DR6A.db -l /dasch/Pipeline/ingest/limvsyear_gsc_DR6A.db
date
update_summary2  -v -l -90 -n 2 -i   -o /dasch/Pipeline/ingest/file_gsc_update_summary2_2019_07_16.log  -u /dasch/Pipeline/ingest/id_gsc_2019_07_16.db
date
#
#
#
echo "atlas"
date
resort_magfiles -r -q atlas  -o /dasch/Pipeline/ingest/resort_magfiles_atlas_DR6A.log -p /dasch/Pipeline/ingest/plates_atlas_DR6A.txt -y /dasch/Pipeline/ingest/countvsyear_atlas_DR6A.db -m /dasch/Pipeline/ingest/countvsmag_atlas_DR6A.db -l /dasch/Pipeline/ingest/limvsyear_atlas_DR6A.db
date
update_summary2  -q atlas -v -l -90 -n 2 -i   -o /dasch/Pipeline/ingest/file_atlas_update_summary2_2019_07_16.log  -u /dasch/Pipeline/ingest/id_atlas_2019_07_16.db
date
#
#
#
echo "kepler"
resort_magfiles -r -q kepler  -o /dasch/Pipeline/ingest/resort_magfiles_kepler_DR6A.log  -p /dasch/Pipeline/ingest/plates_kepler_DR6A.txt -p /dasch/Pipeline/ingest/plates_kepler_DR6A.txt -y /dasch/Pipeline/ingest/countvsyear_kepler_DR6A.db -m /dasch/Pipeline/ingest/countvsmag_kepler_DR6A.db -l /dasch/Pipeline/ingest/limvsyear_kepler_DR6A.db
update_summary2  -q kepler -v -l -90 -n 2 -i   -o /dasch/Pipeline/ingest/file_kepler_update_summary2_2019_07_16.log  -u /dasch/Pipeline/ingest/id_kepler_2019_07_16.db
date
#
#
#
date
