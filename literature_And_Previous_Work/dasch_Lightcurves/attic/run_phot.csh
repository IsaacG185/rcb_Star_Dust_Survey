#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Insert pipleline data in the photometry database
#
#  echo "source run_phot.csh >& ~/backup/2012_XX_XX/phot_XXX.log" | at now
#
#
set catalogs = ${DASCH_CATALOG:h}
#
echo "atlas insert"
date
update_photometry -a -q atlas  -l ${DASCH_SCRIPTS}/photometry.tmp -o ${DASCH_INGEST}/file_atlas_update_photometry.log
date
resort_magfiles -q atlas -o ${DASCH_INGEST}/resort_magfiles_atlas.log   -p ${DASCH_INGEST}/plates_atlas.txt  -g ${DASCH_INGEST}/binhistogram_atlas.txt 
date
update_limiting -q atlas -p -b 2 -l ${DASCH_INGEST}/file_atlas_update_limiting.log -o ${catalogs}/tmpatlaslimiting.dat
rm ${catalogs}/tmpatlaslimiting.tmp
cp ${catalogs}/atlaslimiting.dat ${catalogs}/oldatlaslimiting.dat 
cp ${catalogs}/atlaslimiting.idx ${catalogs}/oldatlaslimiting.idx
mv ${catalogs}/tmpatlaslimiting.dat ${catalogs}/atlaslimiting.dat
mv ${catalogs}/tmpatlaslimiting.idx ${catalogs}/atlaslimiting.idx
date
#
#
#
echo "atlas search_none  (use -t 2019_12_28 for transient flare searches) (use -g ${DASCH_SCRIPTS}/gsc_bin_repeat.list for repeat searches) (use -O to test the Optimized TC rematch algorithm)"
date
mpirun -np 31 search_none  -t 2019_12_28     -q atlas -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_atlas_2019_12_28.log -d ${DASCH_INGEST}/search_none_atlas_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_atlas_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_atlas_2019_12_28.db -m ${DASCH_INGEST}/propermotion_atlas_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_atlas_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_atlas.txt 
#date
#search_none   -t 2019_12_28    -q atlas -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_atlas_2019_12_28.log -d ${DASCH_INGEST}/search_none_atlas_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_atlas_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_atlas_2019_12_28.db -m ${DASCH_INGEST}/propermotion_atlas_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_atlas_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_atlas.txt 

#
#
#
echo "apass insert"
date
update_photometry  -q apass -a  -l ${DASCH_SCRIPTS}/photometry.tmp -o ${DASCH_INGEST}/file_apass_update_photometry.log
date
resort_magfiles -q apass  -o ${DASCH_INGEST}/resort_magfiles_apass.log  -p ${DASCH_INGEST}/plates_apass.txt -g ${DASCH_INGEST}/binhistogram_apass.txt 
date
update_limiting -p -b 2 -q apass -l ${DASCH_INGEST}/file_apass_update_limiting.log -o ${catalogs}/tmpapasslimiting.dat
rm ${catalogs}/tmpapasslimiting.tmp
cp ${catalogs}/apasslimiting.dat ${catalogs}/oldapasslimiting.dat
cp ${catalogs}/apasslimiting.idx ${catalogs}/oldapasslimiting.idx
mv ${catalogs}/tmpapasslimiting.dat ${catalogs}/apasslimiting.dat
mv ${catalogs}/tmpapasslimiting.idx ${catalogs}/apasslimiting.idx
#
#
#
#echo "apass search_none  (use -t 2019_12_28 for transient flare searches) (use -g ${DASCH_SCRIPTS}/gsc_bin_repeat.list for repeat searches) (use -O to test the Optimized TC rematch algorithm)"
#date
#mpirun -np 31 search_none    -t 2019_12_28   -q apass -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_apass_2019_12_28.log -d ${DASCH_INGEST}/search_none_apass_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_apass_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_apass_2019_12_28.db -m ${DASCH_INGEST}/propermotion_apass_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_apass_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_apass.txt 
#date
#search_none       -q apass -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_apass_2019_12_28.log -d ${DASCH_INGEST}/search_none_apass_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_apass_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_apass_2019_12_28.db -m ${DASCH_INGEST}/propermotion_apass_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_apass_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_apass.txt 
#
#
#
echo "gsc2.3.2 insert"
date
update_photometry -a   -l ${DASCH_SCRIPTS}/photometry.tmp -o ${DASCH_INGEST}/file_update_photometry.log
date
resort_magfiles  -o ${DASCH_INGEST}/resort_magfiles.log   -p ${DASCH_INGEST}/plates_gsc.txt  -g ${DASCH_INGEST}/binhistogram_gsc.txt 
date
update_limiting -p -b 2 -l ${DASCH_INGEST}/file_gsc_update_limiting.log -o ${catalogs}/tmpgsclimiting.dat
rm ${catalogs}/tmpgsclimiting.tmp
cp ${catalogs}/gsclimiting.dat ${catalogs}/oldgsclimiting.dat 
cp ${catalogs}/gsclimiting.idx ${catalogs}/oldgsclimiting.idx
mv ${catalogs}/tmpgsclimiting.dat ${catalogs}/gsclimiting.dat
mv ${catalogs}/tmpgsclimiting.idx ${catalogs}/gsclimiting.idx
#
#
#
#echo "gsc search_none  (use -t 2019_12_28 for transient flare searches) (use -g ${DASCH_SCRIPTS}/gsc_bin_repeat.list for repeat searches) (use -O to test the Optimized TC rematch algorithm)"
#date
#mpirun -np 31 search_none    -t 2019_12_28    -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_gsc_2019_12_28.log -d ${DASCH_INGEST}/search_none_gsc_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_gsc_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_gsc_2019_12_28.db -m ${DASCH_INGEST}/propermotion_gsc_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_gsc_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_gsc.txt 
#date
#search_none        -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_gsc_2019_12_28.log -d ${DASCH_INGEST}/search_none_gsc_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_gsc_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_gsc_2019_12_28.db -m ${DASCH_INGEST}/propermotion_gsc_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_gsc_2019_12_28.log -b ${DASCH_INGEST}/binhistogram_gsc.txt 
#
#
#
#echo "gaia insert"
#date
#update_photometry  -a -q gaia  -l ${DASCH_SCRIPTS}/photometry.tmp -o ${DASCH_INGEST}/file_gaia_update_photometry.log
#date
#resort_magfiles -q gaia -o ${DASCH_INGEST}/resort_magfiles_gaia.log   -p ${DASCH_INGEST}/plates_gaia.txt  -g ${DASCH_INGEST}/binhistogram_gaia.txt 
#date
#update_limiting -q gaia -p -b 4 -l ${DASCH_INGEST}/file_gaia_update_limiting.log -o ${catalogs}/tmpgaialimiting.dat
#rm ${catalogs}/tmpgaialimiting.tmp
#cp ${catalogs}/gaialimiting.dat ${catalogs}/oldgaialimiting.dat 
#cp ${catalogs}/gaialimiting.idx ${catalogs}/oldgaialimiting.idx
#mv ${catalogs}/tmpgaialimiting.dat ${catalogs}/gaialimiting.dat
#mv ${catalogs}/tmpgaialimiting.idx ${catalogs}/gaialimiting.idx
#
#
#
echo "atlas_update_summary (use -w for website statistitcs only)"
date
update_summary2  -w -q atlas -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_atlas_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_atlas_2019_12_28.db
#
#
#
echo "apass update_summary (use -w for website statistics only)"
date
update_summary2  -w  -q apass  -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_apass_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_apass_2019_12_28.db 
#
#
#
echo "gsc_update_summary (use -w for website statistitcs only)"
date
update_summary2  -w  -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_gsc_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_gsc_2019_12_28.db
#
#
#
#echo "gaia_update_summary (use -w for website statistitcs only)"
#date
#update_summary2  -w -q gaia -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_gaia_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_gaia_2019_12_28.db
#
#
#
#echo "gsc2.3.2 search_none"
#date
#mpirun -np 21 search_none  -v -l -90 -c 1 -i   -o ${DASCH_INGEST}/search_none_2019_12_28.log -d ${DASCH_INGEST}/search_none_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_gsc_2019_12_28.db -m ${DASCH_INGEST}/propermotion_gsc_2019_12_28.log -a ${DASCH_INGEST}/dradhistogram_gsc.log -b ${DASCH_INGEST}/binhistogram_gsc.txt
#date
#
#
#
#
#
#
#echo "kepler"
#date
#update_photometry -q kepler -a  -l ${DASCH_SCRIPTS}/photometry.tmp -o ${DASCH_INGEST}/file_kepler_update_photometry.log
#date
#resort_magfiles -q kepler  -o ${DASCH_INGEST}/resort_magfiles_kepler.log   -p ${DASCH_INGEST}/plates_kepler.txt 
#date
#update_limiting -p -b 16 -q kepler -l ${DASCH_INGEST}/file_kepler_update_limiting.log -o ${catalogs}/tmpkeplerlimiting.dat
#cp ${catalogs}/keplerlimiting.dat ${catalogs}/oldkeplerlimiting.dat 
#cp ${catalogs}/keplerlimiting.idx ${catalogs}/oldkeplerlimiting.idx 
#mv ${catalogs}/tmpkeplerlimiting.dat ${catalogs}/keplerlimiting.dat
#mv ${catalogs}/tmpkeplerlimiting.idx ${catalogs}/keplerlimiting.idx
#date
#search_none  -q kepler -c 1 -i -l -90.0  -o ${DASCH_INGEST}/search_none_kepler_2019_12_28.log -d ${DASCH_INGEST}/search_none_kepler_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_kepler_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_kepler_2019_12_28.db  -m ${DASCH_INGEST}/propermotion_kepler_2019_12_28.log  -a ${DASCH_INGEST}/dradhistogram_kepler.log 
#date
#update_summary2  -q kepler -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_kepler_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_kepler_2019_12_28.db
#
#  Currently only used for limiting magnitudes.
#
#echo "experimental"
#date
#combine_photometry -v -i  -l -q apass gsc2.3.2 experimental -o  ${DASCH_INGEST}/combine_photometry2019_12_28.log
#date
#update_photometry -q experimental -a  -l ${DASCH_SCRIPTS}/kepler.list -o ${DASCH_INGEST}/file_experimental_update_photometry.log
#date
#resort_magfiles -q experimental  -o ${DASCH_INGEST}/resort_magfiles_experimental.log   -p ${DASCH_INGEST}/plates_experimental.txt 
#date
#update_limiting -q experimental -l ${DASCH_INGEST}/file_experimental_update_limiting.log -o ${catalogs}/tmpexperimentallimiting.dat
#mv ${catalogs}/experimentallimiting.dat ${catalogs}/oldexperimentallimiting.dat 
#mv ${catalogs}/experimentallimiting.idx ${catalogs}/oldexperimentallimiting.idx 
#mv ${catalogs}/tmpexperimentallimiting.dat ${catalogs}/experimentallimiting.dat
#mv ${catalogs}/tmpexperimentallimiting.idx ${catalogs}/experimentallimiting.idx
#date
#search_none  -q experimental -c 1 -i   -o ${DASCH_INGEST}/search_none_experimental_2019_12_28.log -d ${DASCH_INGEST}/search_none_experimental_2019_12_28.db  -p ${DASCH_INGEST}/plate_outlier_experimental_2019_12_28.db -u ${DASCH_INGEST}/id_unmatched_experimental_2019_12_28.db  -a ${DASCH_INGEST}/dradhistogram_experimental.log 
#date
#update_summary2  -q kepler -v -l -90 -n 2 -i   -o ${DASCH_INGEST}/file_experimental_update_summary2_2019_12_28.log  -u ${DASCH_INGEST}/id_experimental_2019_12_28.db
#
#
echo "done"
date
#
#
#
#
#
