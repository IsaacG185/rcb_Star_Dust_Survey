#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  August 29, 2012 Edward J. Los Initial version
#
#  Usage: echo "source run_id_table.csh >& run_id_table.log" | at now
#
#
# GSC
#
echo "CREATING id_table_2013_09_06.db"
cd /dasch/Pipeline/ingest
source /dasch/mysql/id_table.csh
mv id_table.db id_table_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_2013_09_06.db
mv summary_table.db summary_table_2013_09_06.db
date
row 'keplerField == "yes"' < id_table_2013_09_06.db > id_table_gsc_keplerfield_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_gsc_keplerfield_2013_09_06.db
mv summary_table.db summary_table_gsc_keplerfield_2013_09_06.db
#
# APASS
#
echo "CREATING id_table_apass_2013_09_06.db"
cd /dasch/Pipeline/ingest
source /dasch/mysql/id_table_apass.csh
mv id_table_apass.db id_table_apass_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_apass_2013_09_06.db
mv summary_table.db summary_table_apass_2013_09_06.db
row 'keplerField == "yes"' < id_table_apass_2013_09_06.db > id_table_apass_keplerfield_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_apass_keplerfield_2013_09_06.db
mv summary_table.db summary_table_apass_keplerfield_2013_09_06.db
#
# KIC
#
echo "CREATING id_table_kepler_2013_09_06.db"
cd /dasch/Pipeline/ingest
source /dasch/mysql/id_table_kepler.csh
mv id_table_kepler.db id_table_kepler_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_kepler_2013_09_06.db
mv summary_table.db summary_table_kepler_2013_09_06.db
row 'keplerField == "yes"' < id_table_kepler_2013_09_06.db > id_table_keplerfield_2013_09_06.db
source /dasch/mysql/summary_table.csh id_table_keplerfield_2013_09_06.db
mv summary_table.db summary_table_keplerfield_2013_09_06.db
#
# expermental
#
#echo "CREATING id_table_experimental_2013_09_06.db"
#cd /dasch/Pipeline/ingest
#source /dasch/mysql/id_table_experimental.csh
#mv id_table_experimental.db id_table_experimental_2013_09_06.db
#source /dasch/mysql/summary_table.csh id_table_experimental_2013_09_06.db
#mv summary_table.db summary_table_experimental_2013_09_06.db
#row 'keplerField == "yes"' < id_table_experimental_2013_09_06.db > id_table_experimental_keplerfield_2013_09_06.db
#source /dasch/mysql/summary_table.csh id_table_experimental_keplerfield_2013_09_06.db
#mv summary_table.db summary_table_experimental_keplerfield_2013_09_06.db
#
# Done
#
echo "COMPLETED: run_id_table.csh"
