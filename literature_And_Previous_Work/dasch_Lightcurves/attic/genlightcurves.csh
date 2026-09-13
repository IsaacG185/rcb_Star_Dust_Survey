#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# 
#  
rm /dasch/Pipeline/ingest/master_all_M44_total.db
rm /dasch/raid011/Pipeline/ingest/master_all_M44_total.db
cd /dasch/Pipeline
./extract_goodpoints.csh M44_total.list
mv /dasch/Pipeline/ingest/master_all_M44_total.db /dasch/raid011/Pipeline/ingest/master_all_M44_total.db
cd /dasch/Pipeline/ingest
ln -s /dasch/raid011/Pipeline/ingest/master_all_M44_total.db master_all_M44_total.db
cd /dasch/Pipeline
run_lc.csh 'M44'  M44_total.list

rm /dasch/Pipeline/ingest/master_all_kepler.db
rm /dasch/raid011/Pipeline/ingest/master_all_kepler.db
cd /dasch/Pipeline
./extract_goodpoints.csh kepler.list
mv /dasch/Pipeline/ingest/master_all_kepler.db /dasch/raid011/Pipeline/ingest/master_all_kepler.db
cd /dasch/Pipeline/ingest
ln -s /dasch/raid011/Pipeline/ingest/master_all_kepler.db master_all_kepler.db

rm /dasch/Pipeline/ingest/master_all_other.db
rm /dasch/raid011/Pipeline/ingest/master_all_other.db
cd /dasch/Pipeline
./extract_goodpoints.csh other.list
mv /dasch/Pipeline/ingest/master_all_other.db /dasch/raid011/Pipeline/ingest/master_all_other.db
cd /dasch/Pipeline/ingest
ln -s /dasch/raid011/Pipeline/ingest/master_all_other.db master_all_other.db


#rm /dasch/Pipeline/ingest/master_all_3c273_total.db
#rm /dasch/raid011/Pipeline/ingest/master_all_3c273_total.db
#cd /dasch/Pipeline
#./extract_goodpoints.csh 3c273_total.list
#mv /dasch/Pipeline/ingest/master_all_3c273_total.db /dasch/raid011/Pipeline/ingest/master_all_3c273_total.db
#cd /dasch/Pipeline/ingest
#ln -s /dasch/raid011/Pipeline/ingest/master_all_3c273_total.db master_all_3c273_total.db
#
#rm /dasch/Pipeline/ingest/master_all_baade.db
#rm /dasch/raid011/Pipeline/ingest/master_all_baade.db
#cd /dasch/Pipeline
#./extract_goodpoints.csh baade.list
#mv /dasch/Pipeline/ingest/master_all_baade.db /dasch/raid011/Pipeline/ingest/master_all_baade.db
#cd /dasch/Pipeline/ingest
#ln -s /dasch/raid011/Pipeline/ingest/master_all_baade.db master_all_baade.db
#
#
