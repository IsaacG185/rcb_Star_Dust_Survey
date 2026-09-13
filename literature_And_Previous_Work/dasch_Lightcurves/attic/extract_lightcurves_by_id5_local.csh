# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# extract_lightcurves script by Silas Laycock, June 2007
# Based on earlier efforts. 
# Simplified Lightcurve extractor. no fiducial image, just a search catalog and a set of catalogs to search
# All fancy cross-calibration stuff removed.
# The target coordinates are specified in a "key-file" RA,Dec
# (1)  Loop over list of plate catalogs, searching for all lightcurve-stars by indexed star ID.
# (3)  Write lightcurves
#
# Highly optimized version that uses the minimum number of "search" calls.
# Speed is roughly 1 second per lightcurve, which is 200X faster than previous version!
# June 4, 2007: Introduce upper limits for non detected stars that are on a given plate.
# This is done by pre-filtering the target list to identify all targets that lie on the plate.
# Then only this list is searched, and if the star is not found, its magnitude is set to the limiting-magnitude
# in such cases the limiting magnitude value is then set to 99 and its error to 0.
# June 5, 2007: Determine the correct spatial_bin for undetected stars and hence apply the correct limit.
# June 13, 2007: Correct the spatial bin definitions to account for the outermost bin being square.
#              : Spatial bin definitions files have been upgraded so that the correct selections are automatically made.
# July 19 2007: Check for existence of "magcal_local" and process if found. 
#               More flexible identification of spatial_bin definition file PLATE.*.spatial_bins.db
# Aug 8, 2007:  Small modification so that magcal_local=limiting_mag when a star is not detected, previously it was flagged as 0.
#		Added special options to enable the final index command to handle very large master_target_data.db ( setenv TEMPDIR = . )
# Sept 12, 2007: Silas Laycock. Scrutinized flag/placeholder settings. Introduced "detected" flag-column. Where "detected=0" denoes
#                              stars for which we only have an upper limit on the magnitude, this removes ambiguity concerning stars 
#	                       that are detected with magnitude > limiting_mag from those which were not detected at all.
#				AS A RESULT THE PLOTTING/FLAGGING SCRIPT MUST BE UPDATED to read detected=0 (upper-lim), detected=1(measured magnitude)
#
# Sep 20, 2007: Edward J. Los  Merge with extract_lightcurves_by_id_local7.csh and Change TEMPDIR to TMPDIR
# Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
# Nov 14, 2007 Edward J. Los - Eliminate JulianDates.txt and access the MySQL database directly
# Feb 19, 2008 Edward J. Los - Correct Julian Date

set scripts = $DASCH_SCRIPTS

# Here are the command line arguments for the script
set catalog_list = $1        
set target_catalog = $2  
set catalogpath = $3
setenv TMPDIR .

# Error trap to make sure all arguments are specified
if ( $#argv != 3 ) then
  echo "usage: 1. catalog_list 2. target_catalog 3. catalog_path 4. Date_list"
  echo "### 1. list of catalogs to process. Two column format:  path/filename, DATE"
  echo "### 2. Table (starbase format) of objects for which we want to extract lightcurves."
  echo "### 3. Catalog path (directory path to the photometry catalogs, relative or absolute))" 
  exit 
else 
  # report the time and date at start
  echo Begin processing `date`
endif


# Typical Columns appearing in post-MATLAB catalogs
# REF ra dec stdmag color MAG_ISO MAG_APER ISOAREA_IMAGE magcal_aper magcal_aper_err_a magcal_aper_err_b magcal_iso magcal_iso_err_a magcal_iso_err_b magcal_area 
# magcal_area_err_a magcal_area_err_b

column -i $target_catalog src_name REF ra dec >! target_catalog.tmp
index -mb -n target_catalog.tmp ra
index -mi -n target_catalog.tmp dec
index -mi target_catalog.tmp REF

set new = 1
foreach catalog ( `cat $catalog_list | gawk '{print $1}'` )
   set plate = `echo ${catalog:t:r} | gawk -F"_" '{print $1}'`
#   set date = `grep $plate $date_list | gawk '{print $2}'`
   set date = `$DASCH_SCRIPTS/juliandate $plate | gawk '{print $1}'`
   if ( $#date == 0 ) then 
      echo "No date found for" $plate in input list 
      set date = 0 
   endif
   echo Date $date $plate
   set catalog = ${catalogpath}/${catalog}.out.local.db
   if (-e $catalog) then
     # Measure the boundaries of the sky-region covered by the plate/catalog.
     column -i $catalog ra dec | mintable | column -b >! plate.limits 
     column -i $catalog ra dec | maxtable | column -b >> plate.limits 
     set ra_min = `gawk 'NR==1{print $1}' < plate.limits`
     set ra_max = `gawk 'NR==2{print $1}' < plate.limits`
     set dec_min = `gawk 'NR==1{print $2}' < plate.limits`
     set dec_max = `gawk 'NR==2{print $2}' < plate.limits`
     set ra_center = `echo $ra_max $ra_min | gawk '{print ($1+$2)/2}'`
     set dec_center = `echo $dec_max $dec_min | gawk '{print ($1+$2)/2}'`
     echo ra_min $ra_min ra_max $ra_max dec_min $dec_min dec_max $dec_max ra_center $ra_center dec_center $dec_center
  
     # Select stars from the target list whose positions lie on the plate. Using the RA/Dec limits measured above. 
     # This saves a lot of time because there is no point trying to find stars that are not obviously not in the field of view.
     if ( ! -e target_catalog_${plate}.tmp ) then
        row -i target_catalog.tmp "ra > ${ra_min} && ra < ${ra_max} && dec > ${dec_min} && dec < ${dec_max}" | column REF ra dec | sorttable REF >! target_catalog_${plate}.tmp
     else
        echo using existing target_catalog_${plate}.tmp   
     endif
  
     # Search for the selected stars, using the fast, indexed starbase search command. Please make sure all input catalogs are indexed already.
     # In case this is a reprocessing the search result may already exist.
     # Note: setting of the "detected=1" flag for these stars, to later distinguish non detected stars (detected=0) unambiguously.
    if ( ! -e target_data_${plate}.tmp ) then
         if ( ! -e ${catalog}.REF-i ) then
            echo "YOUR CATALOG" $catalog "IS NOT INDEXED by REF COLUMN!"
            echo "Indexing..... (this will take a while)......."
            index -mi $catalog REF
         endif
         search $catalog REF < target_catalog_${plate}.tmp | sorttable REF | column -a detected | compute 'detected=1' >! target_data_${plate}.tmp     
     else
         echo using existing target_data_${plate}.tmp
     endif
  
     # Find out what columns are in the table, and use this "columns" variable to propagate them throughout the script 
     # check for the existence of "magcal_local" because it will need to be used if it is present.
     # The flag-variable "LOCAL" provides an automatic switch to handle data that have been locally recalibrated as well as data that have not. 
     set columns = `head -1 target_data_${plate}.tmp`
     echo "Columns =" $columns	
     set LOCAL = `echo $columns | grep -c "magcal_local"`		
     echo "LOCAL =" $LOCAL
  
     ###########################################
     ### NON-DETECTED STARS AND UPPER LIMITS ###
     ###########################################
     # This means objects in the target list (stars for which you are extracting lightcurves) that were not found on this plate despite being in the field of view.
     # These are important because we can give upper limits based on the local limiting magnitude. For example ecliping binary in mid-eclipse, dwarf nova in quiescence.
     # Retrieve the non-detected stars and calculate how far from the plate center they would have been, and assign limits, placeholders and flags accordingly.
     # Note that at this stage we propagate a lot of columns using the "columns" parameter
     # Note the setting of the "detected=0" flag   
  
     # Retrieve the non-detected stars:
     if ( ! -e target_data_${plate}.2.tmp ) then
       jointable -v1 -j REF target_catalog_${plate}.tmp target_data_${plate}.tmp | sed 's/ra_1/ra/g' | sed 's/dec_1/dec/g' | compute "plate_dist=sqrt( (((ra - $ra_center)*(cos(((dec + $dec_center)/2)/57.29577951))))^2 + (dec - $dec_center)^2); detected=0" | column $columns >! target_data_${plate}.2.tmp
     else
       echo using existing target_data_${plate}.2.tmp  
     endif
  
     # Assign Upper-limits, Flags and Placeholders.
     # This process uses a small fortran program (assign_spatial_bins_2) to assign the correct spatial bin and limiting magnitude to each undetected star.
     # There should be a lookup table for each plate, "PLATE.*.spatial_bins.db" listing the spatial bins used for calibration. 
     # The bin the boundaries in RA/Dec/plate_dist, limiting_mag and matlab fit quality are listed for each bin.
     # These lookup tables were created by "ingest_matlab.csh" and will be used to assign magnitude limits to the undetected stars.
     # if the lookup table is missing, the code below will report the problem and generate error flags in place of magnitudes.
  
       echo processing upper-limits for non-detected stars
       cp target_data_${plate}.tmp target_data_${plate}.all.tmp
       set spatial_bins = `ls ${catalogpath}/${plate}*.spatial_bins.db`
       echo $spatial_bins
       if ( $#spatial_bins > 0 ) then
          echo "Found" $spatial_bins
          column -i $spatial_bins plate_dist_min plate_dist_max ra_min ra_max dec_min dec_max spatial_bin limiting_mag >! limits.tmp
          column -i target_data_${plate}.2.tmp plate_dist ra dec >! x.tmp
          echo x.tmp >! assign_spatial_bin.par
          echo limits.tmp >> assign_spatial_bin.par
          ${scripts}/assign_spatial_bins_2 < assign_spatial_bin.par | gawk '{OFS="\t"}{print $1,$2}' >! values.tmp
          if ( $LOCAL == 1 ) then 
             # These objects were not detected so their magnitudes are set to the local limiting magnitude.
             paste target_data_${plate}.2.tmp values.tmp | compute 'magcal_iso=MAGLIMITVAL; magcal_iso_error=0; spatial_bin=SPATIALBINVAL;limiting_mag=MAGLIMITVAL; magcal_local=MAGLIMITVAL; magcal_local_error=0; npoints_local=0; cal_local=0' | column -b $columns >> target_data_${plate}.all.tmp 	
  	endif
          if ( $LOCAL == 0 ) then
             paste target_data_${plate}.2.tmp values.tmp | compute 'magcal_iso=MAGLIMITVAL; magcal_iso_error=0; spatial_bin=SPATIALBINVAL;limiting_mag=MAGLIMITVAL' | column -b $columns >> target_data_${plate}.all.tmp 	
          endif
       else 
          # This block deals with missing plate-limits file (produced by ingest_matlab.csh) by inserting place-holder values. A potential trouble spot.
  	echo "missing plate-limits file for" $plate
  	if ( $LOCAL == 1 ) then
  	   compute -i target_catalog_${plate}.2.tmp 'magcal_iso=99; magcal_iso_error=0; spatial_bin=99; limiting_mag=0; magcal_local=99; magcal_local_error=0; npoints_local=0; cal_local=0' | column -b $columns >> target_catalog_${plate}.all.tmp
          endif
  	if ( $LOCAL == 0 ) then   
             compute -i target_catalog_${plate}.2.tmp 'magcal_iso=99; magcal_iso_error=0; limiting_mag=0; spatial_bin=99' | column -b $columns >> target_catalog_${plate}.all.tmp 
  	endif
       endif
  
    # Now Finally write all the data extracted from this plate to the master table:
   
     if ( $new == 1 ) then 
         cat target_data_${plate}.all.tmp | setcolumn Plate ${plate} | setcolumn Date ${date} >! master_target_data.db
         set new = 0
     else
         cat target_data_${plate}.all.tmp | setcolumn Plate ${plate} | setcolumn Date ${date} | column -b >> master_target_data.db
     endif

     rm target_catalog_${plate}.tmp
     rm target_data_${plate}.tmp
     rm target_data_${plate}.2.tmp
     rm target_data_${plate}.all.tmp

   else
     echo "Can not find $catalog"
   endif
end

echo "Catalog search completed."
date 
echo "Writing out lightcurves......"

# Now search the master table for each object:
index -M 10 -mb master_target_data.db REF
set l = `column -i target_catalog.tmp | wc -l` 
set n = 3
while ( $n <= $l )
   echo "REF" >! key.tmp
   echo "---" >> key.tmp
   set object = `sed -n $n"p" target_catalog.tmp | gawk '{print $2}'`
   set src_name = `sed -n $n"p" target_catalog.tmp | gawk '{print $1}'`   
   echo $object >> key.tmp
   search master_target_data.db REF < key.tmp >! lc_${src_name}.db 
  @ n = $n + 1
end

rm key.tmp
rm assign_spatial_bin.par
rm limits.tmp
rm plate.limits
rm values.tmp
rm target_catalog.tmp
rm target_catalog.tmp.dec.i
rm target_catalog.tmp.ra.b
rm target_catalog.tmp.REF-i
#rm master_target_data.db
rm master_target_data.db.REF-b 
rm x.tmp
echo "done"
date
