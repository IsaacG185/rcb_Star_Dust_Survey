# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  evaluate_wedge.csh - Evaluate Pickering Wedge results
#
#  May 29, 2008 Edward J. Los
#  Jun 12, 2008 Edward J. Los - remove plotting script
#  Jan 27, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS, remove CAL_FLAG
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS
set binning = 16
set plotstr = ""
set useRaDec = 0
date
cd ${DASCH_MATCH}
set listroot = $imagelist:t

pwd
echo "Processing ${DASCH_SCRIPTS}/$listroot"
foreach plate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 
   set wedgename = ${DASCH_MATCH}/${plate}_wedge.db
   set allobjectsname = ${DASCH_INGEST}/${plate}_allobjects.db
   if ((-e $wedgename) && (-e $allobjectsname)) then
      echo "Processing plate $plate"
      set outputname = ${DASCH_MATCH}/${plate}_evaluatewedge.db
      set tempname  = ${DASCH_MATCH}/${plate}_evaluatewedge.tmp
      set tempname1 = ${DASCH_MATCH}/${plate}_evaluatewedge1.tmp
      set tempname2 = ${DASCH_MATCH}/${plate}_evaluatewedge2.tmp
      set tempname3 = ${DASCH_MATCH}/${plate}_evaluatewedge3.tmp
      set tempname4 = ${DASCH_MATCH}/${plate}_evaluatewedge4.tmp
      set tempname5 = ${DASCH_MATCH}/${plate}_evaluatewedge5.tmp
      set tempname6 = ${DASCH_MATCH}/${plate}_evaluatewedge6.tmp
      set tempname7 = ${DASCH_MATCH}/${plate}_evaluatewedge7.tmp
      set tempname8 = ${DASCH_MATCH}/${plate}_evaluatewedge8.tmp
     
      cp $wedgename $tempname
      cp $allobjectsname $tempname1
#     echo "Get the magnitudes of the dim objects"
      row 'MATCH_NUMBER != 0' < $tempname | column NUMBER MATCH_NUMBER > $tempname2
      column -i $tempname1 NUMBER magcal_local magcal_local_rms AFLAGS BFLAGS > $tempname3
      index -mb -n $tempname2 NUMBER
      index -mb -n $tempname3 NUMBER
      jointable -j NUMBER -n $tempname2 $tempname3 > $tempname4
#     echo "Get the magnitudes of the bright objects"
      row 'MATCH_NUMBER != 0' < $tempname | column MATCH_NUMBER | column -a NUMBER | compute 'NUMBER = MATCH_NUMBER' | column NUMBER > $tempname5
      index -mb -n $tempname5 NUMBER
      index -mb -n $tempname3 NUMBER
      jointable -j NUMBER -n $tempname5 $tempname3 > $tempname6
#     echo "Now join the two together"
      column -a -i $tempname6 MATCH_NUMBER | compute 'MATCH_NUMBER = NUMBER' | column MATCH_NUMBER magcal_local magcal_local_rms AFLAGS BFLAGS > $tempname7
      index -mb -n $tempname4 MATCH_NUMBER
      index -mb -n $tempname7 MATCH_NUMBER
      jointable -j MATCH_NUMBER -n $tempname4 $tempname7 > $tempname8
      column -a -i  $tempname8 magcal_diff | compute 'magcal_diff = magcal_local_2 - magcal_local_1' > $outputname
      
      rm $tempname
      rm $tempname1
      rm $tempname2
      rm $tempname3
      rm $tempname4
      rm $tempname5
      rm $tempname6
      rm $tempname7
      rm $tempname8

      rm ${tempname2}.NUMBER.b
      rm ${tempname3}.NUMBER.b
      rm ${tempname5}.NUMBER.b

      rm ${tempname4}.MATCH_NUMBER.b
      rm ${tempname7}.MATCH_NUMBER.b


   endif
end
date

