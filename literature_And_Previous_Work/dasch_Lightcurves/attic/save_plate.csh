# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Save the intermediate files of a photometry run for further analysis
# 
# Dec 12, 2014 Edward J. Los 

set plate = $1

if ($#argv != 1) then
     echo "usage: save_plate <plate>"
     echo "       where plate is <series><5 digit plate number>"
     exit
  endif
endif
echo "Creating /dasch/raid021/junk/${plate}"
mkdir /dasch/raid021/junk/${plate}
mkdir /dasch/raid021/junk/${plate}/bin9
mkdir /dasch/raid021/junk/${plate}/catalog9bin
mkdir /dasch/raid021/junk/${plate}/catalogall
mkdir /dasch/raid021/junk/${plate}/ingest
mkdir /dasch/raid021/junk/${plate}/match
mkdir /dasch/raid021/junk/${plate}/astrometry
mkdir /dasch/raid021/junk/${plate}/headers
echo "Copying files"
cp /dasch/Pipeline/bin9/${plate}* /dasch/raid021/junk/${plate}/bin9/
cp /dasch/Pipeline/catalog9bin/${plate}* /dasch/raid021/junk/${plate}/catalog9bin/
cp /dasch/Pipeline/catalogall/${plate}* /dasch/raid021/junk/${plate}/catalogall/
cp /dasch/Pipeline/ingest/${plate}* /dasch/raid021/junk/${plate}/ingest/
cp /dasch/Pipeline/match/${plate}* /dasch/raid021/junk/${plate}/match/
cp /dasch/Pipeline/match/match_${plate}* /dasch/raid021/junk/${plate}/match/
cp /dasch/Pipeline/astrometry/${plate}* /dasch/raid021/junk/${plate}/astrometry/
cp /dasch/Pipeline/astrometry/fail_${plate}* /dasch/raid021/junk/${plate}/astrometry/
cp /dasch/Pipeline/headers/${plate}* /dasch/raid021/junk/${plate}/headers/
echo "Making tarball"
cd /dasch/raid021/junk/${plate}/
tar -cvf /dasch/raid021/junk/${plate}.tar .
echo "Creating /dasch/raid021/junk/${plate}.tar.bz2"
ls -l /dasch/raid021/junk/${plate}.tar
bzip2 /dasch/raid021/junk/${plate}.tar
ls -l /dasch/raid021/junk/${plate}.tar.bz2
exit


