#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# setenvironment.csh 
#
#  Used in conjunction with runpipeline to set up multiple hosts on Odyssey
#   
#  mpirun --hostfile /dasch/Pipeline/mpihostalldasch -np 8 runpipeline -g 'source setenvironment.csh' 
#
#   mpirun --hostfile /n/dasch8/Pipeline/Pipeline/mpihostalldasch -np 8 /n/dasch8/Pipeline/Pipeline/runpipeline -g '/n/dasch8/Pipeline/Pipeline/setenvironment.csh' > /n/dasch8/Pipeline/Pipeline/losafter2.log
#
#  Nov 11, 2008 Edward J. Los - Initial version
#  Mar  7, 2016 Edward J. Los - Add dasch12, remove dasch3,dasch4, and dasch5
#
set dasch6base = "/n/dasch6"
set dasch7base = "/n/dasch7"
set dasch8base = "/n/dasch8"
set dasch9base = "/n/dasch9"
set dasch10base = "/n/dasch10"
set dasch11base = "/n/dasch11"
set dasch12base = "/n/dasch12"
set dasch13base = "/n/dasch13"

if ($HOST == "dasch6.rc.fas.harvard.edu") then
   set localmachine = "dasch6"
   set dasch6base = "/dasch6"
else if ($HOST == "dasch7.rc.fas.harvard.edu") then
   set localmachine = "dasch7"
   set dasch7base = "/dasch7"
else if ($HOST == "dasch8.rc.fas.harvard.edu") then
   set localmachine = "dasch8"
   set dasch8base = "/dasch8"
   set dasch11base = "/dasch11"
else if ($HOST == "dasch9.rc.fas.harvard.edu") then
   set localmachine = "dasch9"
   set dasch9base = "/dasch9"
else if ($HOST == "dasch10.rc.fas.harvard.edu") then
   set localmachine = "dasch10"
   set dasch10base = "/dasch10"
else if ($HOST == "dasch12.rc.fas.harvard.edu") then
   set localmachine = "dasch12"
   set dasch12base = "/dasch12"
else if ($HOST == "dasch13.rc.fas.harvard.edu") then
   set localmachine = "dasch13"
   set dasch13base = "/dasch13"
endif 
echo "localmachine $localmachine Host: $HOST"
echo "dasch6base $dasch6base dasch7base $dasch7base dasch8base $dasch8base dasch9base $dasch9base dasch10base $dasch10base dasch11base $dasch11base dasch12base $dasch12base dasch13base $dasch13base "
if (0) then
   cd /dasch
    cp /dasch/raid021/junk/cshrc_dasch /dasch/cshrc_dasch
    diff -b -w -s /dasch/raid021/junk/cshrc_dasch ./cshrc_dasch 
endif
echo "1"
if (0) then
   cd /dasch
   rm  raid016 
   rm  scannerkits1 
   rm  scannerkits2 
   rm  scannerkits3 
   rm  raid017 
   rm  raid018 
   rm  backup 
   ln -s  ${dasch8base}/scanner/backup backup
   rm  common 
   ln -s  ${dasch8base}/scanner/common common
   rm  logindex 
   ln -s  ${dasch8base}/scanner/auxdocs/logindex logindex
   rm  mysql 
   ln -s  ${dasch8base}/scanner/mysql mysql
   rm  perl 
   ln -s  ${dasch8base}/scanner/linux/perl perl
   rm  web 
   ln -s  ${dasch8base}/scanner/web web
   cd /dasch
   rm  catalogs 
   ln -s  /${localmachine}/catalogs catalogs
   rm  ucac3 
   ln -s  /${localmachine}/ucac3 ucac3
   rm  ucac4 
   ln -s  /${localmachine}/ucac4 ucac4
   rm  ucac5 
   ln -s  /${localmachine}/ucac5 ucac5
   rm  scratch 
   ln -s  /${localmachine}/scratch scratch
   rm  install 
   ln -s  ${dasch8base}/install install
   rm  junk 
   ln -s  ${dasch8base}/junk junk
   rm  build 
   ln -s  ${dasch8base}/Pipeline/build build
   rm  apps 
   ln -s  ${dasch8base}/Pipeline/build/apps apps
   rm  cfitsio 
   ln -s  ${dasch8base}/Pipeline/build/cfitsio cfitsio
   rm  wcstools 
   ln -s  ${dasch8base}/Pipeline/build/wcstools wcstools
   rm  Pipeline 
   ln -s  ${dasch8base}/Pipeline/Pipeline Pipeline
   rm  raid019 
   ln -s  ${dasch6base} raid019
   rm  raid020 
   ln -s  ${dasch7base} raid020
   rm  raid021 
   ln -s  ${dasch8base} raid021
   rm  raid022 
   ln -s  ${dasch9base} raid022
   rm  raid023
   ln -s  ${dasch10base} raid023
   rm  raid024 
   ln -s  ${dasch11base} raid024
   rm  raid025
   ln -s  ${dasch12base} raid025
   rm  raid026
   ln -s  ${dasch13base} raid026
   rm  raid027
   ln -s /n/boslfs/LABS/dasch_project/dasch14 raid027
   rm  raid028
   ln -s /n/boslfs/LABS/dasch_project/dasch15 raid028
   rm  raid029
   ln -s /n/boslfs/LABS/dasch_project/dasch16 raid029
   rm  gsc23 
   ln -s  ${dasch9base} gsc23
   rm  raid001 
   ln -s  ${dasch9base}/raid001 raid001
   rm  raid002 
   ln -s  ${dasch9base}/raid002 raid002
   rm  raid003 
   ln -s  ${dasch9base}/raid003 raid003
   rm  mosaic 
   ln -s  ${dasch9base} mosaic
   rm  magnitudes 
   ln -s  ${dasch11base}/mysql/photfiles magnitudes
   rm  magnitudes1 
   ln -s  ${dasch11base}/mysql/photfiles1 magnitudes1
   rm  magnitudes2 
   ln -s  ${dasch11base}/mysql/photfiles2 magnitudes2
   rm  magnitudes3 
   ln -s /n/dasch14/mysql/photfiles3 magnitudes3
endif
echo "2"
if (0) then
   cd /dasch
   if ($localmachine != "dasch8") then
     rm  localPipeline 
     ln -s  /${localmachine}/localPipeline localPipeline
   else 
     echo "No localPipeline on $localmachine"
   endif
   rm  scanner 
   ln -s  ${dasch8base}/scanner scanner
   rm  scannerkits 
   ln -s  ${dasch8base}/scannerkits scannerkits
   rm  wcsfit 
   ln -s  ${dasch8base}/wcsfit wcsfit

endif

echo "3"

if (0) then
   cd /dasch
   ls -l | grep ">"
endif
echo "4"

if (0) then
    cd /dasch/catalogs
    pwd
    #mv apass.dat apass_dr7.dat
    #mv apass.idx apass_dr7.idx
    if ($HOST != "dasch8.rc.fas.harvard.edu") then
        echo "cp /dasch/raid021/Pipeline/catalogs/apass.dat /dasch/catalogs/"
        cp /dasch/raid021/Pipeline/catalogs/apass.dat /dasch/catalogs/
        cp /dasch/raid021/Pipeline/catalogs/apass.idx /dasch/catalogs/
    endif 
    ls -l apass*.dat
    ls -l apass*.idx
endif
echo "5"

if (1) then
    cd /${localmachine}/
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_${localmachine}.txt
    if ($HOST == "dasch8.rc.fas.harvard.edu") then
    cd /dasch11/
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_dasch11.txt 
    df -B 1000000
    date
    cd /n/dasch16/
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_dasch16.txt 
    date
    endif
    if ($HOST == "dasch12.rc.fas.harvard.edu") then
    cd /n/dasch13/
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_dasch13.txt 
    df -B 1000000
    endif
    if ($HOST == "dasch10.rc.fas.harvard.edu") then
    cd /n/dasch14/
    date
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_dasch14.txt 
    endif
    date
    if ($HOST == "dasch9.rc.fas.harvard.edu") then
    date
    cd /n/dasch15/
    pwd
    du -B 1000000 --max-depth=2 > /dasch/raid021/junk/du_dasch15.txt 
    date
    endif


    
endif
echo "6"

if (0) then
   cd /dasch
   rm  magnitudes
   ln -s  ${dasch11base}/mysql/photfiles magnitudes
   rm  magnitudes1 
   ln -s  ${dasch11base}/mysql/photfiles1 magnitudes1
   rm  magnitudes2 
   ln -s  ${dasch11base}/mysql/photfiles2 magnitudes2
   rm  magnitudes3 
    ln -s /n/dasch14/mysql/photfiles3 magnitudes3
    ls -l /dasch/mag* 
endif
echo "7"

if (0) then
    df -B 1000000 /n/dasch12
endif
echo "8"
if (0) then
echo "dasch6"
ls  /dasch/raid019/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid020/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid021/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid022/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid023/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid024/ExposureData/Mosaics/
echo "dasch6"
ls  /dasch/raid025/ExposureData/Mosaics/
echo "raid026"
ls  /dasch/raid026/ExposureData/Mosaics/
echo "dasch13"
ls /n/dasch13
echo "end"
endif
echo "9"
if (0) then
   if ($localmachine == "dasch8") then
    echo  "skipping dasch8"
  else
    echo "$localmachine"
    echo "du output:"
    cd /${localmachine}/localPipeline
    du -B 1000000 
    echo "find output:"
    find . -type f | wc
  endif
endif
if (0) then
    dmesg
endif
echo "10"
if (0) then
    cd $DASCH_SCRATCH
    du -B 1000000 --summarize
    find . -name "*.fit"
endif
if (0) then
    cd /${localmachine}/
    pwd
    #find ./ExposureData/Mosaics -name "*.fit" -exec ls -l {} \;
    find ./ExposureData/Mosaics -type d -exec ls -l {} \;
    if ($HOST == "dasch8.rc.fas.harvard.edu") then
    cd /dasch11/
    pwd
    #find ./ExposureData/Mosaics -name "*.fit" -exec ls -l {} \;
    find ./ExposureData/Mosaics -type d -exec ls -l {} \;
    endif
    if ($HOST == "dasch12.rc.fas.harvard.edu") then
    cd /n/dasch13/
    pwd
    #find ./ExposureData/Mosaics -name "*.fit" -exec ls -l {} \;
    find ./ExposureData/Mosaics -type d -exec ls -l {} \;
    endif
endif
echo "1"
exit
