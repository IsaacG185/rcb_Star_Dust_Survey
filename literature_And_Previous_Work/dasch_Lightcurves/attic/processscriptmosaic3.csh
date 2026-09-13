#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#############################
#  This script prepares a Linux.bat output for processing
#
#  Mar 11, 2011 Edward J. Los - Initial Version
#  Apr 16, 2011 Edward J. Los - Change /media to /dasch
#
###########################

if ($#argv != 2) then
  echo "usage: processscript.csh <yyyy_mm_dd> <disk number>"
  exit
endif

set uDate = $1
set raidArray = $2
set hDate = `echo $uDate | sed 's/_/-/g'`

echo "Date $uDate $hDate raidArray $raidArray"

echo "cd /dasch/mosaic3/ExposureData/Logs/$hDate/"
cd /dasch/mosaic3/ExposureData/Logs/$hDate/
echo "mkdir -p /dasch/TravelDrive/${uDate}/Logs"
mkdir -p /dasch/TravelDrive/${uDate}/Logs
cp * /dasch/TravelDrive/${uDate}/Logs
cd /dasch/TravelDrive/${uDate}/Logs
echo "cp Linux.bat ../process${uDate}.tmp2"
cp Linux.bat ../process${uDate}.tmp2
cd ..
dos2unix -f process${uDate}.tmp2
echo "exit\ndate\ncd /dasch/mosaic3/ExposureData\ncat /proc/diskstats\ncp -r * /dasch/raid0${raidArray}/ExposureData/\ncat /proc/diskstats\ndate\n" > process${uDate}.tmp1
echo "echo runpipeline -d ${hDate} -p 24 -a -w  -2 wd"
echo "cat /proc/diskstats\necho "\""cp -r * /dasch/raid017/ExposureData"\"" | at now\ndate\n" >  process${uDate}.tmp4

echo "cat /proc/diskstats\necho "\""runpipeline -d ${hDate} -p 24 -a -w  -2 wd >& /dasch/scanner/backup/${uDate}/run_xxx.log"\"" | at now\ndate\n" >  process${uDate}.tmp3
cat     process${uDate}.tmp1  process${uDate}.tmp2   process${uDate}.tmp3  >   process${uDate}.bat
rm process${uDate}.tmp1
rm process${uDate}.tmp2
rm process${uDate}.tmp3
rm process${uDate}.tmp4
echo "cp process${uDate}.bat ~"
cp process${uDate}.bat ~
cd ~
#emacs process${uDate}.bat
#at -f process${uDate}.bat now

