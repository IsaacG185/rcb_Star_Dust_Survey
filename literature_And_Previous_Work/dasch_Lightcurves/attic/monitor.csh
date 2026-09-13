#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# This utility monitors the state of runpipline
#
#
if ($#argv != 1) then
  echo "usage: monitor.csh <file list>"
  exit
endif
set filelist = $1
set startdate = `date`

set listsize = `wc -l $1`

echo "$listsize started on $startdate"

while (1) 
  set curdate = `date`
  set filestr = `ls -rt *.xxx | gawk '{print $1}' | head -n 1 | sed 's/.xxx//g'`
  if ($#filestr != 1) then
    echo "Finished at $curdate"
    exit
  endif

  set grepstr = `grep -n $filestr $filelist`
  if ($#grepstr != 1) then
    echo "Finished at $curdate searching for $filestr"
    exit
  endif
  echo "$curdate $grepstr ---- $startdate $listsize"
  #free -m -t
  #vmstat
  sleep 60s
end
