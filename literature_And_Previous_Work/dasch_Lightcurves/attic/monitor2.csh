#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# This utility monitors the state of runpipline
#
#
set startdate = `date`

echo "started on $startdate"

while (1) 
  set curdate = `date`

  set grepstr = `ps -ef | egrep "search|update|resort_magfiles" | grep "log" `
  if ($#grepstr < 6) then
    echo "Finished at $curdate"
    exit
  endif
  echo "$curdate $grepstr ---- $startdate"
  df -B 1000000 | grep "raid008"
  #free -m -t
  #vmstat
  sleep 60s
end
