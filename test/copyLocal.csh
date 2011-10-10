#!/bin/csh

set username = $1
set dir = $2

ls -d -1 /pnfs/cms/WAX/11/store/user/$username/$dir/nuTuple_*.root  > DataFiles.txt

if (! -d ~/nobackup/$dir) then
   mkdir ~/nobackup/$dir
endif

foreach line (`cat DataFiles.txt`)
   echo "${line} to ~/nobackup/$dir/"
   dccp ${line} ~/nobackup/$dir/  
end

rm DataFiles.txt
