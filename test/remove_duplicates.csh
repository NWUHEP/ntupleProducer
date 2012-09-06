#!/bin/csh

#set mypath = "nuTuples_v1_7TeV/VBFHZZ250"
#set mypath = "/eos/uscms/store/user/bpollack/May15/MC/WWJets"
set mypath = "/pnfs/cms/WAX/11/store/user/naodell/March18/MC/VBFZ"

set nfiles = 101
set count = 1
while ($count <= $nfiles)
   set nlines = `ll -h ${mypath}/nuTuple_${count}_*.root | wc -l`
   #echo $nlines
   if ( $nfiles == 0 ) then
    echo "File not found for # " $count
   endif
   if ( $nlines > 1) then
       #echo $count $nlines
       ll -ht ${mypath}/nuTuple_${count}_*.root > listoffiles.txt
       set count2 = 1 
       while ($count2 <= $nlines)
          set fname = `cat listoffiles.txt | head -n $count2 | tail -1 | cut -f9 -d' '`
          echo $fname 
               if ( $count2>1) then
                   echo "To be removed: " $fname 
                   rm $fname 
               endif
        @ count2 ++
        end
   endif
 @ count ++
end
