#!/bin/csh

set nfiles = 402
set count = 1
set name = "/pnfs/cms/WAX/11/store/user/andreypz/Data3/DoubleMu_PromptV4/nuTuple"
while ($count <= $nfiles)
   set nlines = `ll -h ${name}_${count}_*.root | wc -l`
   if ( $nlines > 1) then
       echo $count $nlines
       ll -ht ${name}_${count}_*.root > listoffiles.txt
       set count2 = 1 
       while ($count2 <= $nlines)
          set fname = `cat listoffiles.txt | head -n $count2 | tail -1 | cut -f9 -d' '`
          echo $fname 
               if ( $count2>1) then
                   echo "To be removed: " $fname 
#                  rm $fname 
               endif
        @ count2 ++
        end
   endif
 @ count ++
end

