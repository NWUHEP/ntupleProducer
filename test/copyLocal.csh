#!/bin/csh

set username = andreypz

set dir = nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011B

set outDir = ~/eos/$dir
set inDir = /pnfs/cms/WAX/11/store/user/$username/$dir
echo $inDir
ls -d -1 $inDir/nuTuple_*.root  > DataFiles.txt

if (! -d $outDir) then
   mkdir $outDir
endif

foreach line (`cat DataFiles.txt`)
   echo "${line} to ${outDir}"
   dccp ${line} $outDir  
end

rm DataFiles.txt
