#!/bin/bash
set -o nounset

PROGNAME=$(basename $0)

function usage
{
  cat <<EOF
Find a list of duplicate root files for a dataset at the SE that 
should be removed.

Usage: $PROGNAME -c <crab_dir> [--h | --help]
where options are:
  -c              Mandatory argument, crab project directory
  -v|--verbose    Turn on debug statements (D=false)
  -h|--help       This message

  example: $PROGNAME -c <crab_dir> -v

  This script creates two files in the present directory:

    allfiles.list - all the root files for the dataset present at the SE
    goodfiles.list - root files for successful jobs as found in the crab_fjr_n.xml files  

  and finds the duplicate files from the difference. Note, that at times jobs may finish
  and root files tranferred to the SE successfully, but crab may not immediately know about job
  completion. Those 'most recent' root files will be tagged as duplicate, but they
  are not.   
EOF

  exit 1
}

[ $# -gt 0 ] || usage

crab_dir=""
let "verbose = 0"
let "quiet = 0"
while [ $# -gt 0 ]; do
  case $1 in
    -c)                     shift
                            crab_dir=$1
                            ;;
    -v | --verbose )        let "verbose = 1"
                            ;;
    -q | --quiet )          let "quiet = 1"
                            ;;
    -h | --help )           usage
                            ;;
     * )                    usage
                            ;;
  esac
  shift
done

[ $crab_dir != "" ] && [ -e $crab_dir ] || usage

gflist=goodfiles.list
aflist=allfiles.list

# First of all get the list of goodfile by reading the fjr files
#export PERL5LIB=/afs/cern.ch/user/s/sarkar/public/perl/lib/perl5/site_perl/5.8.8:$PERL5LIB
#perl -w /afs/cern.ch/user/s/sarkar/public/ListGoodOutputFiles_new.pl $project/res > $gflist

[ $quiet -gt 0 ] || echo ">>> Find list of good files from fjr files..."
python find_goodfiles.py -c $crab_dir -q > $gflist 
# Now find the remote directory name
rdir=$(dirname $(head -1 $gflist))
srmp=$(echo $rdir | awk -F= '{print $1}')

# Get list of all files for the project
[ $quiet -gt 0 ] || echo ">>> Find list of all root files at $rdir ..."
srmls $rdir 2> /dev/null | grep '.root$' | awk '{if (NF==2) print $NF}' > $aflist

# Now compare
[ $quiet -gt 0 ] || echo ">>> Following is the list of duplicate files at $rdir ..."
for file in $(cat $aflist)
do
  grep $file $gflist > /dev/null
  [ $? -eq 0 ] && continue

  bname=$(basename $file)
  grep $bname $gflist > /dev/null
  [ $? -eq 0 ] && continue

  echo "$srmp""=""$file"
done

exit 0
