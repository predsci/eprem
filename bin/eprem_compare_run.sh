#!/bin/bash

if [ ! $# -eq 2 ]
then
  echo "USAGE:  eprem_compare_run.sh RUNDIR_REF RUNDIR_NEW"
  exit 1
fi

rundir1=$1
rundir2=$2

if [ ! -d $rundir1 ]
then
  echo "ERROR $rundir1 does not exist!"
  exit 1
fi

if [ ! -d $rundir2 ]
then
  echo "ERROR $rundir2 does not exist!"
  exit 1
fi

trap ctrl_c INT

function ctrl_c() {
  echo "Caught CTRL-C, shutting down!"
  exit 1
}



seperator_str="--------------------------------------------------------------------------------------------------------------------"

echo $seperator_str

echo "Comparing EPREM runs:"
echo "a: ${rundir1}"
echo "b: ${rundir2}"

echo $seperator_str

failed=0

i=0
for fullfn in $(ls ${rundir1}/obs*.nc)
do
  if [ $i -eq 0 ]
  then
    i=1
    fn=${fullfn#${rundir1}/}
    if [ ! -e ${rundir2}/$fn ]
    then
      echo "ERROR! $fn does not exist in ${rundir2}!"
      exit 1
    fi
    epremdiff.py ${rundir1}/$fn ${rundir2}/$fn > temp.txt
    header_str=$(tail -n 4 temp.txt | head -n 1)
    echo "Filename     $header_str"
    echo $seperator_str
  fi
  fn=${fullfn#${rundir1}/}
  if [ ! -e ${rundir2}/$fn ]
  then
    echo "ERROR! $fn does not exist in ${rundir2}!"
    exit 1
  fi
  epremdiff.py ${rundir1}/$fn ${rundir2}/$fn > temp.txt
  result_str=$(tail -n 2 temp.txt | head -n 1)
  passfail_str=$(tail -n 1 temp.txt | head -n 1)
  echo "${fn}    $result_str    $passfail_str"
  if [ "$passfail_str" = "FAIL" ]
  then
    failed=1
  fi
done
rm temp.txt
echo $seperator_str
if [ $failed -eq 1 ]
then
  echo "FAIL"
else
  echo $passfail_str
fi
echo $seperator_str

exit 0
