#!/bin/bash

for file in $(ls *.hdf)
do

  hdf32to64 $file $file

done
