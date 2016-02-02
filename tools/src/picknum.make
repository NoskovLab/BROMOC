#!/bin/bash
if [ "$1" == "intel" ]; then
  pname=$(basename $(echo $0 | sed 's/\(.*\)\..*/\1/' ))
  mv $pname.f90 $pname.f90.tmp
  sed 's/!use ifport/use ifport/g' $pname.f90.tmp > $pname.f90
fi
./make.generic $0 "$1"
if [ "$1" == "intel" ]; then
  mv $pname.f90.tmp $pname.f90
fi


