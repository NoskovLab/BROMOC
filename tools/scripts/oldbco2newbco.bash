#!/bin/bash
file=$1
ofile=$2
numl=$(head -1 $file)
head -1 $file > $ofile
endl=$(($numl+1))
i=$endl
while [ $i -ge 6 ] ; do
  ln2=$(($i-1))
  ln3=$(($i-2))
  sed -n $ln2','$i'p' $file >> $ofile
  sed -n $ln3'p' $file >> $ofile
  i=$(($i-3))
done
sed -n '2,3p' $file >> $ofile
tail -1 $file >> $ofile

