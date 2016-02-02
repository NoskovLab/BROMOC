#!/bin/bash
fitit=fit-it
calcf=calcf
ini=${1:-1}
fin=${2:-3000}
step=${3:-1}
div1=${4:-1000}
rm -f table-conc-mu.dat
for (( i=$ini ; i<=$fin ; i=i+$step )); do
  echo -ne '\b\b\b\b\b\b'
  echo -n $i
  conc=$($calcf $i % $div1)
  $fitit conc-mu.dat 5 $conc > result
  cla=$(grep PA1 result | awk '{print $6}')
  pot=$(grep PA2 result | awk '{print $6}')
  echo $conc $cla $pot >> table-conc-mu.dat
done
echo
echo Done
rm -r result
