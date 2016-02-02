#!/bin/bash
stat=stat
calcf=calcf
div2=${1:-100}
ini=${2:-1}
fin=${3:-1000}
step=${4:-1}

rm -f mu-conc.dat
for (( i=$ini ; i<=$fin ; i=i+$step )); do
if [ -d ./$i ]; then 
  mu=$($calcf $i % -$div2)
  cd ./$i
  cla=$($stat conc-CLA-efpot.dat 2 | grep Media | awk '{print $2}')
  pot=$($stat conc-POT-efpot.dat 2 | grep Media | awk '{print $2}')
  echo $mu $cla $pot >> ../mu-conc.dat
  echo $i
  cd ..
fi
done
