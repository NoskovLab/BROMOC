#!/bin/bash
#example
#getall.bash 100 1000 force
mudatcalc=mudatcalc.bash
getres=getres.bash
calcf=calcf
div1=${1:-100}
div2=${2:-100}
ini=${3:-1} 
fin=${4:-1000}
step=${5:-1}
echo $div1 $div2 $ini $fin $step
force="no"
if [ "$ini" == "force" ]; then force=yes ; ini=1 ; fi
rm -f conc-mu.dat
for (( i=$ini ; i<=$fin ; i+=$step )); do
if [ -d ./$i ]; then
  conc=$($calcf $i % $div1)
  echo ./$i
  cd ./$i
  if [ ! -e mu-conc.dat -o "$force" == "yes" ]; then $mudatcalc $div2 ;fi
  if [ ! -e result -o "$force" == "yes" ]; then $getres $div1; fi
  cla=$(grep PA1 result | awk '{print $6}')
  pot=$(grep PA2 result | awk '{print $6}')
  echo $conc $cla $pot >> ../conc-mu.dat
  cd ..
fi
done

