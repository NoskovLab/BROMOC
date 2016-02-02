#!/bin/bash
#paths
bromoc=bromoc
concbromoc=conc-bromoc
calcf=calcf
from=${1:?"Argument missing"}
to=${2:?"Argument missing"}
step=${3:?"Argument missing"}
boxx=${4:?"Argument missing"}
boxy=${5:?"Argument missing"}
boxz=${6:?"Argument missing"}
div2=${7:-100}
dirname=${PWD##*/}
for (( i=$from ; i<=$to ; i+=$step )); do
  mu=$($calcf $i % -$div2)
  echo $i $mu
  mkdir ./$i
  cd ./$i
  ln -s ../cla-cla.pot  
  ln -s ../cla-pot.pot  
  ln -s ../pot-pot.pot 
  ln -s ../static.pbeq
  ln -s ../rfpar.pbeq
  sed 's/mumumu/'"$mu"'/g' ../efpot.inp > efpot.inp
  $bromoc < efpot.inp > efpot.out
  $concbromoc efpot.btr box $boxx $boxy $boxz > /dev/null 
  rm efpot.btr
  cd ..
done

