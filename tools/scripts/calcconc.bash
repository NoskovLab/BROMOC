#!/bin/bash
#paths
calcmu=calcmu.bash
mudatcalc=mudatcalc.bash
getres=getres.bash
calcf=calcf

#arguments
init=${1:?"Missing argument"}
end=${2:?"Missing argument"}
step=${3:?"Missing argument"}
from=${4:?"Missing argument"}
to=${5:?"Missing argument"}
stp=${6:?"Missing argument"}
boxx=${7:?"Missing argument"}
boxy=${8:?"Missing argument"}
boxz=${9:?"Missing argument"}
div1=${10:-100}
div2=${11:-100}

halfz=$($calcf $boxz % 2)
halfmz=$($calcf $boxz % -2)

# create template input file
cat << EOF > efpot.inp
TITLE name efpot
SYSTEM LX $boxx LY $boxy LZ $boxz temp 300.0 cdie 80.0
PARTICLES
POT  charge  1.0 diffusion 0.196
CLA  charge -1.0 diffusion 0.203
END
BUFFERS
POT  conc cncncn  mu mumumu  voltage 0.00  LZmin $halfmz  LZmax $halfz
CLA  conc cncncn  mu mumumu  voltage 0.00  LZmin $halfmz  LZmax $halfz
END
OPEN read unit 13 name    pot-pot.pot
OPEN read unit 14 name    cla-pot.pot
OPEN read unit 15 name    cla-cla.pot
EFPOT res 0.1
POT POT read unit 13
POT CLA read unit 14
CLA CLA read unit 15
END
CLOSE unit 13
CLOSE unit 14
CLOSE unit 15
COOR gener all
SIMULATION ncycle 750000 ngcmc 1 nbd 0 dt 0.02 nprint 1000
OPEN write unit 20 file name efpot.btr
SIMULATION ncycle 750000 ngcmc 1 nbd 0 dt 0.02 nsave 10 traject iuntrj 20 nprint 1000
CLOSE unit 20
EXIT
EOF

#main
for (( i=$init ; i<=$end ; i+=$step )); do
  conc=$($calcf $i % $div1)
  mkdir ./$i
  cd ./$i
  ln -s ../cla-cla.pot  
  ln -s ../cla-pot.pot  
  ln -s ../pot-pot.pot 
  sed 's/cncncn/'"$conc"'/g' ../efpot.inp > efpot.inp
  $calcmu $from $to $stp $boxx $boxy $boxz $div2
  $mudatcalc $div2
  $getres $div1
  cd ..
done
