#!/bin/bash
calc_mu=calc-mu
calcf=calcf
div1=${1:-100}
dirname=${PWD##*/}
conc=$($calcf $dirname % $div1)
$calc_mu mu-conc.dat 5 $conc > result
