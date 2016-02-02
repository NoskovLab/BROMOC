#!/bin/bash
min () {
if [ $1 -lt $2 ]; then
  echo $1
else
  echo $2
fi
}

max () {
if [ $1 -gt $2 ]; then
  echo $1
else
  echo $2
fi
}

index () { 
local i=$(min $1 $2)
local j=$(max $1 $2)
local n=$3
echo $(($n*($i-1)-($i-1)*($i-2)/2+$j-$i+1)) 
}

filename=dna-kcl-in.pot
output=dna-ions-in.pot
#actualorder='PHO SUGA ADEN CYTO GUAN THYM CLA POT'
#neworder='ADEN CYTO GUAN THYM SUGA PHO KJ CLJ'
#iordernum='6 5 1 2 3 4 8 7'
ordernum='3 4 5 6 2 1 8 7'
iov=( $ordernum )
n=$(head -1 $filename | awk '{print $1}')
lines=$(head -1 $filename | awk '{print $2}')
#nn=$(($n*($n+1)/2))
head -1 $filename > $output
for (( i = 1 ; i <= $n ; i++ )); do
  for (( j = $i ; j <= $n ; j++ )); do
    ii=${iov[$(($i-1))]}
    jj=${iov[$(($j-1))]}
    ind=$(index $ii $jj $n)
    ind2=$(index $i $j $n)
    echo $ii $jj $i $j $n $ind $ind2 
    sed -n $((($ind-1)*$lines+2))','$(($ind*$lines+1))'p' $filename | awk '{print substr($0,1,41)}' > tmptmp1
    sed -n $((($ind2-1)*$lines+2))','$(($ind2*$lines+1))'p' $filename | awk '{print substr($0,42)}' > tmptmp2
    paste tmptmp1 tmptmp2 > tmptmp3
    cat tmptmp2 >> $output
  done
done
rm -f tmptmp1 tmptmp2 tmptmp3


