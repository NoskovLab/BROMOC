#!/bin/bash
nat1=$(head -1 $1)
nat2=$(sed -n $(($nat1+2))'p' $1)
nat=$(($nat1+$nat2))
echo $nat > $2
echo >> $2
#head -$nat1 $1 | tail -$(($nat)) | awk '{print $3"   "$4"   "$5"   "$6}' >> $2
sed -n 2','$(($nat1+1))'p' $1 | awk '{print substr($0, 5, 80)}' >> $2
sed -n $(($nat1+3))',$p' $1 | awk '{print substr($0, 10, 50)}' | sed "h;s/38/K /g;s/39/Cl/g;G;s/\(..\).*\n../\1/" >> $2
