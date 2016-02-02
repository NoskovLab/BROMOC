#!/bin/bash
file=$1
num=$(grep ' - ' $file | wc -l )
lines=( $(grep -n ' - ' $file | sed 's/ //g' | sed 's/:/ /g' | awk '{print $1}') )
names=$(grep -n ' - ' $file | sed 's/ //g' | sed 's/:/ /g' | awk '{print $2}')
j=0
for i in $names; do
if [ $(($j+1)) -lt $num ]; then
  sed -n $((${lines[$j]}+1))','$((${lines[$(($j+1))]}-1))'p' $file > $i'.dat'
else
  sed -n $((${lines[$j]}+1))',$p' $file > $i'.dat'
fi
j=$(($j+1))
done

