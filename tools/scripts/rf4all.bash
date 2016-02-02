#!/bin/bash
calcf=calcf
calc=calc
box=${1:-60x60x60}
div1=${2:-100}
conc=${3:-10,20,40,60,80,100,125,150}

boxvec=( $(echo $box | sed 's/x/ /g') )
dir=$(echo $conc | sed 's/,/ /g')
ncl=( $(for (( i=0; i<=2; i++ )); do echo -n $($calc $($calc ${boxvec[$i]} + 2) x 2)' ' ; done) )
bmm=( $(for (( i=0; i<=2; i++ )); do echo -n $($calcf ${boxvec[$i]} % 2 | sed 's/\.0000000000000000/\.0/g')' ' ; done) )

#create dir
mkdir $box
cd $box

#create static input
cat << EOF > static.inp
IONTYPE
POT  charge  1.0  
CLA  charge -1.0  
END
PHYSICAL_PARAMETERS epsp 2.0  epsw 80.0  temp 300.0  conc cccc  watr 0.0  ionr 0.0
OBJECTS rsphe 0.0 xsphe 0.0 ysphe 0.0 zsphe 0.0 epss 0.0 -
        rcyln 0.0 hcyln 0.0 xcyln 0.0 ycyln 0.0 zcyln 0.0 epsc 0.0 -
        bxmax ${bmm[0]} bymax ${bmm[1]} bzmax ${bmm[2]} bxmin -${bmm[0]} bymin -${bmm[1]} bzmin -${bmm[2]} epsb 0.0
GRID nclx ${ncl[0]}  ncly ${ncl[1]} nclz ${ncl[2]}  dcel 1.5  xbcen 0.0  ybcen 0.0  zbcen 0.0
PBEQ
GRID nclx ${ncl[0]}  ncly ${ncl[1]} nclz ${ncl[2]}  dcel 0.5  xbcen 0.0  ybcen 0.0  zbcen 0.0
PBEQ  phifocus -
      pnclx ${ncl[0]} pncly ${ncl[1]} pnclz ${ncl[2]} pdcel 1.5 pxbcen 0.0 pybcen 0.0 pzbcen 0.0
open write unit 14 file name static.pbeq
write phi  unit 14
STOP
EOF

#create rfpar input
cat << EOF > rfpar.inp
IONTYPE
POT  charge  1.0 radius 0.5
END
PHYSICAL_PARAMETERS epsp 2.0  epsw 80.0  temp 300.0  conc cccc  watr 0.0  ionr 0.0
OBJECTS rsphe 0.0  xsphe 0.0  ysphe 0.0  zsphe 0.0  epss 0.0 -
        rcyln 0.0  hcyln 0.0  xcyln 0.0  ycyln 0.0  zcyln 0.0 epsc 0.0 -
        bxmax ${bmm[0]} bymax ${bmm[1]} bzmax ${bmm[2]} bxmin -${bmm[0]} bymin -${bmm[1]} bzmin -${bmm[2]} epsb 0.0
GRID nclx ${ncl[0]}  ncly ${ncl[1]}  nclz ${ncl[2]}  dcel 0.5  xbcen 0.0  ybcen 0.0  zbcen 0.0
RFPAR box
OPEN write unit 15 file name rfpar.pbeq
WRITE rfpar POT  unit 15
STOP
EOF

# main
for i in $dir; do
  mkdir $i
  cd $i
  conc=$($calcf $i % 100)
  echo $i static
  sed 's/cccc/'"$conc"'/g' ../static.inp > static.inp
  pb-pnp < static.inp > static.out &
  echo $i rfpar
  sed 's/cccc/'"$conc"'/g' ../rfpar.inp > rfpar.inp
  pb-pnp < rfpar.inp > rfpar.out & 
  cd ..
done 


