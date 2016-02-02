pname=potmap-pbeq-read
gfortran -ffree-line-length-none -static -frecord-marker=8 -o $pname'-8bh' $pname'.f90'
gfortran -ffree-line-length-none -static -frecord-marker=4 -o $pname'-4bh' $pname'.f90'
mv  $pname'-8bh' ../bin
mv  $pname'-4bh' ../bin
rm -f ../../bin/$pname'-8bh'
rm -f ../../bin/$pname'-4bh'
ln -s ../tools/bin/$pname'-8bh' ../../bin/$pname'-8bh'
ln -s ../tools/bin/$pname'-4bh' ../../bin/$pname'-4bh'

