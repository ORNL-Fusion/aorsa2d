#!/bin/csh
#
F90 -c fixcon.f90 >& compiler.txt
if ( $status != 0 ) then
  echo "Errors while compiling fixcon.f90"
  exit
endif
rm compiler.txt
#
F90 fixcon.o
if ( $status != 0 ) then
  echo "Errors while loading fixcon.o"
  exit
endif
rm fixcon.o
#
mv a.out ~/bin/$ARCH/fixcon
#
echo "Executable installed as ~/bin/$ARCH/fixcon."
