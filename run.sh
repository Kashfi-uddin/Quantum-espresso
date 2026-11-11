#!/bin/bash

for i in {001..005};
do
cat header.in supercell-$i.in > pwscf-$i.in;
export NUM_OMP_THREADS=2
mpirun -np 18 pw.x -in pwscf-$i.in | tee pwscf-$i.out;
done
