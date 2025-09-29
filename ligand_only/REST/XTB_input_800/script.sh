#!/bin/bash
export OMP_NUM_THREADS=4,1
export OMP_STACKSIZE=4G
#
mkdir ../XTB_output
#
for file in *.xyz; 
do
	base=$(basename "$file" .xyz)
	xtb "$file" --opt crude --alpb water --json
    mv xtbopt.xyz ../XTB_output/${base}_opt.xyz
	mv xtbout.json ../XTB_output/${base}_opt.json
done
