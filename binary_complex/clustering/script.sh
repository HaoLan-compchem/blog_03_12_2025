#!/bin/bash
gmx_mpi trjcat -f ../gromacs/all_processed.xtc ../openmm/all_processed.xtc -o merged.xtc -cat
echo '3 0' | gmx_mpi trjconv -f merged.xtc -s ../openmm/complex_in_solvent.pdb -n index.ndx -o merged_fit.xtc -fit rot+trans -skip 10
echo '13 0' | gmx_mpi cluster -f merged_fit.xtc -s ../openmm/complex_in_solvent.pdb -n index.ndx -cutoff 0.2 -method gromos -wcl 5 -nst 50 -nofit
echo '13 0' | gmx_mpi cluster -f merged_fit.xtc -s ../openmm/complex_in_solvent.pdb -n index.ndx -cutoff 0.2 -method gromos -wcl 5 -nst 50 -nofit -cl clusters.gro