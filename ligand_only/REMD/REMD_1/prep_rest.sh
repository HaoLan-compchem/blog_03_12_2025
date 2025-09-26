export OMP_NUM_THREADS=8
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#echo 'q' | gmx_mpi make_ndx -f gromacs.gro -o index.ndx
#gmx_mpi grompp -f min.mdp -c gromacs_ligand.gro -p gromacs_ligand.top -o min.tpr -r gromacs_ligand.gro -n index.ndx
#gmx_mpi mdrun -v -deffnm min

T_min=300
alpha=1.04729
n_replicas=16

for ((i=6; i<n_replicas; i++)); do
    mkdir replicate_$i
    cd replicate_$i
    cp ../gromacs_ligand.mdp ../nvt.mdp ../min.gro ../gromacs_ligand.top ../index.ndx ./
    sed -i "s/^ref_t.*/ref_t = $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l)) $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l))/" gromacs_ligand.mdp
    sed -i "s/^ref_t.*/ref_t = $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l)) $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l))/" nvt.mdp
    gmx_mpi grompp -f nvt.mdp -c min.gro -p gromacs_ligand.top -o nvt.tpr -r min.gro -n index.ndx -maxwarn 1
    gmx_mpi mdrun -v -deffnm nvt
    gmx_mpi grompp -f gromacs_ligand.mdp -c nvt.gro -p gromacs_ligand.top -o md.tpr -r nvt.gro -n index.ndx -maxwarn 1
    cd ..
done

export OMP_NUM_THREADS=1
mpirun -np 16 gmx_mpi mdrun -multidir replicate_0 replicate_1 replicate_2 replicate_3 replicate_4 replicate_5 replicate_6 replicate_7 replicate_8 replicate_9 replicate_10 replicate_11 replicate_12 replicate_13 replicate_14 replicate_15 -s md.tpr -deffnm md -replex 500
