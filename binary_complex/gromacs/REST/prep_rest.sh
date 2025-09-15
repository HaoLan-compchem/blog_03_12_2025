export OMP_NUM_THREADS=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#echo 'q' | gmx_mpi make_ndx -f gromacs.gro -o index.ndx
#gmx_mpi grompp -f min.mdp -c gromacs_ligand.gro -p gromacs_ligand.top -o min.tpr -r gromacs_ligand.gro -n index.ndx
#gmx_mpi mdrun -v -deffnm min

n_replicas=12
BC_SCALE=6
MIN_LAMBDA=0.6

for ((i=0; i<n_replicas; i++)); do
    mkdir replicate_$i
    cd replicate_$i
    cp ../gromacs_ligand.mdp ../min.gro ../gromacs_ligand.top ../index.ndx ./
    if [[ $i -eq $((n_replicas - 1)) ]]; then
      LAMBDA_VAL=1.0
    else
      LAMBDA_VAL=$(echo "scale=$BC_SCALE; $MIN_LAMBDA^($i/($n_replicas-1))" | bc -l)
    fi

    printf "\ncouple-lambda = %s\n" "$LAMBDA_VAL" >> gromacs_ligand.mdp
    gmx_mpi grompp -f gromacs_ligand.mdp -c min.gro -p gromacs_ligand.top -o md.tpr -r min.gro -n index.ndx -maxwarn 1
    cd ..
done
