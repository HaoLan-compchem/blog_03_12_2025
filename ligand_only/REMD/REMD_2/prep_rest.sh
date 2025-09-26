#export OMP_NUM_THREADS=8
#export GMX_ENABLE_DIRECT_GPU_COMM=1
#export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#echo 'q' | gmx_mpi make_ndx -f gromacs_ligand.gro -o index.ndx
#gmx_mpi grompp -f min.mdp -c gromacs_ligand.gro -p gromacs_ligand.top -o min.tpr -r gromacs_ligand.gro -n index.ndx
#gmx_mpi mdrun -v -deffnm min

#T_min=300
#N=16
#alpha=1.03

#for ((i=0; i<N; i++)); do
    #mkdir replicate_$i
    #cd replicate_$i
    #cp ../gromacs_ligand.mdp ../nvt.mdp ../min.gro ../gromacs_ligand.top ../index.ndx ./
    #sed -i "s/^ref_t.*/ref_t = $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l)) $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l))/" gromacs_ligand.mdp
    #sed -i "s/^ref_t.*/ref_t = $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l)) $(printf "%.0f\n" $(echo "$T_min * $alpha ^ $i" | bc -l))/" nvt.mdp
    #gmx_mpi grompp -f nvt.mdp -c min.gro -p gromacs_ligand.top -o nvt.tpr -r min.gro -n index.ndx -maxwarn 1
    #gmx_mpi mdrun -v -deffnm nvt
    #gmx_mpi grompp -f gromacs_ligand.mdp -c nvt.gro -p gromacs_ligand.top -o md.tpr -r nvt.gro -n index.ndx -maxwarn 1
    #cd ..
#done

#export OMP_NUM_THREADS=1
#mpirun -np 16 gmx_mpi mdrun -multidir replicate_{0..15} -s md.tpr -deffnm md -replex 500

#post-process trajactory at T_min (swapped)
echo '1 0' | gmx_mpi trjconv -s replicate_0/md.tpr -f replicate_0/md.xtc -o replicate_0/md_pbc.xtc -pbc mol -center
echo '1 0' | gmx_mpi trjconv -s replicate_0/md.tpr -f replicate_0/md_pbc.xtc -o replicate_0/md_pbc_fit.xtc -fit rot+trans

#post-process trajactory from replicate index 0 (unswapped)
cat replicate_*/md.log > full_remd.log
demux.pl full_remd.log
gmx_mpi trjcat -f replicate_{0..15}/md.xtc -demux replica_index.xvg -o sorted_traj.xtc
echo '1 0' | gmx_mpi trjconv -s replicate_0/md.tpr -f 0_sorted_traj.xtc -o 0_sorted_md_pbc.xtc -pbc mol -center
echo '1 0' | gmx_mpi trjconv -s replicate_0/md.tpr -f 0_sorted_md_pbc.xtc -o 0_sorted_md_pbc_fit.xtc -fit rot+trans

