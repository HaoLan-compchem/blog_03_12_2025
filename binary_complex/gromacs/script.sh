export OMP_NUM_THREADS=4
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
for i in $(seq 1 10); 
do
	cd replicate_$i
	sed -i "/^gen_seed/c\gen_seed = $i" gromacs.mdp
	sed -i 's/\r//' gromacs.mdp
	gmx_mpi grompp -f gromacs.mdp -c min.gro -p gromacs.top -o md.tpr -r min.gro -n index.ndx -maxwarn 1
	gmx_mpi mdrun -v -deffnm md -nb gpu -bonded gpu -pme gpu -update gpu
	echo '11' | gmx_mpi energy -f md.edr -o potential.xvg
	cd ..
done
