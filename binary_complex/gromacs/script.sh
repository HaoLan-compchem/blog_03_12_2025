export OMP_NUM_THREADS=4
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#run gromacs with .top and .gro files converted from opemm 
#define index.ndx
echo '1 | 13' | gmx_mpi make_ndx -f gromacs.gro -o index.ndx
#write min.mdp configuration and run minimisation
gmx_mpi grompp -f min.mdp -c gromacs.gro -p gromacs.top -o min.tpr -r gromacs.gro -n index.ndx
gmx_mpi mdrun -v -deffnm min
#write gromacs.mdp configuration and run simulation (NPT) 
for i in $(seq 1 10); 
do
	mkdir replicate_$i
	cd replicate_$i
	cp ../gromacs.mdp ../min.gro ../gromacs.top ../index.ndx ./
	sed -i "/^gen_seed/c\gen_seed = $i" gromacs.mdp
	sed -i 's/\r//' gromacs.mdp #might need for WSL2
	gmx_mpi grompp -f gromacs.mdp -c min.gro -p gromacs.top -o md.tpr -r min.gro -n index.ndx -maxwarn 1
	gmx_mpi mdrun -v -deffnm md -nb gpu -bonded gpu -pme gpu -update gpu
	echo '11' | gmx_mpi energy -f md.edr -o potential.xvg
	cd ..
done
