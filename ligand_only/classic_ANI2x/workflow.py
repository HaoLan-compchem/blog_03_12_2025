import sys
import torch
import numpy as np
#from torchani.models import ANI2x
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm.app import PDBFile
#from openmmtorch import TorchForce
from openmmml import MLPotential # this module include functions of torch, torchani and openmmtorch
#
cuda_available = torch.cuda.is_available()
print(f"CUDA Available: {cuda_available}")

if cuda_available:
    num_gpus = torch.cuda.device_count()
    print(f"Number of GPUs: {num_gpus}")
    current_gpu_name = torch.cuda.get_device_name(0)
    print(f"Current GPU Name: {current_gpu_name}")
    print("\nGPU Memory Usage:")
    print(f"  Allocated: {round(torch.cuda.memory_allocated(0) / 1024**3, 1)} GB")
    print(f"  Cached: {round(torch.cuda.memory_reserved(0) / 1024**3, 1)} GB")
else:
    print("PyTorch is running on CPU.")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"\nUsing device: {device}")
#
with open('gaff_system.xml', 'r') as f:
    system_1 = mm.XmlSerializer.deserialize(f.read())
pdb = PDBFile('gaff_ligand_in_solvent.pdb')
#system_1.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin))
#
chains = list(pdb.topology.chains())
ml_atoms = [atom.index for atom in chains[0].atoms()]
#
potential = MLPotential('ani2x', implementation='nnpops')
mixed_system = potential.createMixedSystem(pdb.topology, system_1, ml_atoms)
#
# Create an integrator with a time step of 0.5 fs (ensure NN potential gradient force precision and system stability)
temperature = 300 * unit.kelvin
frictionCoeff = 1 / unit.picosecond
timeStep = 0.5 * unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, frictionCoeff, timeStep)

# Create a simulation and set the initial positions and velocities
simulation = app.Simulation(pdb.topology, mixed_system, integrator)
print("Simulation is running on:", simulation.context.getPlatform().getName())
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)

# Minimize the system
simulation.minimizeEnergy()
with open(f'min_complex.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
