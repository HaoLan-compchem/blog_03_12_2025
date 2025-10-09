import sys
import torch
import numpy as np
#
import openmm as mm
import openmm.unit as unit
import openmm.app as app
from openmm.app import PDBFile
#
from openmmtorch import TorchForce
from torchani.models import ANI2x

#-----------------------------------------------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------------------------------------------#

with open('gaff_system.xml', 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
pdb = PDBFile('gaff_ligand_in_solvent.pdb')
system.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin))
#
chains = list(pdb.topology.chains())
print(chains)
ml_atoms = [atom.index for atom in chains[0].atoms()]
print(ml_atoms)
atomic_numbers = [atom.element.atomic_number for atom in chains[0].atoms()]
print(atomic_numbers)

#-----------------------------------------------------------------------------------------------------------------#

def removeBonds(system, atoms, removeInSet=True, removeConstraints=True):

    atomSet = set(atoms)
    import xml.etree.ElementTree as ET
    xml = mm.XmlSerializer.serialize(system)
    root = ET.fromstring(xml)

    def shouldRemove(termAtoms):
        return all(a in atomSet for a in termAtoms) == removeInSet

    for bonds in root.findall("./Forces/Force/Bonds"):
        for bond in bonds.findall("Bond"):
            bondAtoms = [int(bond.attrib[p]) for p in ("p1", "p2")]
            if shouldRemove(bondAtoms):
                bonds.remove(bond)
    for angles in root.findall("./Forces/Force/Angles"):
        for angle in angles.findall("Angle"):
            angleAtoms = [int(angle.attrib[p]) for p in ("p1", "p2", "p3")]
            if shouldRemove(angleAtoms):
                angles.remove(angle)
    for torsions in root.findall("./Forces/Force/Torsions"):
        for torsion in torsions.findall("Torsion"):
            torsionAtoms = [int(torsion.attrib[p]) for p in ("p1", "p2", "p3", "p4")]
            if shouldRemove(torsionAtoms):
                torsions.remove(torsion)

    if removeConstraints:
        for constraints in root.findall("./Constraints"):
            for constraint in constraints.findall("Constraint"):
                constraintAtoms = [int(constraint.attrib[p]) for p in ("p1", "p2")]
                if shouldRemove(constraintAtoms):
                    constraints.remove(constraint)

    return mm.XmlSerializer.deserialize(ET.tostring(root, encoding="unicode"))


def removeMMInteraction(system, ml_atoms):
    newSystem = removeBonds(system, ml_atoms)
    for force in newSystem.getForces():
        if isinstance(force, mm.NonbondedForce):
            for i in range(len(ml_atoms)):
                for j in range(i):
                    force.addException(ml_atoms[i], ml_atoms[j], 0, 1, 0, True)
        elif isinstance(force, mm.CustomNonbondedForce):
            existing = set(tuple(force.getExclusionParticles(i)) for i in range(force.getNumExclusions()))
            for i in range(len(ml_atoms)):
                a1 = ml_atoms[i]
                for j in range(i):
                    a2 = ml_atoms[j]
                    if (a1, a2) not in existing and (a2, a1) not in existing:
                        force.addExclusion(a1, a2, True)
    return newSystem

mixed_system = removeMMInteraction(system, ml_atoms)
print("number of forces (before nnp)  = ", mixed_system.getNumForces())
print(mixed_system.getForces())

#-----------------------------------------------------------------------------------------------------------------#

class HybridNNP(torch.nn.Module):
    def __init__(self, atomic_numbers, ml_atoms):
        super().__init__()
        self.indices = torch.tensor(ml_atoms, dtype=torch.long, device=device)
        self.atomic_numbers = torch.tensor(atomic_numbers, device=device).unsqueeze(0)
        self.model = ANI2x(periodic_table_index=True)
        self.model.to(self.atomic_numbers.device)
    def forward(self, positions):
        positions = positions.to(self.atomic_numbers.device)
        positions = positions[self.indices]
        positions = positions.unsqueeze(0).float() * 10 # nm -> Å
        result = self.model((self.atomic_numbers, positions))
        energy = result.energies[0] * 2625.5 # Hartree -> kJ/mol
        return energy

#------------------------------------------------------------------------------------------------#

hybrid_nnp = HybridNNP(atomic_numbers, ml_atoms)
torch_module = torch.jit.script(hybrid_nnp)
torchforce = TorchForce(torch_module)
mixed_system.addForce(torchforce)

print("number of forces (after nnp) = ", mixed_system.getNumForces())
print(mixed_system.getForces())

#-------------------------------------------------------------------------------------------------#

# Set up simulation
temperature = 500 * unit.kelvin
frictionCoeff = 1 / unit.picosecond
timeStep = 0.5 * unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, frictionCoeff, timeStep)
#
simulation = app.Simulation(pdb.topology, mixed_system, integrator)
print("Simulation is running on:", simulation.context.getPlatform().getName())
#
simulation.reporters.clear()
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
with open(f'min.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
#
state = simulation.context.getState(getEnergy=True, getForces=True)
openmm_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(openmm_energy)
openmm_force = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
print(openmm_force)
#
simulation.reporters.append(app.DCDReporter('test_ani_mixed.dcd', 100))
simulation.context.setVelocitiesToTemperature(temperature)
reporter = app.StateDataReporter(
    'data.csv', 
    100, 
    step=True, 
    time=True, 
    potentialEnergy=True,
    kineticEnergy=True,
    totalEnergy=True,  
    temperature=True,
    volume=True,
    density=True, 
    speed=True)
simulation.reporters.append(reporter)
simulation.step(400000)
