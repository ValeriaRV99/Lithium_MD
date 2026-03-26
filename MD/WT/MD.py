from dftpy.ions import Ions
from dftpy.field import DirectField
from dftpy.grid import DirectGrid
from dftpy.functional import LocalPseudo, Functional, TotalFunctional
from dftpy.formats import io
from dftpy.optimization import Optimization, OESCF
from dftpy.mpi import sprint
from dftpy.functional.pseudo.psp import PSP
from dftpy.constants import environ
import numpy as np
from dftpy.mpi import MP, sprint
from ase.io import read
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from dftpy.api.api4ase import DFTpyCalculator
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.io.trajectory import Trajectory
import time
from ase.constraints import FixAtoms
import os

mp = MP(parallel=True)
atoms = read('/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/QE/PBE_SC/Li40/rho_ks_gbrv_1.xsf')
atoms = atoms.repeat((4, 4, 1))
ions = Ions.from_ase(atoms)
z_thresh = 42
index = []
for i, atom in enumerate(atoms):
    if atom.symbol=='S' or atom.symbol=='P':
        index.append(atom.index)
    elif atom.position[2] > z_thresh:
        index.append(atom.index)
        
constraint = FixAtoms(indices=index)
atoms.set_constraint(constraint)

# ions, rho_target, _ = io.read_all('/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/QE/PBE_SC/Li40/rho_ks_gbrv_1.xsf')
grid = DirectGrid(ions.cell, mp=mp, ecut=50, full=True)
PP_list = {'S': '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/S_OEPP_PZ.UPF',
          'Li': '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/Li_OEPP_PZ.UPF',
          'P' : '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/P_OEPP_PZ.UPF'}

PSEUDO = Functional(type='PSEUDO', grid=grid, ions=ions, PP_list=PP_list)
XC = Functional(type='XC',name='PBE')
HARTREE = Functional(type='HARTREE')
KE = Functional(type='KEDF', name='WT')
opt_options = {'econv' : 1e-8, 'maxiter':300}
rho_ini = DirectField(grid=grid)
rho_ini[:] = ions.get_ncharges()/ions.cell.volume

evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
opt = Optimization(EnergyEvaluator=evaluator, optimization_options = opt_options, optimization_method = 'CG')
rho = opt.optimize_rho(guess_rho=rho_ini)

calc = DFTpyCalculator(optimizer = opt, evaluator = evaluator, rho = rho)
atoms.calc = calc
MaxwellBoltzmannDistribution(atoms, temperature_K=50)

temperature_K = 800
friction = 0.1
timestep = 1 * units.fs
nsteps_equil = 2000  # 2 ps equilibration
nsteps_prod = 10000   # 10 ps production
step = 0 
interval=1

def printenergy(a=atoms):
    global step, interval
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    sprint(
        "Step={:<8d} Epot={:.5f} Ekin={:.5f} T={:.3f} Etot={:.5f}".format(
            step, epot, ekin, ekin / (1.5 * units.kB), epot + ekin
        )
    )
    step += interval

def check_stop():
    if os.path.isfile('dftpy_stopfile'): exit()


dyn = Langevin(atoms, timestep, temperature_K=temperature_K, friction=friction)
traj = Trajectory('md_trajectory.traj', 'w', atoms)
dyn.attach(check_stop, interval=1)
dyn.attach(printenergy, interval=1)
dyn.attach(traj.write, interval=1)
dyn.run(nsteps_equil)

traj.close()
dyn.observers.clear()  # Clear previous trajectory attachment

traj = Trajectory('md_trajectory_prod.traj', 'w', atoms)
dyn.attach(check_stop, interval=1)
dyn.attach(printenergy, interval=1)
dyn.attach(traj.write, interval=1)
dyn.run(nsteps_prod)


