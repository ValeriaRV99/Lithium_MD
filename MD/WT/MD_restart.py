"""
Restart NVT MD (WT KEDF) from the last frame of an ASE trajectory.

Use after equilibration in MD.py: default input is md_trajectory.traj (end of equil).
For restart from production, set TRAJ_RESTART below to 'md_trajectory_prod.traj@-1'.

NVT_WT did not ship MD_restart.py; this script matches MD.py (PBE, WT, ecut=50, same PP paths).
"""
from dftpy.ions import Ions
from dftpy.field import DirectField
from dftpy.grid import DirectGrid
from dftpy.functional import Functional, TotalFunctional
from dftpy.optimization import Optimization
from dftpy.mpi import MP, sprint
from ase.io import read
from ase.md.langevin import Langevin
from dftpy.api.api4ase import DFTpyCalculator
from ase import units
from ase.io.trajectory import Trajectory
from ase.constraints import FixAtoms
import os

# Last frame of equil trajectory (or 'md_trajectory_prod.traj@-1' to continue prod)
TRAJ_RESTART = "md_trajectory.traj@-1"

mp = MP(parallel=True)
atoms = read(TRAJ_RESTART)

ions = Ions.from_ase(atoms)
z_thresh = 42
index = []
for i, atom in enumerate(atoms):
    if atom.symbol == "S" or atom.symbol == "P":
        index.append(atom.index)
    elif atom.position[2] > z_thresh:
        index.append(atom.index)

constraint = FixAtoms(indices=index)
atoms.set_constraint(constraint)

grid = DirectGrid(ions.cell, mp=mp, ecut=50, full=True)
PP_list = {
    "S": "/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/S_OEPP_PZ.UPF",
    "Li": "/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/Li_OEPP_PZ.UPF",
    "P": "/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/P_OEPP_PZ.UPF",
}

PSEUDO = Functional(type="PSEUDO", grid=grid, ions=ions, PP_list=PP_list)
XC = Functional(type="XC", name="PBE")
HARTREE = Functional(type="HARTREE")
KE = Functional(type="KEDF", name="WT")
opt_options = {"econv": 1e-8, "maxiter": 300}
rho_ini = DirectField(grid=grid)
rho_ini[:] = ions.get_ncharges() / ions.cell.volume

evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
opt = Optimization(
    EnergyEvaluator=evaluator,
    optimization_options=opt_options,
    optimization_method="CG",
)
rho = opt.optimize_rho(guess_rho=rho_ini)

calc = DFTpyCalculator(optimizer=opt, evaluator=evaluator, rho=rho)
atoms.calc = calc

temperature_K = 800
friction = 0.1
timestep = 1 * units.fs
nsteps_prod = 10000
step = 0
interval = 1


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
    if os.path.isfile("dftpy_stopfile"):
        exit()


dyn = Langevin(atoms, timestep, temperature_K=temperature_K, friction=friction)
_prod = "md_trajectory_prod.traj"
_traj_mode = "a" if os.path.isfile(_prod) else "w"
traj = Trajectory(_prod, _traj_mode, atoms)
dyn.attach(check_stop, interval=1)
dyn.attach(printenergy, interval=1)
dyn.attach(traj.write, interval=1)
dyn.run(nsteps_prod)
