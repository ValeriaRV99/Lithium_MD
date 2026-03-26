  
import numpy as np
from ase.build import bulk
from ase.io import read,write
from ase.io.trajectory import Trajectory
import py3Dmol
from qepy.calculator import QEpyCalculator
from ase.build import bulk
from qepy.driver import Driver
from qepy.io import QEInput
import time

from scipy.optimize import minimize
try:
    from mpi4py import MPI
    comm=MPI.COMM_WORLD
except:
    comm=None
    
qe_options = {
    '&control': {
        'calculation': "'scf'",
        'prefix': "'tmp_1'",
        'pseudo_dir': "'/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/'"},
    '&system': {
        'ibrav' : 0,
        'degauss': 0.02,
        'input_dft': "'pbe'",
        'ecutwfc': 80,
        'occupations': "'smearing'",
        'smearing': "'fermi-dirac'" 
    },
    '&electrons': {
        'conv_thr' : 1.0e-8},
     'atomic_species': ['Li 6.94 Li_OEPP_PZ.UPF', 'S  32.06  S_OEPP_PZ.UPF', 'P  30.974  P_OEPP_PZ.UPF'],
      'k_points gamma': "True",
}

atoms = read("/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/QE/interface.vasp")
#atoms *= (1, 1, 1)

qe_options = QEInput.update_atoms(atoms, qe_options = qe_options,  extrapolation=False)
#QEInput().write_qe_input('qe_1.in', qe_options=qe_options)
start_time = time.time()
print('Before SCF')
driver = Driver(qe_options=qe_options, atoms=atoms, logfile='tmp_1.out', comm=comm)
print('SCF started')
#driver.calc_energy()
driver.scf()
density = driver.get_density()
rho_ks = driver.data2field(density)
print(rho_ks.grid.nr)
ions = driver.get_dftpy_ions()
end_time = time.time()
execution_time = end_time - start_time
np.save('ks_pbe_time_1', execution_time)
rho_ks.write('rho_ks_gbrv_1.xsf', ions=ions)

