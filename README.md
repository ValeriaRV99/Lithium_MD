# Lithium_MD

Orbital-free DFT molecular dynamics of a lithium thiophosphate (Li–P–S) interface,
run with [DFTpy](https://gitlab.com/pavanello-research-group/dftpy) driven through ASE,
alongside the Kohn–Sham reference calculation that supplies the starting density.

## Layout

```
MD/WT/          OF-DFT molecular dynamics with the Wang-Teter KEDF
  MD.py         main driver (DFTpy + ASE)
  MD_restart.py restart from an existing trajectory
  run.sh        SLURM submission script
QE/             Quantum ESPRESSO reference calculation
  DEN/          Kohn-Sham density used to build the OF-DFT starting grid
  interface.vasp  interface structure
```

## Setup

As configured in `MD/WT/MD.py`:

| | |
|---|---|
| Kinetic energy functional | Wang–Teter (`KEDF`, `name='WT'`) |
| Exchange–correlation | PBE |
| Pseudopotentials | OEPP (PZ) for Li, S and P |
| Plane-wave cutoff | `ecut = 50` |
| Density optimization | conjugate gradient, `econv = 1e-8`, max 300 iterations |
| Structure | QE density (`QE/DEN/rho_ks_gbrv_1.xsf`), repeated 4×4×1 |
| Constraints | all S and P atoms fixed, plus any atom above z = 42 |
| Dynamics | ASE Langevin, initialized at 50 K and run at 800 K |

## Running

The Kohn–Sham reference first, then the OF-DFT dynamics:

```bash
cd QE/DEN && sbatch run_1.sh
```

```bash
cd MD/WT && sbatch run.sh
```

Both `run.sh` scripts target a SLURM cluster (Rutgers Amarel, `price-pi` partition) and
load `intel/19.1.1` + `mvapich2/2.2`. Adjust the module and partition lines for your system.

## Note on reproducibility

`MD.py` reads its pseudopotentials from `../OF/PP/{Li,S,P}_OEPP_PZ.UPF`, which are **not
included in this repository**. You will need to supply those files, or repoint `PP_list`,
before the driver will run.

## Related

- Wang–Teter kernels: [ValeriaRV99/wt](https://github.com/ValeriaRV99/wt) and
  V. Rios-Vargas, X. Shao, S. B. Trickey, M. Pavanello, *Phys. Rev. B* **110**, 085129 (2024),
  [doi:10.1103/PhysRevB.110.085129](https://doi.org/10.1103/PhysRevB.110.085129)
- Lithium pseudopotentials: [ValeriaRV99/Li](https://github.com/ValeriaRV99/Li)
