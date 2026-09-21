# examples/dnp_liq/sdnp/scalar_dnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/sdnp/scalar_dnp.m`
- Signature: `scalar_dnp()`
- Total lines: 116

## Purpose

Field dependence of the couping factor between 13C of CHCl3 and the electron spin of a nitroxide radical. Further particulars here: Experimental data from: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Field dependence of the couping factor between 13C of CHCl3 and the
- electron spin of a nitroxide radical. Further particulars here:
- Experimental data from:
- Calculation time: seconds
- Spin system
- Coordinates for dipolar Redfield
- Zeeman interactions for CSA/g-aniso Redfield
- Static hyperfine coupling
- Formalism and approximation
- Relaxation theories
- Electron R1 and R2 for empirical T1/T2
- Nuclear R1 and R2 for empirical T1/T2

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `b_vector()`, `create()`, `basis()`, `relaxation()`, `state()`, `R1n()`, `kfigure()`, `expt_data()`, `set()`, `kxlabel()`, `kylabel()`.
