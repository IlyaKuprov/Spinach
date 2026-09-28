# examples/spin_chemistry/singlet_yield_anisotropy_3.m

- Signature: `singlet_yield_anisotropy_3()`

## Purpose

Calculate singlet-yield anisotropy for a model radical-pair reaction using the Haberkorn recombination model. The source comments give a run time of minutes on an NVidia Titan V card and hours on a CPU.

## Physical / mathematical content

- Set the magnetic field to `50e-6` T and the isotopes to `{'E','E','14N','14N','1H','1H','1H'}`.
- Specify five hyperfine-coupling tensors with `mt2hz`: one between spins 1 and 3, one between spins 2 and 4, and the same tensor between spin 2 and each of spins 5, 6, and 7.
- Set the Zeeman scalars to `{2.0023 2.0025 0 0 0 0 0}`.
- Use `inter.chem.rp_theory='haberkorn'`, radical-pair electrons `[1 2]`, and recombination rates `[1e6 1e6]`.

## Numerical / algorithmic content

- Use the `sphten-liouv` formalism with `IK-0` approximation and `bas.inter_level=5`.
- Set `parameters.grid='leb_1ang_rank_63'`, `parameters.spins={'E'}`, `parameters.tol=1e-2`, `parameters.verbose=0`, and `parameters.sum_up=0`.
- GPU arithmetic is not enabled in this file: `% sys.enable={'gpu'};` is commented out.

## Implementation structure

- Create the spin system with `create(sys,inter)` and apply the basis with `basis(spin_system,bas)`.
- Run `[yields,grid]=powder(spin_system,@rydmr,parameters,'labframe')`.
- Plot `cell2mat(yields)` against `grid.betas`, labelling the axes `beta spherical angle, radians` and `singlet yield`.
