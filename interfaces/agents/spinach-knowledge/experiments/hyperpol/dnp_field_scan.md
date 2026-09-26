# experiments/hyperpol/dnp_field_scan.m

- Signature: `dnp=dnp_field_scan(spin_system,parameters,H,R,K)`

## Purpose

Magnetic field scan steady-state DNP experiment. Returns the steady-state population of the user-specified state as a fun- ction of magnetic field. Syntax: dnp=dnp_field_scan(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.mw_pwr -microwave power, Hz
- parameters.mw_frq -microwave frequency offset from
- the free electron frequency at
- the reference B0 field, Hz
- parameters.fields -a vector of magnetic field off-
- sets from the reference B0 field,
- Tesla
- parameters.rho0 -equilibrium state at the reference
- B0 field
- parameters.coil -coil state vector or a horizon-
- tal stack thereof
- parameters.mw_oper -microwave irradiation operator
- parameters.ez_oper -Lz operator on the electrons
- parameters.method -'backslash' to use Matlab's
- linear equation solver, 'gmres'
- to use ILU preconditioned GMRES
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- dnp -an array of steady state expectation values for
- the states specified in parameters.coil at each
- of the fields supplied
- Note: the relaxation superoperator should NOT be thermalized
- for this type of calculation.
- Note: thermal equilibrium state and relaxation superoperator are
- assumed to be the same at all fields in the sweep -DO NOT
- USE with broad magnetic field sweep experiments.

## Implementation structure

- Magnetic field scan steady-state DNP experiment. Returns the
- steady-state population of the user-specified state as a fun-
- ction of magnetic field. Syntax:
- dnp=dnp_field_scan(spin_system,parameters,H,R,K)
- parameters.mw_pwr - microwave power, Hz
- parameters.mw_frq - microwave frequency offset from
- the free electron frequency at
- the reference B0 field, Hz
- parameters.fields - a vector of magnetic field off-
- sets from the reference B0 field,
- Tesla
- parameters.rho0 - equilibrium state at the reference
