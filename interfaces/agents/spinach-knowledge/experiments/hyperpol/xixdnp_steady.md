# experiments/hyperpol/xixdnp_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/xixdnp_steady.m`
- Signature: `dnp=xixdnp_steady(spin_system,parameters,H,R,K)`
- Total lines: 173

## Purpose

TPPM DNP and its special case X-inverse-X (XiX) DNP experiment from (https://doi.org/10.1021/jacs.1c09900), steady state ver- sion. Call from powder context. Syntax: dnp=xixdnp_steady(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.coil -detection state vector
- parameters.pulse_dur -pulse duration, seconds
- parameters.phase -phase of each second pulse in radians
- parameters.nloops -number of XiX/TPPM DNP blocks, the
- calculation is faster when this is
- an integer power of 2
- parameters.shot_spacing -delay between microwave irradiation
- periods, seconds
- parameters.addshift -shift of the centre of the field
- profile, Hz
- parameters.el_offs -microwave resonance offsets, a vector
- of frequencies in Hz
- Output:
- dnp -steady state observable on the de-
- tection state vector as a function
- of microwave resonance offset

## Implementation structure

- TPPM DNP and its special case X-inverse-X (XiX) DNP experiment
- from (https://doi.org/10.1021/jacs.1c09900), steady state ver-
- sion. Call from powder context. Syntax:
- dnp=xixdnp_steady(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `propagator()`, `clean_up()`, `ismember()`, `gpuArray()`, `ppower()`, `gather()`, `steady()`, `dnp()`, `ismatrix()`, `isfield()`, `isscalar()`.
