# experiments/hyperpol/solid_effect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/solid_effect.m`
- Signature: `answer=solid_effect(spin_system,parameters)`
- Total lines: 219

## Purpose

Solid effect DNP experiment, computed using the large-scale formalism described in (http://dx.doi.org/10.1039/C2CP23233B). The system is re- stricted to one electron and one nucleus type, but the number of nuc- lei may be very large. Syntax: answer=solid_effect(spin_system,parameters)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `find()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.mw_pwr -microwave power in rad/s
- parameters.theory -level of theory. Set to 'exact' for
- the electron rotating frame calcula-
- tion or to any of the following six
- options for the average Hamiltonian
- theory calculation on top of the el-
- ectron + nuclear rotating frame:
- 'ah_first_order', 'ah_second_order',
- 'ah_third_order', 'kb_first_order',
- 'kb_second_order', 'kb_third_order'.
- See average.m function for the mea-
- ning of these options.
- parameters.nuclear_frq -nuclear Zeeman frequency in rad/s
- parameters.calc_type -set to 'time_dependence' to get the
- time dependence of the longitudinal
- magnetization and to 'steady_state'
- to get the asymptotic longitudinal
- magnetization.
- parameters.time_step -if 'time_dependence' is set in the
- calc_type parameter, sets the time
- step, seconds.
- parameters.n_steps -if 'time_dependence' is set in the
- calc_type parameter, sets the num-
- ber of time steps.

## Outputs

- answer -with the 'time_dependence' calculation type, the
- function returns the observables detected using
- the coil states specified at each point in time;
- with the 'steady_state' option specified, the
- function returns the steady state values detec-
- ted using the coil states specified.
- Note: this function generates its own Liouvillian and should be
- called directly, without a context wrapper.

## Implementation structure

- Solid effect DNP experiment, computed using the large-scale formalism
- described in (http://dx.doi.org/10.1039/C2CP23233B). The system is re-
- stricted to one electron and one nucleus type, but the number of nuc-
- lei may be very large. Syntax:
- answer=solid_effect(spin_system,parameters)
- parameters.mw_pwr -microwave power in rad/s
- parameters.theory -level of theory. Set to 'exact' for
- the electron rotating frame calcula-
- tion or to any of the following six
- options for the average Hamiltonian
- theory calculation on top of the el-
- ectron + nuclear rotating frame:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `hamiltonian()`, `assume()`, `orientation()`, `operator()`, `average()`, `isfield()`, `state()`, `cell2mat()`, `equilibrium()`, `relaxation()`, `evolution()`, `unit_state()`, `ismember()`, `cellfun()`.
