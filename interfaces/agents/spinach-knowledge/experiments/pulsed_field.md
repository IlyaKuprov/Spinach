# experiments/pulsed_field.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pulsed_field.m`
- Signature: `answer=pulsed_field(spin_system,parameters,H,R,K) %#ok<INUSD>`
- Total lines: 225

## Purpose

Magnetisation dynamics under a time-dependent magnetic field along the Z axis of the laboratory frame with spin-phonon relaxation, as measured in pulsed-field magnetometry of molecular magnets. The field profile is replaced by a staircase; on each stair the Hamiltonian is constant, the spin-phonon dissipator is rebuilt in the eigenbasis of that Hamiltonian, and the density matrix is propagated in that eigenbasis by a symmetric split: exact coherent phases for half a stair, the dissipative step to second order in the dissipator times the stair width, and the phases again. The dissipator times the stair width must be small; the coherent part is treated exactly for any stair width. The dissipator is applied as Hilbert space matrix products (see phonon_oper.m), so the cost of a stair is cubic in the dimension of the Hilbert space.

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.field_prof -function handle returning the field
- in Tesla at a time in seconds
- parameters.hzeeman -Zeeman operator per Tesla, rad/s/T,
- Hilbert space, supplied by the con-
- text when 'zeeman_op' is requested
- in parameters.needs
- parameters.timestep -stair width, seconds
- parameters.nsteps -number of stairs
- parameters.coil -Hilbert space observable operator or
- a cell array of them
- parameters.phonon_x -spin-phonon coupling operator, see
- rlx_phonon.m
- parameters.phonon_i0 -phonon spectral density prefactor,
- see rlx_phonon.m
- parameters.phonon_alpha -phonon spectral density exponent,
- 1 or above, see rlx_phonon.m
- parameters.nout -number of stairs between recorded
- observable values
- H -Hamiltonian received from the context function, Hilbert
- space, containing the Zeeman term at sys.magnet=1 Tesla;
- the function removes that term and adds the field of
- each stair itself
- R -relaxation superoperator received from the context
- function; ignored, the spin-phonon dissipator is built
- here at every stair
- K -kinetics superoperator received from the context
- function; ignored

## Outputs

- answer.t -column of recording times, seconds
- answer.field -column of field values at those times, Tesla
- answer.obs -matrix of observable expectation values, one
- column per coil, at the recording times
- Note: the sequence works in zeeman-hilb formalism under the crystal
- and powder contexts, which assemble the anisotropic part of
- the Hamiltonian; the liquid context drops that part, and with
- it the crystal field of a giant spin. The context must be
- called with the labframe assumption set, so that H and the
- Zeeman operator are built consistently; the powder context
- must be called with parameters.sum_up=false because the
- answer is a structure; additional rotating frames (parame-
- ters.rframes) and frequency offsets (parameters.offset) are
- not supported because the field operator is added in the
- laboratory frame. The temperature of the phonon bath is
- inter.temperature.
- Note: sys.magnet must be 1 Tesla, so that parameters.hzeeman is
- the Zeeman operator per Tesla; the Hamiltonian received from
- the context then contains the Zeeman term at 1 Tesla, which
- this function removes before adding the field on each stair.
- The initial state is the thermal equilibrium of the field-
- free Hamiltonian at the temperature of the phonon bath.
- Note: the Hamiltonian on each stair uses the field at the midpoint
- of the stair; answer.field is the profile evaluated at the
- recording times, which are the ends of the recorded stairs.

## Implementation structure

- Magnetisation dynamics under a time-dependent magnetic field along
- the Z axis of the laboratory frame with spin-phonon relaxation, as
- measured in pulsed-field magnetometry of molecular magnets. The
- field profile is replaced by a staircase; on each stair the Hamil-
- tonian is constant, the spin-phonon dissipator is rebuilt in the
- eigenbasis of that Hamiltonian, and the density matrix is propa-
- gated in that eigenbasis by a symmetric split: exact coherent
- phases for half a stair, the dissipative step to second order in
- the dissipator times the stair width, and the phases again. The
- dissipator times the stair width must be small; the coherent part
- is treated exactly for any stair width. The dissipator is applied
- as Hilbert space matrix products (see phonon_oper.m), so the cost

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `grumble()`, `phonon_oper()`, `report()`, `eig()`, `trace()`, `num2str()`. The dissipator is applied through the anonymous function `dissip`, rebuilt on every stair from the dressed operator in the eigenbasis.
