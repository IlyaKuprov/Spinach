# kernel/utilities/rotor_stack.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rotor_stack.m>

## Purpose

Returns a rotor stack of Liouvillians or Hamiltonians for the traditional style calculation of MAS (magic angle spinning) dynamics.

## Behaviour

- Syntax: `[L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)`.
- Calls `grumble` to enforce consistency of the inputs, then applies the requested assumptions throughout the rotor and frame pipeline via `assume`.
- Obtains the Hamiltonian (and its spherical-tensor interaction blocks `Q`) with `hamiltonian`, then applies transmitter offsets with `frqoffset`.
- Derives the rotor axis orientation from `parameters.axis` via `cart2sph`, converting the polar angle with `rotor_theta=pi/2-rotor_theta`.
- Computes rotor phases with `fourdif(spin_system,2*parameters.max_rank+1,1)`, yielding `2*max_rank+1` rotor ticks.
- For each rotor tick (parallelised with `parfor`), and for each spherical rank `r` in `Q`, builds Wigner rotations:
  - `D_rot2lab=wigner(r,+rotor_phi,+rotor_theta,0)` and `D_lab2rot=wigner(r,0,-rotor_theta,-rotor_phi)` for the rotor axis tilt.
  - `D_initial=wigner(r,orientation(1),orientation(2),orientation(3))` for the initial crystallite orientation at rotor phase zero.
  - `D_rotor=wigner(r,0,0,rotor_phases(n))` for the rotor rotation.
- Composes the rotations depending on `parameters.masframe`:
  - `'magnet'` (initial orientation in the lab frame; three-angle powder grids required): `D=D_rot2lab*D_rotor*D_lab2rot*D_initial`.
  - `'rotor'` (initial orientation in the rotor frame; two-angle powder grids required): `D=D_rot2lab*D_rotor*D_initial`.
  - Any other value raises the error `'unknown MAS frame.'`.
- Accumulates each block as `L{n}=L{n}+D(k,m)*Q{r}{k,m}` over `k` and `m` from 1 to `2*r+1`.
- Applies interaction representations: for each entry of `parameters.rframes`, calls `rotframe(spin_system,C{k},(L{n}+L{n}')/2,parameters.rframes{k}{1},parameters.rframes{k}{2})`, where `C{k}` is the carrier operator obtained from `carrier(spin_system,parameters.rframes{n}{1})`.
- Cleans up each block with `clean_up(spin_system,L{n},spin_system.tols.liouv_zero)`.
- Relaxation and chemical kinetics are not included.
- Validation (`grumble`) enforces:
  - `parameters.axis` must be present and a row vector of three real numbers.
  - `parameters.offset` must be specified, non-empty, numeric, and have the same number of elements as `parameters.spins`.
  - `parameters.spins` must be a non-empty cell array of strings, all referring to isotopes present in the system.
  - `parameters.max_rank` must be a non-negative real integer.
  - `parameters.orientation` must be a row vector of three real numbers.
  - `assumptions` must be a character string.
  - `parameters.rframes` must be a cell array whose elements are cell arrays of exactly two sub-elements: a character string naming an isotope present in the system, and a number.
  - `parameters.masframe` must be a character string equal to `'rotor'` or `'magnet'`.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object.
- `parameters` — structure with subfields:
  - `axis` — spinning axis, a normalised 3-element vector.
  - `offset` — nonempty numeric array of transmitter offsets in Hz, with one element per spin in `parameters.spins`; a cell array is rejected.
  - `spins` — cell array of the spins the offsets refer to, e.g. `{'1H','13C'}`.
  - `max_rank` — maximum harmonic rank to retain in the solution (increase until convergence is achieved; approximately equal to the number of spinning sidebands in the spectrum).
  - `rframes` — rotating frame specification, e.g. `{{'13C',2},{'14N',3}}` requests second order rotating frame transformation with respect to carbon-13 and third order with respect to nitrogen-14. When this option is used, the assumptions on the respective spins should be laboratory frame.
  - `orientation` — orientation of the spin system at rotor phase zero, a vector of three Euler angles in radians.
  - `masframe` — the frame in which the rotations are applied: `'magnet'` (initial orientation in the lab frame; three-angle powder grids required) or `'rotor'` (initial orientation in the rotor frame; two-angle powder grids required).
- `assumptions` — assumption set used in generating the Hamiltonian and validating numerical rotating frames, regardless of prior assumptions on the input object. The spins in `parameters.rframes` must remain in the laboratory frame under this set; already-rotating spins are rejected by `rotframe.m`. Numerical frames for carrier-free `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0` components are not implemented. See `assume.m`.

Outputs:

- `L` — cell array of Hamiltonian or Liouvillian matrices, one for each tick of the rotor (`2*max_rank+1` elements).
- `rotor_phases` — rotor phases at each tick, in radians.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=rotor_stack.m>
