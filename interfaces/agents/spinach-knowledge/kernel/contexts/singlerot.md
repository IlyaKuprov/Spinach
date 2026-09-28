# kernel/contexts/singlerot.m

- Signature: `[answer,sph_grid]=singlerot(spin_system,pulse_sequence,...`

## Purpose

Single-angle spinning context. In Liouville space, it passes a Fokker–Planck evolution generator to a user-supplied pulse-sequence function. In Hilbert space, it passes a stack of spin Hamiltonians, one for each rotor phase.

## Physical / mathematical content

- The Hamiltonian at each rotor phase combines an isotropic term with orientation-dependent terms rotated using Wigner matrices. The rotor phases number `2*parameters.max_rank+1`.
- In Liouville space, the enlarged rotor-phase and spin space includes a rotor-turning generator proportional to `2*pi*parameters.rate`, alongside relaxation and kinetics generators projected into that space.

## Numerical / algorithmic content

- The code evaluates the pulse sequence at each orientation in the spherical grid, using a parallel powder loop unless serial execution is requested.
- For Liouville-space formalisms, the pulse sequence receives the assembled Fokker–Planck generator; for Hilbert-space formalisms, it receives the Hamiltonian stack. The result is either a grid-weighted sum or a cell array of orientation-specific results, according to `parameters.sum_up`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.
- parameters.axis -spinning axis, given as a normalized
- 3-element vector; this is the direction
- around which the rotor is turning
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.max_rank -maximum rotor harmonic rank to retain
- in the solution (increase till conver-
- gence is achieved, a good guess value
- is the number of spinning sidebands
- expected in the spectrum)
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N',3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.grid -spherical grid file name; see grids
- directory in the kernel. Two-angle
- grids should be used in Liouville
- space and three-angle grids in Hil-
- bert space.
- parameters.needs -a cell array of character strings spe-
- cifying additional requirements that
- the sequence has:
- 'iso_eq' -thermal equilibrium state
- of the isotropic Hamiltonian will be
- placed into parameters.rho0
- parameters.sum_up -when set to 1 (default), returns the
- powder average. When set to 0, returns
- individual answers for each point in
- the powder as a cell array.
- parameters.* -additional subfields may be required
- by the pulse sequence -check its do-
- cumentation page
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- The main function applies defaults, validates inputs, constructs the Hamiltonian components and dissipators, and loads the spherical grid. It then prepares formalism-specific rotor-phase operators and calls the pulse sequence for each grid orientation.
- Local `defaults` and `grumble` functions set optional fields and enforce input constraints, respectively.
