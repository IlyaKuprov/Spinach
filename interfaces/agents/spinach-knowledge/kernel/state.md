# kernel/state.m

- MATLAB implementation: [kernel/state.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/state.m)

## Purpose

Create a density operator or state vector from Spinach's human-readable spin/operator descriptions. The result is a Hilbert-space density matrix, a Liouville-space state vector, or (in wavefunction formalism) a product-basis wavefunction, depending on the configured formalism and call form.

## Calls

`rho=state(spin_system,states,spins,method)`

For wavefunction formalism, supply the projection quantum number for every spin and omit the other optional arguments:

`psi=state(spin_system,mz)`

For example, `psi=state(spin_system,[-1/2 1/2 0])` specifies the projections for a `{'1H','1H','14N'}` system.

## Inputs and state construction

- In the general call, `states` names operator components and `spins` selects the target spins. A string isotope selector sums the requested one-spin operator over every matching spin, e.g. `state(spin_system,'Lz','13C')`; the selectors `'electrons'`, `'nuclei'`, and `'all'` are also documented. A numeric spin-index vector such as `[1 2 4]` likewise requests a sum over those sites.
- To form a product operator, pass matching cell arrays, e.g. `states={'Lz','L+'}; spins={1,2}`. The operators act on the listed sites as a product, not a sum; cell-array spin indices must be positive integers and cannot repeat.
- Supported operator labels are `'E'` (identity), `'Lz'`, `'Lx'`, `'Ly'`, `'L+'`, `'L-'`, `'Tl,m'` (irreducible spherical tensor, integer `l,m`), and `'CTx'`, `'CTy'`, `'CTz'`, `'CT+'`, `'CT-'` (central-transition operators in the Zeeman basis). The source routes these descriptions through `human2opspec` and constructs the corresponding basis representation.
- In `sphten-liouv`, `method` may be `'exact'` (default; correct normalisation), `'cheap'` (faster for large systems, but deliberately unnormalised), or `'chem'` (exact state weighted by concentrations in `inter.chem.concs`). The method choice is ignored by the Zeeman Hilbert and Liouville formalisms; those modes do not provide these shortcuts/chemical-kinetics weighting.

## Output and limitations

The output variable is conventionally `rho` (or `psi` in the wavefunction example); its representation is determined by the formalism. In `sphten-liouv`, a requested state must be present unambiguously in the configured basis; absent or ambiguous basis descriptors raise an error. Wavefunction mode requires one projection value per spin. Do not use the `'cheap'` result when correctly normalised state amplitudes are required.

Source and credits: [Spinach Wiki: state.m](https://spindynamics.org/wiki/index.php?title=state.m). D. Savostyanov, Luke Edwards, and Ilya Kuprov (contact details are in the source file).