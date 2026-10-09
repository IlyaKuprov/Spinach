# kernel/contexts/gridfree.m

- Signature: `answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)`
- Source: [kernel/contexts/gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/gridfree.m)
- Wiki: [gridfree.m](https://spindynamics.org/wiki/index.php?title=gridfree.m)

## Contract

`gridfree` builds a Fokker–Planck Liouvillian for magic-angle spinning and stochastic Liouville equation (SLE) simulations, then calls the pulse-sequence handle with `(spin_system,parameters,H,R,K)`. The sequence's return value is the context output; the context itself does not prescribe its shape.

The spin subspace has dimension `spn_dim=size(H,1)`, equal to the compiled terminal offset. The SLE orientation subspace has dimension `spc_dim`, obtained from the spatial operators, and the combined Fokker–Planck dimension is `spn_dim*spc_dim`. The spatial basis is the Wigner-D-function basis used by the SLE operators, not a sampled spherical `parameters.grid` (that field is rejected here). The context reports the powder average of the pulse-sequence result. It is restricted to the Liouville formalisms `zeeman-liouv` and `sphten-liouv`. Any `parameters.rframes` field is rejected: numerical rotating-frame transformations are not supported in this SLE context. Use `singlerot()` when those corrections are required.

## Spin and orientation inputs

`assumptions` is applied before constructing the Hamiltonian. `parameters.spins` lists channel spins in channel order; the documented example is `{'1H','13C'}`. `parameters.offset` supplies the matching transmitter offsets in Hz and defaults to zero offsets when omitted. `parameters.add_terms`, when present, is a cell array of {c,A} pairs; each `c*A` is added to the isotropic spin Hamiltonian after frequency offsets.

For spinning, `parameters.rate` is in Hz: positive values denote JEOL rotation and negative values Varian/Bruker rotation. `parameters.axis` is a normalised three-component axis; the implementation normalises it before forming the rotation term, which is `2*pi*rate` times the axis-weighted SLE angular-momentum operators. `parameters.max_rank` truncates the Wigner-D ranks. Increase it until convergence; the source notes that slower spinning requires higher ranks and gives the number of spinning sidebands as an approximate starting scale.

For rotational diffusion, `parameters.tau_c` is in seconds. A scalar specifies isotropic diffusion; a symmetric positive-definite 3-by-3 correlation-time tensor specifies anisotropic diffusion, with rotational diffusion tensor `inv(6*tau_c)`. For this SLE context the correlation times belong in `parameters.tau_c`, not `inter.tau_c` (the latter is for Redfield theory).

If the sequence requests `iso_eq`, the context constructs thermal equilibrium from the isotropic Hamiltonian and replaces any supplied `parameters.rho0`. Initial and detection states are placed in the `D[0,0,0]` component before the sequence call. The context also sets `parameters.spc_dim` and `parameters.spn_dim` for the sequence.

## Example from the source documentation

`parameters.spins={'1H','13C'}` shows the channel-list form. The source specifies rate, axis, offsets, rank truncation, and correlation times by the fields above; it does not provide a complete runnable parameter set.

Additional isotropic terms are checked against `bas.offsets(end)`, the compiled spin dimension across all substances.

## State-dependent chemistry boundary

This context rejects a function handle returned by `kinetics` with `Spinach:gridfree:stateDependentKinetics`. Multi-reactant or callback-rate reaction records require a custom pulse sequence using `step`/`iserstep`, rather than static context assembly; see `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m`. Constant matrix kinetics remain supported.
