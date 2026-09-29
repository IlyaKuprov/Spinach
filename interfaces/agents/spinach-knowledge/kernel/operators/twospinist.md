# kernel/operators/twospinist.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/twospinist.m
Wiki: https://spindynamics.org/wiki/index.php?title=twospinist.m

## Purpose and inputs

Build one two-spin irreducible spherical tensor component. The signature is `T=twospinist(spin_system,spin_a,spin_b,indices,type)`. `spin_a` and `spin_b` are distinct, valid one-based spin indices. `indices=[L,M]` selects rank `L=1` or `L=2`; the implemented projections are `M=-1,0,+1` for rank 1 and `M=-2,-1,0,+1,+2` for rank 2. `type` is a required character input in this function signature.

## Components and coefficients

Let `P(X,Y)` denote an individual `operator(spin_system,{'X','Y'},{spin_a spin_b},type,'csc')` term: its first label is assigned to `spin_a`, its second to `spin_b`. The source combines these terms with these coefficients and signs:

- Rank 1: `T(1,+1)=-(1/2)*(P(L+,Lz)-P(Lz,L+))`; `T(1,0)=-(1/8)*(P(L+,L-)-P(L-,L+))`; `T(1,-1)=-(1/2)*(P(L-,Lz)-P(Lz,L-))`.
- Rank 2: `T(2,+2)=+(1/2)*P(L+,L+)`; `T(2,+1)=-(1/2)*(P(Lz,L+)+P(L+,Lz))`; `T(2,0)=+sqrt(2/3)*(P(Lz,Lz)-(1/4)*(P(L+,L-)+P(L-,L+)))`; `T(2,-1)=+(1/2)*(P(Lz,L-)+P(L-,Lz))`; `T(2,-2)=+(1/2)*P(L-,L-)`.

These are the literal combination factors in this function; `operator` constructs each labelled product. Rank or projection combinations outside the listed cases raise an error.

## Formalism, basis order, and action

In Zeeman formalisms, `operator` embeds the spin-local factors into the full system by Kronecker products in spin-system order; `spin_a` and `spin_b` identify the factors receiving the first and second labels. In `sphten-liouv`, the action is represented in the supplied spherical-tensor basis, using its left/right product tables; this function does not reorder that basis. The output shape follows the selected formalism: Hilbert operators are `d`-by-`d`, where `d=prod(spin_system.comp.mults)`; Zeeman Liouville superoperators are `d^2`-by-`d^2`; and the `sphten-liouv` matrix is sized by the number of rows in `spin_system.bas.basis`.

For Liouville calculations, `type='left'` acts as `O*rho`, `type='right'` as `rho*O`, `type='comm'` as `O*rho-rho*O`, and `type='acomm'` as `O*rho+rho*O`. In Hilbert calculations `type` is ignored and the operator itself is returned. The function signature still requires `type`; it does not supply a default argument.

## Implementation notes

The consistency check rejects equal spin indices, noninteger or out-of-range spin indices, noncharacter `type`, non-real or non-two-element `indices`, unsupported ranks, and noninteger projections. For a supported rank, a switch selects the explicit expression above. This routine constructs an operator or superoperator; it is not a time propagator.
