# kernel/optimcon/trapdiff.m

Signature: `[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`


This routine returns directional derivatives of the trapezium-product interval propagator `expm(-i*((HL+HR)/2 + i*dt*(sqrt(3)/12)*[HL,HR])*dt)` with respect to the left and right endpoint coefficients. Here `HL=Hd{1}+cL*Hc` and `HR=Hd{2}+cR*Hc`; `[HL,HR]` denotes the commutator. The `sqrt(3)/12` coefficient and propagator form are retained from the source documentation.

`Hd` is a two-element cell array of same-size square drift matrices, ordered left edge then right edge; `Hc` is a matching-size square control operator or superoperator. `dt` must be a positive real scalar; `cL` and `cR` must be real numeric scalars. The guards do not impose finiteness on these scalar inputs. `DL` and `DR` are matrix derivatives, each the same size as the propagator, with respect to `cL` and `cR` respectively.

The routine builds endpoint-specific directions: `H_dir_L=Hc/2+1i*dt*sqrt(3)/12*(Hc*Hd{2}-Hd{2}*Hc)` and `H_dir_R=Hc/2+1i*dt*sqrt(3)/12*(Hd{1}*Hc-Hc*Hd{1})`. It passes each direction to `dirdiff` and selects the second returned derivative entry. It has no waveform, freeze or phase-cycle mask.

The trapezium-product quadrature is attributed to Iserles and Nørsett, Corollary 3.3, [doi:10.1098/rsta.1999.0362](https://doi.org/10.1098/rsta.1999.0362). The derivative method cites Goodwin and Kuprov, Eq. 16, [doi:10.1063/1.4928978](https://doi.org/10.1063/1.4928978).

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/trapdiff.m)
[Spinach Wiki](https://spindynamics.org/wiki/index.php?title=trapdiff.m)
