# examples/microfluidics/plain_reaction.m

## Model

A homogeneous cycloaddition concentration model with two reactants, two competing product channels, and an inert fifth solvent component; there is no spin dynamics, flow, or diffusion. Explicit reaction records assign rates 0.5 and 0.1 L/(mol*s) to products 3 (endo) and 4 (exo), respectively, retaining the original numerical channels and plot labels. For concentrations A and B, the implemented rates are `dA/dt=dB/dt=-(k1+k2)*A*B`, `dP3/dt=k1*A*B`, `dP4/dt=k2*A*B`, and `dS/dt=0`. Initial concentrations are `[0.6; 0.5; 0; 0; 18.1] mol/L`.

## Integration and output

The concentration trajectory is advanced for 20 seconds in 200 steps with `step` and the `LG4` integrator. A concentration-versus-time plot shows components 1–4 in mol/L and omits solvent. The file header describes runtime as seconds; that is a source estimate, not a timing measurement here.

The five species are spin-free unit-coordinate blocks. A ghost seed is created and traced out with `kill_spin` because `create` requires a non-empty isotope input; no spin physics is introduced. This concentration-only construction does not change the three solvent protons in the shared `dac_reaction` molecular definition. `kinetics` compiles the records, `unit_state` supplies initial populations, and the existing LG4 stepping uses the concentration-dependent generator. Additive product unit arrival is shared equally between reactants instead of being assigned solely to B, so finite-step concentration histories need not equal the old asymmetric allocation.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_reaction.m
