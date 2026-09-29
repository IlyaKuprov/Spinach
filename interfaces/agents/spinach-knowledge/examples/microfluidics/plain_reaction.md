# examples/microfluidics/plain_reaction.m

## Model

A homogeneous cycloaddition concentration model with two reactants, two competing product channels, and an inert fifth solvent component; there is no spin dynamics, flow, or diffusion. The source comments assign `k1=0.5` to the exo channel and `k2=0.1` to the endo channel. The source comment labels the constants `mol/(L*s)`, but with concentrations in mol/L the bilinear terms `k*A*B` require `L/(mol*s)` for both constants. For concentrations A and B, the implemented rates are `dA/dt=dB/dt=-(k1+k2)*A*B`, `dP3/dt=k1*A*B`, `dP4/dt=k2*A*B`, and `dS/dt=0`. Initial concentrations are `[0.6; 0.5; 0; 0; 18.1] mol/L`.

## Integration and output

The concentration trajectory is advanced for 20 seconds in 200 steps with `step` and the `LG4` integrator. A concentration-versus-time plot shows components 1–4 in mol/L and omits solvent. The file header describes runtime as seconds; that is a source estimate, not a timing measurement here.

The labels need care: this file's plot legend calls components 3 and 4 endo and exo, respectively, while the `k1`/`k2` comments associate those channels with exo and endo in the opposite order. The paired flow example titles components 3 and 4 exo and endo.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_reaction.m
