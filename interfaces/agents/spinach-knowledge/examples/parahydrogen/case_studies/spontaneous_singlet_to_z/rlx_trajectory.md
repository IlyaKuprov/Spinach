# examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m

- Signature: `rlx_trajectory()`
- Source: [examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m)

## Purpose

Calculates the relaxation-driven evolution of spin order for a two-proton parahydrogen molecule coordinated to a nickel cage. The source describes a large chemical-shift anisotropy, credits work for the Gloggler group, and says the paper link is forthcoming; it also points to a Mathematica worksheet. No paper identifier is given in the source.

## Spin model and mechanism

The only explicit spins are two `1H` nuclei. The singlet initial state evolves under a Hamiltonian built from the supplied anisotropic Zeeman tensors and a secular Redfield relaxation model. In this setup, the non-identical Zeeman tensors remove the equivalence that protects the singlet from coupling to other pair-spin orders; the plotted channels show how transverse exchange order, longitudinal two-spin order, and total longitudinal order change during the modeled relaxation. Coordinates are supplied for the dipolar tensor. The source does not include explicit nickel spins, a hydrogenation reaction, or a chemical-exchange kinetic model: this is a spin-dynamics trajectory, not a measured hyperpolarisation trace or a modeled ALTADENA transport / SABRE catalyst-exchange cycle.

## Coded parameters and observables

The field is `sys.magnet=18.7893` T. The coordinate triples passed for the two protons are `(0.1399, 16.6491, 21.2341)` and `(-1.1060, 18.4823, 21.4751)`; their length unit is not stated in the source. The trace-subtracted Zeeman matrices are, for spins 1 and 2 respectively:

- `[[61.897647, 13.533565, 1.607280], [11.111659, 47.381098, -25.571782], [-2.442233, -20.906587, 28.555729]]`
- `[[39.335394, 18.002885, -24.966389], [20.627469, 55.372869, 5.145218], [-20.402191, 0.581529, 43.347552]]`

The relaxation settings are `redfield`, zero equilibrium, secular retention, and `tau_c=500e-12` s (500 ps). The basis is `sphten-liouv` with no approximation, and the code makes the `nmr` assumption. Its generator is assembled as `hamiltonian(spin_system)+1i*relaxation(spin_system)`.

The initial density operator is the two-spin singlet. The source defines `A` as the sum of the two opposite-spin raising/lowering terms (the plot labels it `L+S- + L-S+`), `B=LzSz`, and `C=Lz+Sz`; the propagation detects `[A/2, 4*B, C]`. Thus the plotted channels are transverse exchange order divided by two, four times longitudinal two-spin order, and total longitudinal order. Multichannel evolution is called with a duration of `2e-3` s and 1000 points. The plotting code separately constructs `linspace(0,2,1001)` and labels it seconds, so the page preserves both coded settings rather than reconciling that displayed axis with the propagation duration. The plotted values are `real(traj)`.

## Interpretation boundary

The result is a calculation of model spin-order trajectories under the stated Hamiltonian and relaxation assumptions. It is not an experimental observation, a simulated chemical reaction, or evidence for a quantitative hyperpolarisation yield. The source does not report convergence tests or a literature DOI.
