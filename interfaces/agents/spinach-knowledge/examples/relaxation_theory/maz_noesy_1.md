# examples/relaxation_theory/maz_noesy_1.m

- MATLAB implementation: [examples/relaxation_theory/maz_noesy_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/maz_noesy_1.m)

Source: [examples/relaxation_theory/maz_noesy_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/maz_noesy_1.m)

## Purpose and spin system

Simulates a NOESY spectrum for 15N-labelled methylaziridine, including scalar relaxation of the first kind associated in the source with nitrogen-centre inversion modulating J-couplings. The example cites [the study associated with this effect](https://doi.org/10.1002/ange.201410271) and estimates minutes of calculation time. The eight spins are seven protons and one 15N nucleus at 11.75 T: H3, H2, H4, N, HN, and three methyl protons.

## Interactions and relaxation model

The source supplies vacuum-DFT shielding tensors, vacuum-DFT scalar couplings, and coordinates in angstrom. It replaces each tensor's isotropic component with source values labelled as experimental: [1.3, 1.7, 1.9, 0.0, 0.1, 1.2, 1.2, 1.2] for H3, H2, H4, N, HN, and the three methyl protons, respectively; it retains the tensor anisotropy. The scalar couplings are specified for H3–H2 (-11.6), H3–H4 (5.5), H3–HN (7.5), H2–HN (12.7), H4–each methyl proton (4.5), N–H3 (4.4), N–H4 (5.2), and N–HN (44.8). The source labels these values as vacuum DFT but does not state their units.

Redfield relaxation uses tau_c={25e-12} s (25 ps), zero equilibrium, and kite retention (rlx_keep='kite'). The second configured mechanism is SRFK, with srfk_tau_c={[1.0, 1e-3]}. Modulation depths are set to 15.0 for the H3–HN, H2–HN, and H4–HN pairs (srfk_mdepth{1,5}, {2,5}, {3,5}), encoding the selected J-modulation channels described in the source. The units of the SRFK parameter pair and modulation depths are not specified in the example. The basis is sphten-liouv with IK-2, scalar-coupling connectivity, and proximity level 3; the source also sets interaction and proximity cutoffs to 2.0 and 4.0 without giving units and disables Krylov propagation.

## NOESY preparation, acquisition, and output

The sequence uses protons, offset 500, two sweeps of 1400, mixing time 2.0 s, 256 points per dimension, and zero-fill to 1024 per dimension; the axes are labelled ppm. The initial state is proton Lz, and liquid(spin_system,@noesy,parameters,'nmr') performs the simulation. Cosine apodisation is applied to both cosine and sine FIDs; their F2 transforms are combined as a States signal, Fourier transformed along F1, and the negative real spectrum is plotted. This is the modelled 2D spectrum from the stated inputs. The source cites a published effect and uses experimental isotropic components, but it does not supply an experimental overlay, numerical agreement assessment, or rate benchmark.
