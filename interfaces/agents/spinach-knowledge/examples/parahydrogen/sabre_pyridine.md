# examples/parahydrogen/sabre_pyridine.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/sabre_pyridine.m

## Purpose and experimental context

The source calls this a SABRE experiment simulation and says it is set to reproduce Figure 3b of Atkinson et al., “Spontaneous Transfer of Para hydrogen Derived Spin Order to Pyridine at Low Magnetic Field,” *Journal of the American Chemical Society* 131 (37), 13362–13368 (2009), [doi:10.1021/ja903601p](https://doi.org/10.1021/ja903601p). The paper describes reaction of [Ir(COD)(PCy3)(py)]BF4 with dihydrogen and pyridine to form fac,cis-[Ir(PCy3)(py)3(H)2]BF4. In the dihydride, para-H2-derived hydride order transfers at low field to other NMR-active nuclei; exchange of the two pyridine ligands trans to hydride with free pyridine carries the polarisation into solution. The abstract reports experimental signal enhancements above 100-fold. The MATLAB header separately credits the simulation to Eibe Duecker and Christian Griesinger.

In SABRE, para-H2 supplies a proton singlet and reversible substrate binding/exchange at the catalyst transfers spin order to free pyridine without hydrogenating the pyridine. This is distinct from ALTADENA, where hydrogenation-created spin order is transported through a field change. The script is a reduced coherent spin-dynamics model, not a simulation of the complete SABRE exchange cycle: it contains five pyridine protons and two hydride protons only, initialises the hydrides directly in a singlet, and includes no iridium spins, binding/exchange kinetics, para-H2 replenishment, or relaxation superoperator. Removing spins with `decouple` is an idealised switch in the spin Hamiltonian, not a chemical exchange model. The later low-to-high field ramp is a programmed field-cycle step, not proof of a full ALTADENA or SABRE reaction sequence.

## Spin model and sequence

All seven nuclei are `1H`. In proton order 1–7, the Zeeman scalar shifts are 8.54, 7.44, 7.86, 7.44, 8.54, −23.5, and −23.5 ppm; the first five represent pyridine and the final pair the hydrides. The nonzero scalar couplings are in Spinach's Hz convention: pyridine pairs 1–2 and 4–5, 4.88 Hz; 1–4 and 2–5, 1.00 Hz; 1–3 and 3–5, 1.84 Hz; 1–5, −0.13 Hz; 2–3 and 3–4, 7.67 Hz; and 2–4, 1.37 Hz. Hydride couplings are 1–6, 1.12 Hz; 1–7, 1.02 Hz; and 6–7, 7.00 Hz. The basis is spherical-tensor Liouville space without truncation or approximation.

The model starts with a singlet on spins 6–7 at 25 mT and evolves the coupled seven-spin system for 2.5 s. It then disconnects the hydrides and evolves the five-spin substrate for another 2.5 s. The field is raised exponentially from 25 mT to 7.05 T over 5 s in 1024 steps, followed by 1 s at high field. These are coherent propagations; no physical relaxation is specified.

## Detection and interpretation

A proton `pi/2` pulse about y is followed by a simulated FID with 1024 samples, zero-filled to 4096 and exponentially apodised with factor 6 before Fourier transformation. The code sets `offset=2400`, `sweep=600`, and `axis_units='kHz'`; the two numeric acquisition settings are not given units in the source, so they are left as source values rather than relabelled. The plotted spectrum is a calculated signal from the idealised model, not an experimental hyperpolarisation measurement; apodisation is not a relaxation model.
