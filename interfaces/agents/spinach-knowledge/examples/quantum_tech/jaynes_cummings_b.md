# examples/quantum_tech/jaynes_cummings_b.m

- Signature: `jaynes_cummings_b()`
- Source: [`examples/quantum_tech/jaynes_cummings_b.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/jaynes_cummings_b.m)

## Purpose

A finite-dimensional time-domain Jaynes–Cummings model of one electron spin coupled to one cavity mode. The script initialises transverse spin coherence with the cavity in its vacuum state, propagates in Spinach's cavity device context, and plots spin coherence, a cavity-field quadrature, and selected cavity-level populations.

## Physical model and spin selection

The source declares `sys.isotopes={'E','C5'}`: Spinach's generic electron-spin isotope `E` and a cavity oscillator truncated to five levels. It sets `sys.magnet=0.33` T and chooses the cavity frequency from `-sys.magnet*spin('E')/(2*pi)`, so the cavity is resonant with the modeled electron transition in this Spinach frequency convention. The spin–cavity exchange is `2.828e6` Hz. This is a single-spin, single-mode Jaynes–Cummings example, not a specified vacancy, dopant, or other defect model.

The sequence selects `parameters.spins={'E'}`; its offset is `5e6` Hz, sweep is `1e8` Hz, and it requests `251` points. That is an electron-spin selection in the model, not a simulated defect-specific EPR spectrum: this source does not specify defect nuclear isotopes, hyperfine tensors, g-tensor anisotropy, orientation averaging, or a measured spectrum. It uses the `sphten-liouv` formalism with `bas.approximation='none'` and calls `device(...,'cavity')`.

## Initial state and plotted observables

The initial density operator combines the spin `Lx` component paired with cavity state `BL1` and a half-weight electron-spin `E` component paired with the same cavity state. In the source's cavity-level labelling, `BL1` is the empty mode. Detection projects onto the spin `Lx`, the cavity quadrature `(C-A)/2i`, and the populations of `BL1`, `BL2`, and `BL3`. The plotted time axis has `251` samples from 0 to 2.5 μs.

The figures therefore represent simulated spin/cavity coherence and low-lying cavity-state populations for this specified model. The source sets no independent drive-amplitude parameter, cavity linewidth, or relaxation term; it supplies the cavity-device context through `device(...,'cavity')` and the listed sequence parameters. These traces should not be read as measured defect spectra, device-fidelity results, or evidence of a particular experimental Rabi frequency. The source comment's runtime estimate is not repeated as a verified runtime claim.
