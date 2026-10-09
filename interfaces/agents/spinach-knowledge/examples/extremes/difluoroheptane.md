# examples/extremes/difluoroheptane.m

- Signature: `difluoroheptane()`
- Source: [`examples/extremes/difluoroheptane.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/difluoroheptane.m)

## Purpose and molecular spin model

The source describes the `19F` NMR spectrum of anti-3,4-difluoroheptane using explicit time-domain evolution in Liouville space. The MATLAB isotope list contains 23 sites: seven `12C` entries, fourteen `1H` spins and two `19F` spins. Thus the source header's “16 spins” refers to the magnetically active proton/fluorine spins; the listed carbon isotopes are spin-zero `12C` sites. The script sets `sys.magnet=11.7464` under a “Magnet induction” comment; the source does not state a unit for that assignment.

The model specifies scalar Zeeman shifts and scalar couplings. For the basis it selects `sphten-liouv`, `IK-0`, and `inter_level=1`, with a manually specified three-block basis, two `S3` symmetry groups over spins `[14 15 16]` and `[21 22 23]`, proton zero-quantum order, and projection 1. The script enables greedy parallelisation. These are the source's explicit choices; the example does not use a Hilbert-space propagator. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Acquisition and processing

This is an acquisition/FID calculation rather than an explicitly programmed RF-pulse sequence. It selects `19F` observation, sets the initial state to `state(...,'L+','19F')`, leaves decoupling empty, and passes `offset=-86700`, `sweep=300`, `npoints=512`, and `zerofill=2048` to `liquid(spin_system,@acquire,parameters,'nmr')`. The source sets `axis_units='ppm'` and `invert_axis=1`; it does not annotate units for offset or sweep. The resulting FID receives exponential apodisation with parameter 6, is Fourier transformed, and is displayed with `plot_1d`. The receiver uses the same operator description with `coil_state` instead.

## Resource note and scope

The source warns that the run needs 32 CPU cores, 128 GB of RAM, and a Titan V or later, and estimates minutes on that setup. This is the source's hardware/runtime note, not a general performance guarantee. The page preserves the 16-active-spin description and the model's stated numerical field input without adding an unsupported tesla unit. No DOI or external literature source is supplied in the script.
