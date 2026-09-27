# examples/parahydrogen/sabre_pyridine.m

- Signature: `sabre_pyridine()`

## Purpose

Models the pyridine SABRE experiment associated with Eibe Duecker and Christian Griesinger, aiming to reproduce Figure 3b ([http://dx.doi.org/10.1021/ja903601p](http://dx.doi.org/10.1021/ja903601p)). The source gives a calculation time of minutes.

## Physical / mathematical content

The seven-spin model comprises five pyridine protons and two hydride protons, with the source's chemical shifts and scalar couplings. It begins with a singlet on the hydrides in a 25 mT polarization field. After 2.5 s of evolution, the hydride pair is decoupled and the system evolves for another 2.5 s; the field is then raised exponentially to 7.05 T over 5 s (1024 steps), followed by 1 s of high-field evolution.

## Numerical / algorithmic content

The script propagates the density operator under the Hamiltonian during the low-field, post-decoupling, and field-ramp stages. It then applies a `pi/2` y pulse and simulates a proton FID with a 2400 kHz offset, 600 kHz sweep, 1024 points, and 4096-point zero filling. The FID is exponentially apodised (factor 6) before Fourier transformation.

## Implementation structure

The script defines proton shifts and couplings, builds the spherical-tensor Liouville basis with no approximation, constructs the singlet initial state and Hamiltonian, performs the low-field, decoupling, field-ramp, and high-field stages, then sets acquisition and plotting parameters.
