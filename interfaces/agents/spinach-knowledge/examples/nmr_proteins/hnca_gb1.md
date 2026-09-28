# examples/nmr_proteins/hnca_gb1.m

- Signature: `hnca_gb1()`

## Purpose

Simulates an HNCA spectrum of GB1 protein, assuming that only the backbone is 13C,15N-labelled. The source notes a calculation time of minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

The example imports GB1 protein data from `2N9K.pdb` and `2N9K.bmrb`, selecting the minimal backbone with `options.select='backbone-minimal'` and deleting spins without shifts. It sets the magnetic field to 14.1 and simulates the `hnca` liquid-state NMR sequence for `15N`, `13C`, and `1H`.

## Numerical / algorithmic content

The simulation uses the `sphten-liouv` formalism with the `IK-1` approximation and scalar-coupling connectivity. It enables the `greedy` option and explicitly disables `krylov`; the commented `'gpu'` option is not enabled. The four acquired components receive squared-cosine apodisation in all three dimensions. FFTs with zero filling and `fftshift` are applied along F3, F2, and F1; conjugate combinations form the absorption components before the final spectrum is plotted.

## Implementation structure

1. Import the protein data with `protein('2N9K.pdb','2N9K.bmrb',options)`; set `sys.magnet=14.1`, `sys.tols.inter_cutoff=2.0`, and `sys.tols.prox_cutoff=4.0`.
2. Configure the basis with `bas.inter_level=4` and `bas.prox_level=1`, then build the spin system using `create` and `basis`.
3. Set spin order to `{'15N','13C','1H'}`, sweep widths to `[2800 5000 3000]`, offsets to `[-7200 8600 5100]`, acquisition points to `[128 128 128]`, zero-fill sizes to `[256 256 256]`, and axis units to `ppm`.
4. Run `fid=liquid(spin_system,@hnca,parameters,'nmr')`, apodise `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`, and assemble the three-dimensional Fourier-domain spectrum.
5. Plot `-real(spectrum)` with `plot_3d`.