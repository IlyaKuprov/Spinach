# examples/esr_liq_pulsed/endor_benzoquinone.m

- Call: `endor_benzoquinone()` (no arguments; settings are defined in the function).
- Source: [`examples/esr_liq_pulsed/endor_benzoquinone.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/endor_benzoquinone.m)
- Experimental reference: [Figure 2, DOI 10.1002/mrc.1260280313](https://doi.org/10.1002/mrc.1260280313)

## Spin system

The liquid-state continuous-wave ENDOR example models the 2-methoxy-1,4-benzoquinone radical with one electron and six `1H` spins. It sets `sys.magnet=0.33` and an isotropic electron Zeeman scalar of 2.004577. The source does not annotate the field or the raw coupling arguments with units. Electron–proton scalar couplings are passed to `mt2hz` with values `[0.08, 0.08, 0.08, -0.059, -0.364, -0.204]`; the first three proton spins (indices 2–4) are grouped under `S3` symmetry. The basis uses `sphten-liouv`, `approximation='none'`.

## Simulation and output

It calls `liquid(spin_system,@endor_cw,parameters,'esr')` with offset 0, sweep `50e6`, 1024 points, zero filling to 4096, proton channel `{'1H'}`, axis units MHz, and derivative 1. It mean-centres the FID, applies Kaiser apodisation with parameter 20, Fourier-transforms it, and opens a figure plotting the negative spectrum magnitude with `plot_1d`. No output file is written by the example; the FID and spectrum are transient variables. The source estimates a calculation time of seconds.

## Scope and caveats

This is a CW ENDOR example, not the pulsed Mims sequence used by the neighboring methyl, nitroxide, and phenyl examples. The source supplies no explicit relaxation settings and does not save numerical data. Treat 0.33 and the coupling inputs as the code values shown: the source does not state their units in comments.
