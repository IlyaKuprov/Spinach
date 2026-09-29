# examples/esr_liq_pulsed/endor_phenyl.m

- Call: `endor_phenyl()` (no arguments; settings are defined in the function).
- Source: [`examples/esr_liq_pulsed/endor_phenyl.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/endor_phenyl.m)
- Experimental reference: Kasai, Hedaya & Whipple, “Electron spin resonance study of phenyl radicals isolated in an argon matrix at 4°K,” *J. Am. Chem. Soc.* **91** (1969), 4364–4368, [DOI 10.1021/ja01044a008](https://doi.org/10.1021/ja01044a008).

## Spin system

The source specifies one electron and five ring protons (two ortho, two meta, one para), with isotropic electron g-factor 2.0024 and experimental proton couplings: ortho 17.4 G, meta 5.9 G, para 1.9 G. In code these are converted to mT arguments for `mt2hz` (1.74, 0.59, and 0.19, respectively), using the phenyl g-factor. The ortho and meta pairs are equivalent and use direct-product `S2 x S2` symmetry; the para proton is not included in a symmetry pair. The source sets `sys.magnet=0.33` without annotating that field value's unit. It uses `sphten-liouv` and `approximation='none'`.

## Simulation and output

The liquid-state Mims ENDOR call is `liquid(spin_system,@endor_mims,parameters,'esr')`. Parameters set offset 0, 512 points, sweep `3e8`, `tau=100e-9` (100 ns), zero filling to 4096, spin channel `{'E'}`, and axis units MHz. The source mean-centres the FID, applies Kaiser apodisation with parameter 6, Fourier-transforms, then plots the spectrum magnitude and labels the axis “Nuclear frequency, MHz”. It opens a figure but does not save the FID or spectrum to a file. The source estimates a calculation time of seconds.

## Scope and caveats

The experimentally reported G values and their mT conversions are both retained above; do not read the code's converted inputs as if they were reported directly in G. The source does not specify units for `sys.magnet` or `parameters.sweep`, and it provides no explicit relaxation model.
