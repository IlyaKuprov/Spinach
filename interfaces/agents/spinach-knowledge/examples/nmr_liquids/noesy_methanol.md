# examples/nmr_liquids/noesy_methanol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noesy_methanol.m`
- Signature: `noesy_methanol()`
- Total lines: 88

## Purpose

NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker, CSA tensors from DFT. Note the presence of cross-peaks between 13C doublet components. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-14: Spin system properties (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/methanol.log'), {{'H','1H'},{'C','13C'}},[31.8 182.4],[])`.
- Lines 16-17: Remove the OH proton; implemented by `sys.isotopes(end)=[]`.
- Lines 21-22: Put all chemical shifts on resonance; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2 3 4],[0 0 0 0])`.
- Lines 24-25: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 27-28: Assign J-couplings; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 36-37: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 47-48: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 51-52: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 55-56: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 65-66: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 68-69: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 72-73: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 76-77: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 79-80: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 82-83: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13-14: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/methanol.log'), {{'H','1H'},{'C','13C'}},[31.8 182.4],[])`.
- Lines 17: computes `sys.isotopes(end)` using `sys.isotopes(end)=[]`.
- Lines 18: computes `inter.coordinates(end)` using `inter.coordinates(end)=[]`.
- Lines 19: computes `inter.zeeman.matrix(end)` using `inter.zeeman.matrix(end)=[]`.
- Lines 22: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2 3 4],[0 0 0 0])`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 28: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 29: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=141`.
- Lines 30: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=141`.
- Lines 31: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=141`.
- Lines 32: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-11`.
- Lines 33: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=-11`.
- Lines 34: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=-11`.
- Lines 37: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 38: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 41: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 42: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 43: computes `inter.temperature` using `inter.temperature=298`.

## Implementation structure

- NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker,
- CSA tensors from DFT. Note the presence of cross-peaks between 13C
- doublet components.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Remove the OH proton
- Put all chemical shifts on resonance
- Magnet field
- Assign J-couplings
- Basis set
- Relaxation theory parameters
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
