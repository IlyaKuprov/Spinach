# examples/nmr_solids/fslghetcor_mas_gly.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/fslghetcor_mas_gly.m`
- Signature: `fslghetcor_mas_gly()`
- Total lines: 91

## Purpose

FSLG-HETCOR of alpha-glycine powder under MAS. Calculation time: hours on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-12: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 14-15: 400 MHz spectrometer; implemented by `sys.magnet=9.4`.
- Lines 17-18: Isotropic alpha-glycine chemical shifts; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Ignore interactions below 200 Hz; implemented by `sys.tols.inter_cutoff=2*pi*200`.
- Lines 34-38: Use GPU arithmetic sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Start with Lz of 1H, detect in quadrature on 13C; implemented by `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 45-46: Experiment setup; implemented by `parameters.spins={'1H','13C'}`.
- Lines 62-63: Simulation; implemented by `fid=singlerot(spin_system,@fslghetcor,parameters,'nmr')`.
- Lines 65-66: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 69-70: F2 Fourier transform; implemented by `f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 73-74: States quadrature and F1 Fourier transform; implemented by `f1_states=real(f1_cos)+1i*real(f1_sin)`.
- Lines 77-82: Custom sweep and offset in the F1 direction: 2 pulses per block in frequency-switched Lee-Goldburg; 1/sqrt(3) is the scaling of the chemical shift in the effective Hamiltonian, it may vary in the actual experiment; sqrt(2/3) accounts for the fact that the effective field is at the magic angle.; implemented by `parameters.sweep(1)=sqrt(3)*parameters.hi_pwr/(2*parameters.nblocks*sqrt(2/3))`.
- Lines 85-86: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11-12: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 15: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 29: computes `bas.level` using `bas.level=4`.
- Lines 32: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2*pi*200`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 43: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','13C')`.
- Lines 46: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 47: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 48: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 49: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 50: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 51: computes `parameters.offset` using `parameters.offset=[2e3 1e4]`.
- Lines 52: computes `parameters.hi_pwr` using `parameters.hi_pwr=83e3`.
- Lines 53: computes `parameters.cp_pwr` using `parameters.cp_pwr=[60e3 50e3]`.

## Implementation structure

- FSLG-HETCOR of alpha-glycine powder under MAS.
- Calculation time: hours on NVidia Tesla A100,
- much longer on CPU
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Basis set
- Ignore interactions below 200 Hz
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Start with Lz of 1H, detect in quadrature on 13C

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
