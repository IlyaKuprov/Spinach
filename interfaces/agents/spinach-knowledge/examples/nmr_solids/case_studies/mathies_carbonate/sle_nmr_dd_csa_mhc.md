# examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m`
- Signature: `sle_nmr_dd_csa_mhc()`
- Total lines: 101

## Purpose

Water protons in the unit cell of monohydrocalcite, inc- luding slow isotropic rotational diffusion and MAS. Fur- ther details in: Calculation time: minutes, seconds with a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: 400 MHz NMR; implemented by `sys.magnet=9.4`.
- Lines 16-17: Read CASTEP file; implemented by `props=c2spinach('mhc.magres')`.
- Lines 19-20: Drop C, O, and Ca atoms; implemented by `drop_mask=ismember(props.symbols,{'C','O','Ca'})`.
- Lines 25-26: Keep two protons; implemented by `sys.isotopes{1}='1H'`.
- Lines 29-31: Convert shielding tensors into shift using the parametrisation of Huang et al. ACIE 2021; implemented by `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 34-35: Get coordinates; implemented by `inter.coordinates{1}=props.std_geom(1,:)`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-46: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 49-50: Experiment parameters; implemented by `parameters.rate=10000`.
- Lines 61-62: Isotropic rotational diffusion correlation times; implemented by `tau_c=1e-6*[0.10 1.00 10.0 100.0 1000.0]`.
- Lines 64-65: Wigner function ranks; implemented by `max_rank=[2 3 5 7 13]`.
- Lines 67-68: Get figure going; implemented by `kfigure(); hold on`.
- Lines 70-71: Loop over tau_c; implemented by `for m=1:numel(tau_c)`.
- Lines 73-74: Coorelation time and rank; implemented by `parameters.tau_c=tau_c(m)`.
- Lines 77-78: Time-domain acquisition; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 80-81: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 83-84: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.

### Control flow inferred from the code

- Line 71: `for` loop over `m=1:numel(tau_c)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 17: computes `props` using `props=c2spinach('mhc.magres')`.
- Lines 20: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'C','O','Ca'})`.
- Lines 21: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 22: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 23: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 26: computes `sys.isotopes{1}` using `sys.isotopes{1}='1H'`.
- Lines 27: computes `sys.isotopes{2}` using `sys.isotopes{2}='1H'`.
- Lines 31: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 32: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4}`.
- Lines 35: computes `inter.coordinates{1}` using `inter.coordinates{1}=props.std_geom(1,:)`.
- Lines 36: computes `inter.coordinates{2}` using `inter.coordinates{2}=props.std_geom(4,:)`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 46: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 50: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 51: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 52: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.

## Implementation structure

- Water protons in the unit cell of monohydrocalcite, inc-
- luding slow isotropic rotational diffusion and MAS. Fur-
- ther details in:
- Calculation time: minutes, seconds with a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop C, O, and Ca atoms
- Keep two protons
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Get coordinates
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `create()`, `basis()`, `state()`, `kfigure()`, `tau_c()`, `max_rank()`, `gridfree()`, `apodisation()`, `fftshift()`, `plot_1d()`, `ylim()`, `klegend()`, `kylabel()`.
