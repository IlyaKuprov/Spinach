# examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m`
- Signature: `hnco_ubiquitin_expt()`
- Total lines: 53

## Purpose

Experimental HNCO spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton)

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Data loading and truncation; implemented by `load('hnco_ubiquitin_expt.mat','fid')`.
- Lines 12-13: Apodisation; implemented by `fid=apodisation([],fid,{{'cos'},{'cos'},{'cos'}})`.
- Lines 15-16: F3 processing; implemented by `fid=fft(fid,256,1); fid=fid*exp(-1i*0.75)`.
- Lines 18-19: F2 processing; implemented by `fid=real(fid(:,1:2:end,:))+1i*real(fid(:,2:2:end,:))`.
- Lines 22-23: F1 processing; implemented by `fid=real(fid(:,:,1:2:end))-1i*real(fid(:,:,2:2:end))`.
- Lines 26-27: Window shifting; implemented by `spectrum=permute(spectrum,[3 2 1])`.
- Lines 32-33: Baseline correction; implemented by `spectrum=spectrum-spectrum(end,end,end)`.
- Lines 35-36: Water signal elimination; implemented by `spectrum(:,:,1:40)=0`.
- Lines 38-39: Parameters; implemented by `spin_system.inter.magnet=11.7395`.
- Lines 48-49: Plotting; implemented by `spin_system.sys.disable={'colorbar'}; spin_system.sys.output=1; kfigure()`.

### Key state/data transformations

- Lines 10: computes `fid` using `fid=fid(1:64,:,:)`.
- Lines 24: computes `spectrum` using `spectrum=fft(fid,256,3)`.
- Lines 36: computes `spectrum(:,:,1:40)` using `spectrum(:,:,1:40)=0`.
- Lines 39: computes `spin_system.inter.magnet` using `spin_system.inter.magnet=11.7395`.
- Lines 40: computes `spin_system.sys.disable` using `spin_system.sys.disable={'colorbar'}`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'15N','13C','1H'}`.
- Lines 42: computes `parameters.offset` using `parameters.offset=[-5870 22164 3653]`.
- Lines 43: computes `parameters.sweep` using `parameters.sweep=[2000 1500 3000]`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=[256 256 256]`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=[64 64 64]`.
- Lines 46: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Experimental HNCO spectrum of human ubiquitin.
- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Data loading and truncation
- Apodisation
- F3 processing
- F2 processing
- F1 processing
- Window shifting
- Baseline correction
- Water signal elimination
- Parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `fid()`, `apodisation()`, `fftshift()`, `spectrum()`, `kfigure()`, `plot_3d()`.
