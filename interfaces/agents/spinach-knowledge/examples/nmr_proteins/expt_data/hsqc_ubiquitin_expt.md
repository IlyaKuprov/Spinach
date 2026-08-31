# examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m`
- Signature: `hsqc_ubiquitin_expt()`
- Total lines: 43

## Purpose

Experimental HSQC spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton)

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Load the file; implemented by `load('hsqc_ubiquitin_expt.mat','fid')`.
- Lines 17-18: Phasing; implemented by `phase_f1=exp(-1i*0.7)`.
- Lines 20-21: Apodisation; implemented by `fid.pos=apodisation([],fid.pos,{{'cos'},{'cos'}})*phase_f1`.
- Lines 24-25: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 28-29: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 31-32: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 34-35: Plotting; implemented by `spectrum=flip(flip(spectrum,1),2)`.

### Key state/data transformations

- Lines 10: computes `spin_system.inter.magnet` using `spin_system.inter.magnet=11.7395`.
- Lines 11: computes `parameters.zerofill` using `parameters.zerofill=[1024 1024]`.
- Lines 12: computes `parameters.sweep` using `parameters.sweep=[2000 4000]`.
- Lines 13: computes `parameters.offset` using `parameters.offset=[-5870 3753]`.
- Lines 14: computes `parameters.spins` using `parameters.spins={'15N','1H'}`.
- Lines 15: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 18: computes `phase_f1` using `phase_f1=exp(-1i*0.7)`.
- Lines 21: computes `fid.pos` using `fid.pos=apodisation([],fid.pos,{{'cos'},{'cos'}})*phase_f1`.
- Lines 22: computes `fid.neg` using `fid.neg=apodisation([],fid.neg,{{'cos'},{'cos'}})*phase_f1`.
- Lines 25: computes `f1_pos` using `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 26: computes `f1_neg` using `f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1)`.
- Lines 29: computes `fid` using `fid=f1_pos+conj(f1_neg)`.
- Lines 32: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 36: computes `spin_system.sys.disable` using `spin_system.sys.disable={'colorbar'}`.
- Lines 37: computes `spin_system.sys.output` using `spin_system.sys.output=1`.

## Implementation structure

- Experimental HSQC spectrum of human ubiquitin.
- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Load the file
- Phasing
- Apodisation
- F2 Fourier transform
- Form States signal
- F1 Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `apodisation()`, `fftshift()`, `conj()`, `flip()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
