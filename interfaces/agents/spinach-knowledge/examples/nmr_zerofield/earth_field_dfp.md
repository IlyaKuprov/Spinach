# examples/nmr_zerofield/earth_field_dfp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/earth_field_dfp.m`
- Signature: `earth_field_dfp()`
- Total lines: 102

## Purpose

Earth's field NMR Simulation for 2,6-difluoropyridine; replicates simulated spectra in Figure 7 of without the weighted addition of the uncoupled 1H signal. Calculation time: seconds.

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Earth's field; implemented by `sys.magnet=2*pi*2312.45/spin('1H')`.
- Lines 18-19: Isotopes and labels; implemented by `sys.isotopes={'1H', '1H', '1H', '19F', '19F', '14N'}`.
- Lines 22-23: Chemical shifts; implemented by `inter.zeeman.scalar={6.98 8.06 6.98 -70.69 -70.69 0}`.
- Lines 25-26: J-couplings (experimental); implemented by `inter.coupling.scalar=cell(6,6)`.
- Lines 42-43: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-47: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 52-53: Sequence parameters; implemented by `parameters.sweep=1e4`.
- Lines 63-67: This needs a GPU sys.enable={'gpu'};; implemented by `R_14N=[0 10 100 500 1e3 1e4 1e5 1e6]`.
- Lines 66-67: 14N relaxation rates, Hz; implemented by `R_14N=[0 10 100 500 1e3 1e4 1e5 1e6]`.
- Lines 69-70: Get the figure going; implemented by `kfigure(); scale_figure([1.0 2.0])`.
- Lines 72-73: Loop over 14N relaxation rates; implemented by `for n=1:numel(R_14N)`.
- Lines 75-76: Set relaxation rates; implemented by `inter.r1_rates={0.22 0.22 0.22 0.22 0.22 R_14N(n)}`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 83-84: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 86-87: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 90-93: Frequency axis; implemented by `axis_hz=sweep2ticks(parameters.offset, parameters.sweep, parameters.zerofill)`.

### Control flow inferred from the code

- Line 73: `for` loop over `n=1:numel(R_14N)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=2*pi*2312.45/spin('1H')`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H', '1H', '1H', '19F', '19F', '14N'}`.
- Lines 20: computes `sys.labels` using `sys.labels= {'H3', 'H4', 'H5', 'F2', 'F6', 'N1' }`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.98 8.06 6.98 -70.69 -70.69 0}`.
- Lines 26: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6,6)`.
- Lines 27: computes `inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'N1')}` using `inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'N1')} = 37.35`.
- Lines 28: computes `inter.coupling.scalar{idxof(sys,'F6'),idxof(sys,'N1')}` using `inter.coupling.scalar{idxof(sys,'F6'),idxof(sys,'N1')} = 37.35`.
- Lines 29: computes `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'N1')}` using `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'N1')} = (1.03+0.71)/2`.
- Lines 30: computes `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'N1')}` using `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'N1')} = (1.03+0.71)/2`.
- Lines 31: computes `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F2')}` using `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F2')} = -2.47`.
- Lines 32: computes `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F6')}` using `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F6')} = -2.47`.
- Lines 33: computes `inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F2')}` using `inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F2')} = 8.08`.
- Lines 34: computes `inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F6')}` using `inter.coupling.scalar{idxof(sys,'H4'),idxof(sys,'F6')} = 8.08`.
- Lines 35: computes `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F2')}` using `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'F2')} = 1.29`.
- Lines 36: computes `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F6')}` using `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'F6')} = 1.29`.
- Lines 37: computes `inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'F6')}` using `inter.coupling.scalar{idxof(sys,'F2'),idxof(sys,'F6')} = -12.23`.
- Lines 38: computes `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'H4')}` using `inter.coupling.scalar{idxof(sys,'H3'),idxof(sys,'H4')} = 7.92`.
- Lines 39: computes `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'H4')}` using `inter.coupling.scalar{idxof(sys,'H5'),idxof(sys,'H4')} = 7.92`.

## Implementation structure

- Earth's field NMR Simulation for 2,6-difluoropyridine; replicates
- simulated spectra in Figure 7 of
- without the weighted addition of the uncoupled 1H signal.
- Calculation time: seconds.
- Earth's field
- Isotopes and labels
- Chemical shifts
- J-couplings (experimental)
- Basis set
- Relaxation theory
- Sequence parameters
- This needs a GPU

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `idxof()`, `kfigure()`, `scale_figure()`, `R_14N()`, `create()`, `basis()`, `liquid()`, `fftshift()`, `sweep2ticks()`, `xlim()`, `kxlabel()`, `set()`.
