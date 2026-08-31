# examples/nmr_liquids/pa_difluoroheptane_syn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_difluoroheptane_syn.m`
- Signature: `pa_difluoroheptane_syn()`
- Total lines: 127

## Purpose

Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a manual basis set specification as a merger of Lie algebras of the user-specified structral fragments followed by symmetry fac- torisation and conservation law screening. See our paper: for further information. Calculation time: minutes, faster with a GPU.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet induction; implemented by `sys.magnet=11.7464`.
- Lines 19-22: Isotopes; implemented by `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 24-25: Chemical shifts; implemented by `inter.zeeman.scalar=cell(1,23)`.
- Lines 43-44: J-couplings; implemented by `inter.coupling.scalar=cell(23)`.
- Lines 80-81: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 92-96: GPU is useful here sys.enable={'gpu'};; implemented by `sys.disable={'zte'}`.
- Lines 95-96: ZTE not useful here; implemented by `sys.disable={'zte'}`.
- Lines 98-99: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 102-103: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 114-115: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 117-118: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 120-121: Fourier transform; implemented by `spectrum=real(fftshift(fft(fid,parameters.zerofill)))`.
- Lines 123-124: Plotting; implemented by `plot_1d(spin_system,spectrum,parameters)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=11.7464`.
- Lines 20-22: computes `sys.isotopes` using `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(1,23)`.
- Lines 26: computes `inter.zeeman.scalar{14}` using `inter.zeeman.scalar{14}= 1.0189`.
- Lines 27: computes `inter.zeeman.scalar{15}` using `inter.zeeman.scalar{15}= 1.0189`.
- Lines 28: computes `inter.zeeman.scalar{16}` using `inter.zeeman.scalar{16}= 1.0189`.
- Lines 29: computes `inter.zeeman.scalar{21}` using `inter.zeeman.scalar{21}= 1.0189`.
- Lines 30: computes `inter.zeeman.scalar{22}` using `inter.zeeman.scalar{22}= 1.0189`.
- Lines 31: computes `inter.zeeman.scalar{23}` using `inter.zeeman.scalar{23}= 1.0189`.
- Lines 32: computes `inter.zeeman.scalar{11}` using `inter.zeeman.scalar{11}= 4.6138`.
- Lines 33: computes `inter.zeeman.scalar{18}` using `inter.zeeman.scalar{18}= 4.6138`.
- Lines 34: computes `inter.zeeman.scalar{17}` using `inter.zeeman.scalar{17}= 0.0000`.
- Lines 35: computes `inter.zeeman.scalar{10}` using `inter.zeeman.scalar{10}= 0.0000`.
- Lines 36: computes `inter.zeeman.scalar{8}` using `inter.zeeman.scalar{8}= 2.0878`.
- Lines 37: computes `inter.zeeman.scalar{9}` using `inter.zeeman.scalar{9}= 1.8445`.
- Lines 38: computes `inter.zeeman.scalar{13}` using `inter.zeeman.scalar{13}= 1.7132`.
- Lines 39: computes `inter.zeeman.scalar{19}` using `inter.zeeman.scalar{19}= 1.7132`.
- Lines 40: computes `inter.zeeman.scalar{20}` using `inter.zeeman.scalar{20}= 1.7018`.

## Implementation structure

- Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a
- manual basis set specification as a merger of Lie algebras of
- the user-specified structral fragments followed by symmetry fac-
- torisation and conservation law screening. See our paper:
- for further information.
- Calculation time: minutes, faster with a GPU.
- Magnet induction
- Isotopes
- Chemical shifts
- J-couplings
- Basis set
- GPU is useful here

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `false()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `plot_1d()`.
