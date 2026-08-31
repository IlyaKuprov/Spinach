# tests/kernel/test_finite_difference_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_finite_difference_suite.m`
- Signature: `result=test_finite_difference_suite()`
- Total lines: 170

## Purpose

Tests finite-difference and spectral differentiation helpers. Syntax: result=test_finite_difference_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Announce the test target; implemented by `fprintf('TESTING: Finite-difference and spectral differentiation functions\n')`.
- Lines 21-24: State the differentiation target of the test; implemented by `result=new_test_result('kernel/finite_difference_suite', 'Finite-difference and spectral differentiation functions', 'differentiation and pseudomodulation helpers must r…`.
- Lines 26-27: Three-point centred finite-difference weights at zero are exact and familiar; implemented by `w=fdweights(0,[-1 0 1],2)`.
- Lines 35-36: Five-point wall finite-difference matrix differentiates quadratics exactly on a unit grid; implemented by `x=(1:7)'; f=x.^2`.
- Lines 43-44: Savitzky-Golay differentiation recovers a cubic exactly on a uniform grid; implemented by `xg=(-5:5)'`.
- Lines 51-52: Matrix signal columns must be processed independently; implemented by `sg_mat=sgolaydiff([f 2*f],1,7,3)`.
- Lines 56-57: Savitzky-Golay differentiation requires samples down the rows; implemented by `try`.
- Lines 66-67: Savitzky-Golay differentiation requires an odd window length; implemented by `try`.
- Lines 76-77: Savitzky-Golay differentiation suppresses deterministic high-frequency noise; implemented by `grid=(0:100)'*(2*pi/100)`.
- Lines 87-88: Periodic finite-difference and Laplacian matrices annihilate constants; implemented by `Dp=fdmat(8,5,1,'pbc')`.
- Lines 99-100: Fourier spectral differentiation is exact for a represented sine wave; implemented by `N=16; [grid,D1]=fourdif(N,1); [~,D2]=fourdif(N,2)`.
- Lines 110-111: FFT differentiation kernel must reproduce the same spectral derivative; implemented by `kern=fftdiff(1,N,2*pi/N); kern=kern(:)`.
- Lines 115-116: Pseudomodulation zeroth harmonic with zero amplitude is the input spectrum; implemented by `pm_field=(0:31)'*(2*pi/32)`.
- Lines 122-123: Pseudomodulation first harmonic tends to half the modulation amplitude times the derivative; implemented by `pm_amp=1e-4`.
- Lines 129-130: Pseudomodulation second harmonic tends to minus one sixteenth of amplitude squared times the second derivative; implemented by `pm_second=pseudomodulation(pm_field,sin(2*pm_field),pm_amp,2)`.
- Lines 135-136: Pseudomodulation requires spectra down the rows; implemented by `try`.
- Lines 145-146: Pseudomodulation requires a column field axis; implemented by `try`.
- Lines 155-156: Directional derivative of a commuting matrix exponential has a closed form; implemented by `spin_system.sys.output='hush'`.

### Key state/data transformations

- Lines 22-24: computes `result` using `result=new_test_result('kernel/finite_difference_suite', 'Finite-difference and spectral differentiation functions', 'differentiation and pseudomodulation helpers must r…`.
- Lines 27: computes `w` using `w=fdweights(0,[-1 0 1],2)`.
- Lines 36: computes `x` using `x=(1:7)'; f=x.^2`.
- Lines 37: computes `D` using `D=fdmat(7,5,1,'wall')`.
- Lines 44: computes `xg` using `xg=(-5:5)'`.
- Lines 45: computes `f` using `f=2*xg.^3-xg.^2+3*xg-1`.
- Lines 46: computes `df` using `df=6*xg.^2-2*xg+3`.
- Lines 47: computes `sg_der` using `sg_der=sgolaydiff(f,1,7,3)`.
- Lines 52: computes `sg_mat` using `sg_mat=sgolaydiff([f 2*f],1,7,3)`.
- Lines 77: computes `grid` using `grid=(0:100)'*(2*pi/100)`.
- Lines 78: computes `clean` using `clean=sin(grid)`.
- Lines 79: computes `noisy` using `noisy=clean+0.03*sin(37*grid)+0.02*cos(29*grid)`.
- Lines 80: computes `raw_der` using `raw_der=fdvec(noisy,3,1)/(grid(2)-grid(1))`.
- Lines 82: computes `raw_err` using `raw_err=norm(raw_der-cos(grid),2)`.
- Lines 83: computes `sg_err` using `sg_err=norm(sg_der-cos(grid),2)`.
- Lines 88: computes `Dp` using `Dp=fdmat(8,5,1,'pbc')`.
- Lines 91: computes `Lfd` using `Lfd=fdlap([5 4],[1 2],3)`.
- Lines 94: computes `Kfd` using `Kfd=fdkup([4 4 4],[4 4 4],eye(3),3)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks finite-difference weights, finite-difference matrices,
- Fourier differentiation, Laplacians, FFT differentiation kernels,
- pseudomodulation, and matrix-exponential directional derivatives
- against exact simple cases.

## Implementation structure

- Tests finite-difference and spectral differentiation helpers. Syntax:
- result=test_finite_difference_suite()
- result -regression test result with explanatory messages
- The test checks finite-difference weights, finite-difference matrices,
- Fourier differentiation, Laplacians, FFT differentiation kernels,
- pseudomodulation, and matrix-exponential directional derivatives
- against exact simple cases.
- Announce the test target
- State the differentiation target of the test
- Three-point centred finite-difference weights at zero are exact and familiar
- Five-point wall finite-difference matrix differentiates quadratics exactly on a unit grid
- Savitzky-Golay differentiation recovers a cubic exactly on a uniform grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `fdweights()`, `test_close()`, `fdmat()`, `fdvec()`, `sgolaydiff()`, `test_true()`, `contains()`, `fdlap()`, `fdkup()`, `fourdif()`, `fourlap()`, `fftdiff()`, `kern()`, `pseudomodulation()`, `dirdiff()`.
