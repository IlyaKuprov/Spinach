# tests/kernel/test_tensor_vector_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_tensor_vector_suite.m`
- Signature: `result=test_tensor_vector_suite()`
- Total lines: 136

## Purpose

Tests tensor, vector, distribution, and relaxation utilities. Syntax: result=test_tensor_vector_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_corr_system()`, `local_tensor_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Announce the test target; implemented by `fprintf('TESTING: Tensor, vector, and relaxation utilities\n')`.
- Lines 21-24: State the utility target of the test; implemented by `result=new_test_result('kernel/tensor_vector_suite', 'Tensor, vector, and relaxation utilities', 'Tensor and vector utility helpers must match exact small analytical cas…`.
- Lines 26-27: Check cubic Hermite interpolation on an exactly representable parabola; implemented by `spline_grid=linspace(0,1,6)`.
- Lines 34-35: Check skew normal density in the zero-skew normal-distribution limit; implemented by `pdf_grid=-2:2`.
- Lines 42-43: Check phantom-to-Fokker-Planck state embedding; implemented by `spin_state=[1;2]`.
- Lines 51-52: Check painted-image extraction from a Fokker-Planck vector; implemented by `paint_obs=fpl2phan(fpl_obs,[1;0],[2 2])`.
- Lines 57-58: Check spatial averaging of a Fokker-Planck vector back to spin space; implemented by `rho_obs=fpl2rho(fpl_obs,[2 2])`.
- Lines 64-65: Check isotropic rotational correlation coefficients and rates; implemented by `corr_system=local_corr_system()`.
- Lines 77-78: Check tensor isotropic-shift replacement without changing anisotropy; implemented by `orig_tensor=diag([1 2 6])`.
- Lines 89-90: Check direct coupling tensor extraction from both storage directions; implemented by `spin_system=local_tensor_system()`.
- Lines 98-99: Check g-tensor conversion from Zeeman interaction scaling; implemented by `spin_system.inter.zeeman.ddscal={4*eye(3)}`.
- Lines 107-108: Check isotropic Zeeman offset extraction in Hz; implemented by `spin_system.inter.zeeman.matrix={2*pi*20*eye(3)}`.

### Key state/data transformations

- Lines 22-24: computes `result` using `result=new_test_result('kernel/tensor_vector_suite', 'Tensor, vector, and relaxation utilities', 'Tensor and vector utility helpers must match exact small analytical cas…`.
- Lines 27: computes `spline_grid` using `spline_grid=linspace(0,1,6)`.
- Lines 28: computes `spline_obs` using `spline_obs=herm_spline(0,0,1,2,spline_grid)`.
- Lines 29: computes `spline_ref` using `spline_ref=spline_grid.^2`.
- Lines 35: computes `pdf_grid` using `pdf_grid=-2:2`.
- Lines 36: computes `pdf_obs` using `pdf_obs=snormpdf(pdf_grid,1,2,0)`.
- Lines 37: computes `pdf_ref` using `pdf_ref=exp(-0.5*((pdf_grid-1)/2).^2)/(2*sqrt(2*pi))`.
- Lines 39-40: computes `1e-14,1e-14, 'Azzalini skew normal with alpha` using `1e-14,1e-14, 'Azzalini skew normal with alpha=0 is the ordinary normal density')`.
- Lines 40: computes `'Azzalini skew normal with alpha` using `'Azzalini skew normal with alpha=0 is the ordinary normal density')`.
- Lines 43: computes `spin_state` using `spin_state=[1;2]`.
- Lines 44: computes `phantom` using `phantom=[1 3;2 4]`.
- Lines 45: computes `fpl_obs` using `fpl_obs=phan2fpl(phantom,spin_state)`.
- Lines 46: computes `fpl_ref` using `fpl_ref=kron(phantom(:),spin_state)`.
- Lines 52: computes `paint_obs` using `paint_obs=fpl2phan(fpl_obs,[1;0],[2 2])`.
- Lines 58: computes `rho_obs` using `rho_obs=fpl2rho(fpl_obs,[2 2])`.
- Lines 59: computes `rho_ref` using `rho_ref=mean(phantom(:))*spin_state`.
- Lines 65: computes `corr_system` using `corr_system=local_corr_system()`.
- Lines 66: computes `[weights,rates,states]` using `[weights,rates,states]=corrfun(corr_system,2,3,3,3,3)`.

### Local helper functions

- Line 117: `local_corr_system()` — `function spin_system=local_corr_system()`. Minimal spin_system for tensor extractors
  - Representative operation: `spin_system.bas.formalism='sphten-liouv'`.
  - Representative operation: `spin_system.bas.basis=[1 0;0 1;1 1;0 0]`.
- Line 127: `local_tensor_system()` — `function spin_system=local_tensor_system()`.
  - Representative operation: `spin_system.comp.nspins=2`.
  - Representative operation: `spin_system.comp.isotopes={'1H','13C'}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks Hermite splines, skew-normal density in the normal
- limit, Fokker-Planck vector reshaping helpers, correlation-function
- coefficients, tensor isotope-shift helpers, and small spin-system
- tensor extractors.

## Implementation structure

- Tests tensor, vector, distribution, and relaxation utilities. Syntax:
- result=test_tensor_vector_suite()
- result -regression test result with explanatory messages
- The test checks Hermite splines, skew-normal density in the normal
- limit, Fokker-Planck vector reshaping helpers, correlation-function
- coefficients, tensor isotope-shift helpers, and small spin-system
- tensor extractors.
- Announce the test target
- State the utility target of the test
- Check cubic Hermite interpolation on an exactly representable parabola
- Check skew normal density in the zero-skew normal-distribution limit
- Check phantom-to-Fokker-Planck state embedding

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `herm_spline()`, `test_close()`, `snormpdf()`, `phan2fpl()`, `phantom()`, `fpl2phan()`, `fpl2rho()`, `local_corr_system()`, `corrfun()`, `test_true()`, `isequal()`, `logical()`, `shift_iso()`, `local_tensor_system()`, `get_coupling()`.
