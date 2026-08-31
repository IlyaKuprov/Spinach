# etc/data_processing/autophase.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/data_processing/autophase.m`
- Signature: `[spec,cheb_coeffs]=autophase(spec,guess)`
- Total lines: 98

## Purpose

Chebyshev phase corrector for 1D NMR spectra. Views the phase profile across the spectral window as a slowly va- rying function and approximates it with a linear combi- nation of low-order Chebyshev polynomials. Syntax: [spec,cheb_coeffs]=autophase(spec,guess)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `objective()`, `apply_phases()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(spec,guess)`.
- Lines 34-35: Scale the spectrum by its standard deviation; implemented by `scaling_factor=std(spec); spec=spec/scaling_factor`.
- Lines 37-41: Shift the 4-norm from imaginary to real as much as possible; implemented by `options=optimoptions('fminunc','Display','off','GradObj','off', 'MaxIterations',Inf,'FunctionTolerance',1e-12, 'OptimalityTolerance',1e-12,'StepTolerance',1e-12, 'MaxFun…`.
- Lines 44-45: Apply the correction; implemented by `spec=apply_phases(spec,cheb_coeffs)`.
- Lines 47-48: Undo spectrum scaling; implemented by `spec=scaling_factor*spec`.

### Key state/data transformations

- Lines 35: computes `scaling_factor` using `scaling_factor=std(spec); spec=spec/scaling_factor`.
- Lines 38-41: computes `options` using `options=optimoptions('fminunc','Display','off','GradObj','off', 'MaxIterations',Inf,'FunctionTolerance',1e-12, 'OptimalityTolerance',1e-12,'StepTolerance',1e-12, 'MaxFun…`.
- Lines 42: computes `cheb_coeffs` using `cheb_coeffs=fminunc(@(x)objective(spec,x),guess,options)`.
- Lines 45: computes `spec` using `spec=apply_phases(spec,cheb_coeffs)`.

### Local helper functions

- Line 53: `objective()` — `function obj=objective(spec,phis)`. Apply the phases
  - Representative operation: `spec=apply_phases(spec,phis)`.
  - Representative operation: `obj=norm(imag(spec),4)-norm(real(spec),4)`.
- Line 64: `apply_phases()` — `function spec=apply_phases(spec,phis)`. Get Chebyshev polynomials
  - Representative operation: `cheb=ones(numel(phis),numel(spec))`.
  - Representative operation: `cheb(2,:)=linspace(-1,1,numel(spec))`.
- Line 82: `grumble()` — `function grumble(spec,guess)`. Human beings are born with different capacities. If
  - Representative operation: `if (~isnumeric(spec))||(~isvector(spec))||any(~isfinite(spec(:)))|| (std(spec(:))==0)`.
  - Representative operation: `(std(spec(:))==0)`.

## Parameters / inputs

- spec -1D NMR spectrum, a complex vector
- guess -initial guess for the Chebyshev polynomi-
- al coefficients, radians. [phi 0 0] is a
- good start, where phi is the zero-order
- phase correction guess.

## Outputs

- spec -phased NMR spectrum, a column vector
- coeffs -Chebyshev polynomial coefficients of the
- phase profile across the spectrum with
- the window treated as a [-1,1] interval

## Implementation structure

- Chebyshev phase corrector for 1D NMR spectra. Views the
- phase profile across the spectral window as a slowly va-
- rying function and approximates it with a linear combi-
- nation of low-order Chebyshev polynomials. Syntax:
- [spec,cheb_coeffs]=autophase(spec,guess)
- spec -1D NMR spectrum, a complex vector
- guess -initial guess for the Chebyshev polynomi-
- al coefficients, radians. [phi 0 0] is a
- good start, where phi is the zero-order
- phase correction guess.
- spec -phased NMR spectrum, a column vector
- coeffs -Chebyshev polynomial coefficients of the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `std()`, `optimoptions()`, `fminunc()`, `objective()`, `apply_phases()`, `cheb()`, `phi_mults()`, `spec()`, `isvector()`, `any()`, `isrow()`, `guess()`.
