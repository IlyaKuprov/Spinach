# tests/kernel/test_dynamic_voitlander.m

## Purpose

Regression test for `voitlander()`, the adaptive Voitlander spherical-triangle orientation integration kernel. The test verifies that a constant isotropic transition integrates to the full spherical-triangle area times a Lorentzian line, that per-vertex moment–Jacobian products are averaged correctly, that singular Jacobian vertices are rejected safely, and that local branch renumbering does not affect the result.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_voitlander.m

## Behaviour

- Announces the test target with `TESTING: Voitlander spherical-triangle integration`.
- Builds an isotropic spin-half electron (`sys.magnet=1`, `sys.isotopes={'E'}`, Zeeman matrix `2.0023*eye(3)`, `bas.formalism='zeeman-hilb'`, `bas.approximation='none'`) via `test_spin_system`.
- Sets field-swept EPR parameters: `mw_freq=9.5e9`, `fwhm=2e-3`, `window=[0.33 0.35]`, `npoints=9`, `tm_tol=0`, `rspt_order=Inf`, `int_tol=1e9`, `pp_tol=(window(2)-window(1))/(2*(npoints-1))`, `orientation=[0 0 0]`, `rho0=-state(spin_system,'Lz','E')`.
- Obtains Zeeman, coupling, and microwave Hamiltonians with `hamiltonian(assume(...))`, symmetrises `Ic` and `Iz`, and builds `Hmw=(L+ + L-)/2`.
- Finds the isotropic transition with `eigenfields` at orientation `[0 0 0]` and checks:
  - `isscalar(tf)` — one allowed EPR transition in the window.
  - `tj` equals `spin_system.tols.freeg/2.0023` within `1e-10`/`1e-12` — the field-sweep Jacobian reduces to the free-g over effective-g ratio.
- Defines the positive-octant spherical triangle with vertices `r1=[1;0;0]`, `r2=[0;1;0]`, `r3=[0;0;1]`, computes first-level subdivision midpoints with `sphtrsubd` and sub-triangle areas `area_a`–`area_d` with `sphtarea`.
- Sets `parameters.b_axis=tf+linspace(-2*tw,2*tw,npoints)` to keep all field points inside the internal six-width support window.
- Packages the triangle as a struct array with identical `tf`, `tm`, `tw`, `pd`, `ti`, `tj` at all three vertices and integrates with `voitlander`.
- Compares the result to the analytic reference `sphtarea(r1,r2,r3)*pd*line_shape`, where `line_width=tw/2` and `line_shape=(tm*tj/(pi*line_width))./(1+((b_axis-tf)/line_width).^2)`; tolerance `1e-10`/`1e-12`.
- Checks the spectrum is finite, real, and non-negative (`abs(imag)<1e-14`, `isfinite(real)`, `real>=0`).
- Vertex product averaging check: perturbs vertex `tm`/`tj` to `(tm,4*tj)`, `(2*tm,2*tj)`, `(4*tm,tj)` so per-vertex products remain `tm*tj`; the expected amplitude uses `prod_area=(4*(2*area_a+2*area_b+2*area_c+area_d)-4*sphtarea(r1,r2,r3))/3`, verifying that the Simpson–Richardson correction averages finite per-vertex products directly.
- Singular Jacobian rejection check: sets vertex Jacobians to `Inf`, `tj`, `2*tj`; the expected amplitude uses `sing_area=(4*(area_a+area_b+(4/3)*area_c+area_d)-(3/2)*sphtarea(r1,r2,r3))/3`, verifying that non-finite per-vertex Jacobian products are rejected before the Simpson–Richardson correction and do not contaminate the integral.
- Builds a minimal two-level avoided-crossing system (`spin_system.sys.output='hush'`, `bas.basis=zeros(2,1)`, `comp.mults=2`, `Iz=diag([-1 1])`, `Ic=[0 0.1;0.1 0]`, `Hmw=[0 1;1 0]`, `rho0=diag([1 0])`, `mw_freq=1/(2*pi)`, `window=[-0.8 0.8]`, `fwhm=0.01`, `npoints=9`) with two resonance roots from `eigenfields`.
- Compares integration of the high root (`high_root=2`) using stable branch labels `ti_ref=ti(high_root,:)` against locally renumbered labels `ti_local=[ti(high_root,1:2) 1]` at the other two vertices; the two spectra must agree within `1e-12`/`1e-12`, verifying that field-continuation matching ignores local branch ordinal changes at window edges.

## Inputs and outputs

```matlab
result = test_dynamic_voitlander()
```

- **Output**: `result` — regression test result structure with explanatory messages, created via `new_test_result('kernel/dynamic_voitlander', ...)` and accumulated through `test_true` and `test_close` checks.
- **Input**: none.

## References

- `voitlander` — adaptive Voitlander spherical-triangle integration kernel under test.
- `eigenfields` — resonance root and transition property extraction.
- `sphtrsubd`, `sphtarea` — spherical triangle subdivision and area computation.
- `hamiltonian`, `assume`, `state`, `orientation`, `test_spin_system`, `new_test_result`, `test_true`, `test_close` — Spinach kernel utilities used by the test.
