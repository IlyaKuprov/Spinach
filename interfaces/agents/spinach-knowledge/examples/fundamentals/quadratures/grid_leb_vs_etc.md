# examples/fundamentals/quadratures/grid_leb_vs_etc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/grid_leb_vs_etc.m`
- Signature: `grid_leb_vs_etc()`
- Total lines: 68

## Purpose

Heuristic vs Lebedev spherical quadrature bake-off, illus- trating the fact that, well... heuristic grids suck.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-10: Evaluate Lebedev grid; implemented by `load('../../../kernel/grids/leb_2ang_rank_29.mat', 'alphas','betas','gammas','weights')`.
- Lines 13-14: Evaluate Repulsion grid with Voronoi weights; implemented by `[alphas,betas,gammas]=repulsion(numel(alphas),3,10000)`.
- Lines 21-22: Evaluate ZCWn grid with Voronoi weights; implemented by `[alphas,betas,gammas,weights]=grid_fibon('zcwn',302)`.
- Lines 25-26: Evaluate Igloo grid with Voronoi weights; implemented by `[alphas,betas,gammas,weights]=grid_igloo(17)`.
- Lines 29-30: Evaluate Stoll grid with Voronoi weights; implemented by `[alphas,betas,gammas,weights]=grid_trian('stoll',9)`.
- Lines 33-34: Evaluate ASG grid with Voronoi weights; implemented by `[alphas,betas,gammas,weights]=grid_trian('asg',9)`.
- Lines 37-38: Evaluate SOPHE grid with Voronoi weights; implemented by `[alphas,betas,gammas,weights]=grid_trian('sophe',9)`.
- Lines 41-42: Plot the profiles; implemented by `kfigure(); hold on`.
- Lines 54-55: Residual cosmetics; implemented by `kxlabel('spherical rank')`.

### Key state/data transformations

- Lines 11: computes `perf_leb` using `perf_leb=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 14: computes `[alphas,betas,gammas]` using `[alphas,betas,gammas]=repulsion(numel(alphas),3,10000)`.
- Lines 15: computes `[~,~,~,weights]` using `[~,~,~,weights]=voronoisphere([sin(betas').*cos(gammas')`.
- Lines 18: computes `weights` using `weights=weights/(4*pi)`.
- Lines 19: computes `perf_rep` using `perf_rep=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 22: computes `[alphas,betas,gammas,weights]` using `[alphas,betas,gammas,weights]=grid_fibon('zcwn',302)`.
- Lines 23: computes `perf_zcw` using `perf_zcw=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 27: computes `perf_igl` using `perf_igl=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 31: computes `perf_esp` using `perf_esp=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 35: computes `perf_asg` using `perf_asg=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.
- Lines 39: computes `perf_sop` using `perf_sop=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm')`.

## Implementation structure

- Heuristic vs Lebedev spherical quadrature bake-off, illus-
- trating the fact that, well... heuristic grids suck.
- Evaluate Lebedev grid
- Evaluate Repulsion grid with Voronoi weights
- Evaluate ZCWn grid with Voronoi weights
- Evaluate Igloo grid with Voronoi weights
- Evaluate Stoll grid with Voronoi weights
- Evaluate ASG grid with Voronoi weights
- Evaluate SOPHE grid with Voronoi weights
- Plot the profiles
- Residual cosmetics

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `grid_test()`, `repulsion()`, `voronoisphere()`, `grid_fibon()`, `grid_igloo()`, `grid_trian()`, `kfigure()`, `ylim()`, `set()`, `xlim()`, `kxlabel()`, `kylabel()`, `klegend()`, `scale_figure()`.
