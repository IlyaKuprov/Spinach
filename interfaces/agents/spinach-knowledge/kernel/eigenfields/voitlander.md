# kernel/eigenfields/voitlander.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/eigenfields/voitlander.m) | [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=voitlander.m)

- Signature: `spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`

## Purpose and input data

This routine integrates field-swept EPR transition contributions over the supplied spherical triangle. `triangle` has three vertices. At each vertex, `xyz` is a unit Cartesian direction column; `tf`, `tm`, `tw`, `pd`, and `tj` are real column data indexed by transitions, and `ti` carries transition identity rows. The transition counts need not be treated as a positional identity across vertices: the integrator matches roots by the two level indices in `ti` when available, then uses field continuation for roots that lack that identity information.

The spin-system and Hamiltonian inputs supply the eigensystem calculation at newly generated directions. `Ic` and `Iz` are Hermitian matrices; `Qc` and `Qz` are cell arrays of anisotropic components; `Hmw` is a perturbation operator matching the Hamiltonian dimension (or a column vector of that length). The source describes `Iz` and `Qz` as Zeeman terms normalised to 1 Tesla. The routine forms the midpoint orientation as `[0, pi/2-elev, phi]` from each Cartesian midpoint, using MATLAB spherical-coordinate angles in radians, and calls `eigenfields` there.

## Field axis, weights, and integral

`parameters.b_axis` is a finite real field-axis vector in Tesla. It is used as supplied for the returned spectrum; this function neither constructs nor normalises the axis. The source does not require it to be increasing or uniformly spaced. The separate `parameters.window` input is a two-value real window, not an axis normalisation. Transition fields and linewidths are passed to the convolution with `b_axis`, so they must use compatible field units; the source checks finiteness/shape and positivity where applicable, but performs no unit conversion. No frequency-offset grid or conversion to Hz is created here.

For each matched transition, the triangle integrator takes the arithmetic mean of the three vertex linewidths `tw`. It forms the three vertex intensity products `tm * tj * pd`, discards non-finite products, and averages the remaining values for the transition amplitude. It then evaluates the Lorentzian-convolved triangle contribution on `b_axis` and multiplies it by the spherical triangle area from `sphtarea`. Thus the vertex factors are averaged over the corners and the triangle's spherical area supplies its integration weight; this routine does not divide by a grid-weight sum or by the full-sphere area. `tm`, `tj`, and `pd` are source-provided factors, not a grid generated here.

## Adaptive spherical subdivision

The code bisects the three spherical edges, evaluates the three midpoint eigensets, and forms four child triangles. Let `spec_dir` be the direct triangle estimate and `spec_sub` the sum over those four child estimates. The corrected estimate is `(4*spec_sub-spec_dir)/3`. If `norm(spec_dir-spec_sub,2)` exceeds the absolute tolerance `parameters.int_tol`, the routine recursively integrates the four children and sums their results; otherwise it returns the corrected estimate. The recursion therefore adapts the supplied triangle, rather than sampling a new orientation grid. The four recursive branches may be dispatched asynchronously by MATLAB parallel futures.

## Inputs, constraints, and output

- `parameters.int_tol` must be a positive real scalar; the source does not explicitly test it for finiteness. It is used in the direct-versus-subdivided 2-norm test. `parameters.mw_freq` is required as a real scalar named by the source as the resonance frequency; this file does not specify its unit. `parameters.window` must contain two real values. `parameters.pp_tol` and `parameters.tm_tol` must each be real scalars.
- `parameters.fwhm` must be a finite positive scalar. For the `zeeman-hilb` formalism, `parameters.rspt_order` is required and must be a non-negative integer or `Inf`.
- Each `triangle(n).xyz` is a finite real three-vector with unit norm; transition fields are finite real columns, and the source checks matching transition-data dimensions. The function also validates the supplied structures and Hamiltonian sizes.
- `spec` has the same dimensions as `parameters.b_axis` and is the integral over this triangle, not a full-sphere-normalised spectrum.

No random grid or random seed is used in this routine. For fixed input data, the subdivision and matching rules are source-defined; the file does not claim bitwise reproducibility across MATLAB or parallel execution environments.
