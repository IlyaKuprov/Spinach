# examples/fundamentals/correlation_function_5.m

- MATLAB implementation: [examples/fundamentals/correlation_function_5.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/correlation_function_5.m)

- Signature: `correlation_function_5()`

## Purpose

Estimates the rotational correlation `G(k,m,p,q)=<R(k,m)*R(p,q)>` by Monte Carlo sampling, with `R` the accumulated 3D Cartesian rotation matrix. The source describes a run time of minutes. This example plots the sampled estimate; it does not compare it with a separate analytical curve.

## Sampling and convention

The script sets `sigma_iso=0.2`, selects `(k,m,p,q)=(2,3,2,3)`, and uses `npoints=1e6` with `nlags=300`. It draws a `3×npoints` array of standard-normal angle components and starts with `R(:,:,1)=eye(3)`. For `n=2,...,npoints`, it forms a skew-symmetric generator from the three source matrices

```matlab
J1 = [ 0  1  0; -1  0  0;  0  0  0];
J2 = [ 0  0  1;  0  0  0; -1  0  0];
J3 = [ 0  0  0;  0  0  1;  0 -1  0];
R_gen = angles(1,n)*J1 + angles(2,n)*J2 + angles(3,n)*J3;
```

The accumulated rotation is updated on the right as `R(:,:,n)=R(:,:,n-1)*expm(sigma_iso*R_gen)`. Thus the first angle column is generated but not used in an update.

## Correlation estimate and output

The selected matrix-element traces are passed to `xcorr(...,nlags,'normalized')`. The returned correlation and lag arrays are each transformed with `ifftshift`; the plot uses entries `1:nlags` of those shifted arrays, plotting `real(cf_mc)` against lag in points. The computed correlation is multiplied by `1/3` before plotting. Thus the displayed result is the real part of this finite Monte Carlo estimate under the source's normalisation and index selection.

## Scope

The angle sequence is random and the script sets no seed. It supplies no uncertainty estimate or analytical comparison, and the plotted lag axis is in points rather than a physical time unit. The source does not further explain the normalisation convention behind `xcorr(...,'normalized')` or the factor `1/3`.
