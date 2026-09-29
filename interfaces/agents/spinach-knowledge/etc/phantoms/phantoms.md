# etc/phantoms/phantoms.m

- MATLAB implementation: [etc/phantoms/phantoms.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/phantoms/phantoms.m)

## Purpose

Loads one of the built-in MRI phantom data sets and returns voxelwise relaxation-rate arrays, proton density, physical dimensions and voxel counts.

## Use

```matlab
[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)
```

`ph_name` must be a character array and must name one of the case-sensitive entries below. There is no default; an unrecognised name raises an error.

| `ph_name` | Input and processing | Returned sampling / dimensions |
|---|---|---|
| `boob` | Loads `mri_lab_boob.mat` (variables `boob`, `nvoxels`). Maps its tissue labels to the values in the breast table below. | Assumes 1 mm isotropic voxels; `dims=npts.*[1e-3 1e-3 1e-3]` metres. |
| `brain-highres` | Loads `mri_lab_brain.mat` (`VObj`); `R1Ph=1./VObj.T1`, `R2Ph=1./VObj.T2`, `PDPh=VObj.Rho`. Infinite R1/R2 entries are replaced by zero. | Retains stored sampling. `dims=npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`. |
| `brain-medres` | Same MRiLab fields and infinity handling; takes every second element along each array dimension. | `dims=2*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`. |
| `brain-lowres` | Same MRiLab fields and infinity handling; takes every fourth element along each array dimension. | `dims=4*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`. |

For all cases, `npts=size(PDPh)` and `dims` is a three-element row vector in metres. The brain resolution fields are used in Y, X, Z order by the source. No additional reorientation is performed by this function.

## Breast phantom tissue values

The following source-assigned values are for the `boob` phantom; the source labels these relaxation values as being for 3 T. Units for R1/R2 are not stated in the source, so retain them as supplied rather than assuming a unit. PD is the assigned proton-density value.

| Source label | Tissue | PD | R1 | R2 |
|---:|---|---:|---:|---:|
| -1.0 | Air | 0.00 | 0.00 | 0.00 |
| -2.0 | Skin | 0.40 | 0.82 | 6.49 |
| -4.0 | Muscle | 0.50 | 1.11 | 34.48 |
| +1.1 | Fibroconnective, high water | 0.75 | 0.69 | 18.39 |
| +1.2 | Fibroconnective, medium water | 0.70 | 1.00 | 25.00 |
| +1.3 | Fibroconnective, low water | 0.65 | 1.50 | 35.00 |
| +2.0 | Transitional | 0.60 | 0.90 | 8.00 |
| +3.1 | Fatty, low fat | 0.55 | 2.73 | 14.70 |
| +3.2 | Fatty, medium fat | 0.50 | 2.67 | 16.54 |
| +3.3 | Fatty, high fat | 0.45 | 2.61 | 18.39 |

`R1Ph`, `R2Ph`, and `PDPh` are returned as cubes in MATLAB array layout. For the brain data, the input MAT file's T1/T2 units are inherited by the reciprocals; the source does not specify a unit conversion. Infinite rates are zeroed, but the implementation does not otherwise sanitise the arrays.

## Source link

[Spinach Wiki: phantoms.m](https://spindynamics.org/wiki/index.php?title=phantoms.m)
