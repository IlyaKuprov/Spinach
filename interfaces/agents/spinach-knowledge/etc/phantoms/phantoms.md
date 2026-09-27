# etc/phantoms/phantoms.m

- Signature: `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)`

## Purpose

Loads a built-in MRI phantom and returns voxelwise longitudinal and transverse relaxation rates, proton density, physical dimensions, and voxel counts. Supported names are `boob`, `brain-highres`, `brain-medres`, and `brain-lowres`; other names raise an error.

## Data and processing

For `boob`, the function loads the bundled breast phantom and maps its tissue labels to proton-density, R1, and R2 values specified for 3 T. It assumes isotropic 1 mm voxels. For the three brain entries, it loads the bundled MRiLab brain phantom, sets R1=1/T1, R2=1/T2, and PD from the MRiLab proton-density field, and replaces infinite relaxation rates with zero. The medium- and low-resolution variants subsample the data by factors of 2 and 4, respectively; the high-resolution variant uses the stored sampling. Physical dimensions account for the corresponding sampling interval.

The result arrays follow MATLAB's voxel-array layout. The function requires `ph_name` to be a character array.

## Inputs

- `ph_name` — name of a supported phantom (see above).

## Outputs

- `R1Ph` — cube of longitudinal relaxation rates.
- `R2Ph` — cube of transverse relaxation rates.
- `PDPh` — cube of proton-density values.
- `dims` — three-element row vector of physical dimensions in metres.
- `npts` — three-element row vector of voxel counts.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=phantoms.m).
