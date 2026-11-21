# Periodic Cartesian secular dipolar field

`ddf_periodic(npts,spacing,demag_factor)` returns the longitudinal multiplier `1/3-kz^2/|k|^2` in FFT order on the physical reciprocal lattice. The uniform mode receives `1/3-demag_factor` instead. Each spacing is the corresponding period divided by its voxel count. This boundary model applies to repeated patterns, not isolated objects inside a simulation box.
