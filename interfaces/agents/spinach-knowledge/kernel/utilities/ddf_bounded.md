# Isolated Cartesian secular dipolar field

`ddf_bounded(npts,spacing)` integrates magnetisation in each rectangular source voxel at target cell centres. Signed solid angles of the two axial faces give the longitudinal demagnetising field; the coincident cell receives the spherical-cavity correction. Noncubic self cells generally have a nonzero response. The `(2*npts-1)` embedding contains every distinct source–target displacement: zero extension represents vacuum, not periodically repeated samples. Do not add a sample-shape zero-mode term.
