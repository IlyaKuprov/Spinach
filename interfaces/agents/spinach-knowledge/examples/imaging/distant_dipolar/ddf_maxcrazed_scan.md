# MAXCRAZED acquisition and refinement

Use the common scan to vary the maximum substep and voxel counts independently while keeping the timing, coil, chemical shifts, and frequency bins fixed. The inter-lobe gap and the recovery from the returning lobe to acquisition are independent times, in seconds. Unencoded directions have uniform initial conditions, RF, detection, and reflecting walls. Their solution remains constant along these directions, so only the constant Fourier mode is retained there. Along encoded directions the supplied voxel counts determine resolution.
