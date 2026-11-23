# Distant dipolar field of a one-axis grating

[ddf_field.m](../../../../../kernel/utilities/ddf_field.m) takes equally spaced Cartesian magnetisation samples in A/m and returns the rotating-frame field at the same positions in tesla. All nonuniform harmonics of a one-axis grating share one angular factor. The spatial mean instead sees the axial demagnetising factor of the macroscopic sample. A spherical sample cancels the field of a uniform magnetisation, not the field of every grating. This bulk, periodic calculation is not a finite-geometry dipole sum.
