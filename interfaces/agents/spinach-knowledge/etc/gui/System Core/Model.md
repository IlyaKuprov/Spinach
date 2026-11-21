# etc/gui/System Core/Model.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/System%20Core/Model.cpp)

The GUI model is not MATLAB's spin_system structure. Atoms, interactions and reference frames are mutable objects with change notifications. Positive interaction keys identify physical terms; negative keys are used for calculated chemical bonds. A negative key is not a negative coupling. Adding or removing atoms can change the bonds and the extrema used by the renderer. Do not treat a grid row number as a persistent external spin ID.

A tensor carries a matrix, its frame, its interaction kind, its unit and the participating atoms. When exporting a tensor, resolve its frame through the parent chain; an anisotropic matrix expressed in a child frame is not already in the lab basis. Frame 1 is the lab frame. Re-rooting must leave a consistent tree for observers as well as preserve physical tensor orientations. Molecular susceptibility is a separate optional model tensor in cgs-ppm at the molecular origin, not a spin pair or another entry in the interaction grid.
