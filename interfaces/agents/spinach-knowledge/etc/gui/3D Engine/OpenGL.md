# etc/gui/3D Engine/OpenGL.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/3D%20Engine/OpenGL.cpp)

The viewport uses Windows OpenGL and display lists for atoms, bond-like interactions and tensor ellipsoids. Coordinates determine atom positions; isotope data determine radii and colours. Opaque atom spheres use filled polygons with blending and polygon smoothing disabled for their draw pass. Transparent selection outlines and tensor wire overlays are separate draws. Do not remove a tensor outline to fix atom facet seams.

Display scale is not a physical tensor unit. The tensor sliders adjust scale logarithmically; camera zoom, principal-axis visibility and NMR/EPR visibility modes change the view, not the model. Validate rendering with isotropic and anisotropic cases and compare a rotated tensor to its principal axes. A picture alone does not establish that an exported tensor has the correct sign or units.
