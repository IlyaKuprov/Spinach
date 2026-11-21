# etc/gui/System Core/orientation.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/System%20Core/orientation.cpp)

Orientation converts among a direction-cosine matrix, Euler angles, an angle and axis, and a quaternion. The Euler-to-quaternion code constructs rotations about z, y and z from alpha, beta and gamma and multiplies them as q_gamma q_beta q_alpha. The internal trigonometric calculations take radians; the dialog's degree fields must be converted at the UI boundary.

For a rotation round trip compare matrices acting on a vector, not literal angle triples or quaternion components. Euler triples need not be unique, and q and -q describe the same rotation. Use a non-axial, non-diagonal test matrix to distinguish a transpose or reversed composition from a physically equivalent parametrisation.
