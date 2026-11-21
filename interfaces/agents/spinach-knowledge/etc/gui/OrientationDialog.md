# etc/gui/OrientationDialog.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/OrientationDialog.cpp)

The Magnitude and Orientation dialog exposes the Cartesian tensor, principal values and eigenvectors as well as alternative orientation representations. Editing a principal value changes the tensor; changing the viewing plane does not. Compare a reconstructed matrix after edits, since eigenvector signs and the choice of axes in a degenerate eigenspace are not unique. Orientation export can preview and copy selected quantities, including Euler angles and the Wigner matrix, independently of exporting the full spin system.
