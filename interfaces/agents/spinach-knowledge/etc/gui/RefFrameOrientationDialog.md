# etc/gui/RefFrameOrientationDialog.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/RefFrameOrientationDialog.cpp)

This dialog edits the rotation of a reference frame, not the eigenvalues of one interaction. Changing the frame affects lab-frame orientations of its descendants and attached tensors. For a meaningful check use an anisotropic tensor attached to a child of a rotated parent: changing the parent should change the observed orientation while the child's local tensor values remain the comparison point. A test with just an isotropic tensor cannot distinguish a correct rotation from no rotation.
