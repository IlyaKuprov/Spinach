# etc/gui/MainForm.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/MainForm.cpp)

Open and Save refer to the SpinXML document. Import is a different action: the selected file filter chooses a scientific reader, so a .log suffix alone does not distinguish Gaussian from ORCA. The reader populates a new model. Check the result after the import, filtering and shielding-reference dialogs have finished, not just after choosing the file.

The atoms and interactions grids edit the model. The reference-frame tree controls coordinate systems; changing the lab frame is not a camera movement. The visualization controls set display flags and scale factors, not physical coupling constants. Tensor size sliders use logarithmic scaling; scalar interaction cutoffs use the range of the corresponding interaction kind. Their numeric thresholds are in slider tooltips. Double-clicking a colour label changes the rendered interaction colour without recolouring the header.
