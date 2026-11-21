# etc/gui/IO Formats/EasySpin.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/EasySpin.cpp)

This writer produces an EPR simulation script with Sys, Exp and Opt and a selected garlic, chili or pepper call. It emits Sys.S=1/2: do not assume an arbitrary multielectron model is faithfully represented. Hyperfine and quadrupolar parameters must remain associated with the same nucleus in nucspinadd. Field-valued hyperfine data and frequency-valued data take different conversion paths; changing an output label is not a unit conversion. Inspect the generated experiment, strain, broadening and dynamics settings before treating the output as a measurement prediction.
