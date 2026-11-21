# etc/gui/IO Formats/COSMOS.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/COSMOS.cpp)

The ATOMS section supplies geometry. Within CALCULATED_PROPERTIES the reader interprets CS-TENSOR and J-COUPLINGS, producing shift tensors in ppm and scalar couplings in Hz. This differs from readers that create chemical shieldings requiring a reference conversion. Recognizing the RDC-TENSOR, ORDER-TENSOR, QC-TENSOR, DD-TENSOR or DISTANCES tag does not mean the corresponding physical interaction is imported.
