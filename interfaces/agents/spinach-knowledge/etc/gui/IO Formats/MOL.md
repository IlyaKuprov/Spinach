# etc/gui/IO Formats/MOL.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/MOL.cpp)

The MOL reader uses the V2000 atom block to recover elements and positions. Fields may abut for long signed coordinates, so a parser change should include a fixed-width case rather than testing only whitespace-separated values. The explicit bond block is not the authority for the GUI's display bonds. Use SpinXML, not MOL, when the objective is to retain spin interactions and reference frames.
