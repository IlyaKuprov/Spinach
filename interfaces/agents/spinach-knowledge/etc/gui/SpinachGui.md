# etc/gui/SpinachGui.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/SpinachGui.cpp)

The Windows executable starts a single-threaded-apartment WinForms message loop. With no argument it creates an empty editor; with an argument it passes the first filename to MainForm. This is a document-opening convenience, not a command-line simulation or batch conversion interface. For a quantum chemistry output, use Import and explicitly choose the originating program.
