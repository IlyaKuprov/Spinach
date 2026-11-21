# etc/gui/SpinachGui.vcxproj

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/SpinachGui.vcxproj)

Build the Windows Forms application, not a MATLAB MEX target. The project uses managed C++ (C++/CLI), MSVC v143, .NET Framework 4.8, the Windows 10 SDK and the Eigen headers beside the project. A Windows OpenGL context is needed for the molecular viewport. No separate SpinachGUI DLL or slider assembly is required.

Both solution configurations target x64. Intermediate objects and binaries go under the distribution's build directory. Release also copies the executable to the distribution's release directory; Debug does not. Building needs neither MATLAB nor an input molecule. Starting an interactive window and running a generated Spinach simulation are separate checks.
