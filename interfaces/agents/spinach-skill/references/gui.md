# Build and use the Windows spin-system editor

## Build the editor

From a VS2022 Developer Command Prompt at the repository root, with Desktop C++, C++/CLI support, the .NET Framework 4.8 targeting pack and Windows SDK installed:

~~~bat
msbuild etc\gui\SpinachGui.sln /t:Rebuild /m /p:Configuration=Debug /p:Platform=x64
msbuild etc\gui\SpinachGui.sln /t:Rebuild /m /p:Configuration=Release /p:Platform=x64
release\SpinachGui.exe
~~~

Check the exit status of each build. The immediate outputs are under build/Debug and build/Release; the release copy is for distribution. The application needs an interactive Windows desktop for visual inspection. Do not use MATLAB's MEX compiler, add Eigen to the MATLAB path, build a DLL, or require a separate slider control. The MATLAB test runner is not a GUI test runner.

## Smoke-test a clean build without a sample collection

Create a text file named water.xyz with these contents:

~~~text
3
water, coordinates in angstroms
O  0.000  0.000  0.000
H  0.957  0.000  0.000
H -0.240  0.927  0.000
~~~

1. Start the executable and use **Import**, selecting XYZ. Expect three atoms and the specified coordinates. Geometry and drawn bonds do not imply measured J couplings. If offered orphan removal, leave it off for this geometry-only input.
2. Edit a label and save as SpinXML. Reopen with **Open**, not the XYZ importer. Verify the label, element and coordinates; a successful file write alone does not check the reader.
3. For tensor work, use a nucleus and an interaction kind appropriate to its isotope. Edit a tensor with three distinct principal values, rotate its frame and inspect the principal axes. Change the display scale, then reopen the orientation editor: the tensor values should not have changed. For a frame change compare lab-frame matrices, not Euler angles alone.
4. Use bitmap export or print to file when a picture is required. A bitmap cannot be reopened as an editable model.

## From quantum chemistry output to a calculation

For electronic-structure data choose Gaussian, ORCA or CASTEP. Select the reader for the **program**, not just the filename suffix. Pick one independently identifiable atom, one nonzero diagonal, one off-diagonal and, where available, a coupling pair from the original output. Check their isotope, frame, unit and imported values. When filtering is offered, distinguish exclusion from the model from hiding a term in Visualization. Record the threshold and the shielding reference if either affects the model.

Use **Spinach export** for a MATLAB spin-system function. Check the isotope ordering, empty electron coordinate slots, coupling units and the lab-frame Zeeman matrices in the generated file. A quadrupole tensor should be traceless, and a Gaussian susceptibility should have a centre at the origin. The editable SpinXML writer does not preserve the separate molecular susceptibility field; keep the calculation output. Use a valid MATLAB function name when saving. Run the generated function with Spinach on the MATLAB path before claiming it executes; creating a system and basis does not by itself calculate an experimental observable.

For EasySpin, SIMPSON and SpinEvolution, use the corresponding exporter and inspect the destination-specific channel list, units and experiment settings. Do not infer simulator acceptance from a successful preview or file write.

## Where to investigate discrepancies

If coordinate magnitudes are wrong, inspect the selected reader and its input units before the renderer. If principal values match but the axes differ, compare the local tensor, frame chain and lab-frame rotation. If the saved scientific values agree but the view differs, isolate the atom, bond or tensor display pass, selection state and visibility controls. Use the source-matched knowledge leaves under spinach-knowledge/etc/gui for the specific component; the MATLAB interfaces/spinxml code is a separate implementation.

The editor accepts up to 4096 atoms; electron spin sites do not count towards this limit.
