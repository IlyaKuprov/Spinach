# cce_fixture

[MATLAB source](../../../../../kernel/cce/cce_fixture.m)

Syntax: `[sys,inter,fixture]=cce_fixture(name,revision)`

The supported input is ('nv_four_c13','6f7282b6'): an E3 observed spin and four fixed 13C sites at 0.05 T. The immutable physical data derive from MICCoMpy/PyCCE commit d74d3cb683bf3caa07231654314e554f440833b4 (DOI 10.1002/adts.202100254). These identifiers identify a physical dataset, not a transient implementation revision.

sys contains only the magnetic field and isotopes. inter contains full source hyperfine tensors, zero-field splitting, explicit bath dipoles in Hz, electron g tensors, and nuclear shifts in ppm. The source gyromagnetic ratios are mapped to Spinach's canonical input convention. The caller explicitly selects any observed-Sz projection and execution resources. No transition, cluster order, bath state, or numerical policy is hidden in the reader.

fixture contains name, revision, source (repository, commit, DOI), units, the four-row bath-coordinate array in Angstrom, and the original literal 41-point total echo-time row in seconds spanning 0–0.5 ms. Coordinates are metadata, not inter.coordinates: all pair tensors are already supplied, and adding coordinate dipoles would count the interaction twice. This is not a natural-abundance diamond calculation. Unknown names and revisions are rejected rather than replaced with approximate data. Calls consume no random draws and have no mutable cached state.
