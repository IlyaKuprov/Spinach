# kernel/utilities/isnucleus.m

`verdict=isnucleus(spin_spec)` returns a logical scalar distinguishing a valid nuclear specification from particles and abstract modes. The character input is first checked by the numerical `spin` API; invalid labels or unavailable simulation data therefore retain their lookup diagnostic.

Canonical nuclear labels begin with a mass number, including explicit isomer labels such as `99Tc_m`. Electron/positron, neutron/antineutron, muon, antiproton, and hyperon labels, and the ghost/cavity/phonon/transmon specifications do not. `1H` remains the proton's nuclear row. This naming convention permits the classifier to recognise newly supported physical particles without accessing metadata or maintaining a second particle inventory.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isnucleus.m>

Wiki: <https://spindynamics.org/wiki/index.php?title=isnucleus.m>
