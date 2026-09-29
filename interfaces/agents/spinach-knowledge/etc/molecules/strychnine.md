# etc/molecules/strychnine.m

- MATLAB implementation: [etc/molecules/strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/strychnine.m)

## Purpose

`strychnine` assembles a Spinach spin-system/interactions pair for strychnine. It supplies isotope selections, scalar shifts, Cartesian coordinates and available scalar couplings; it is a molecular data set, not a simulation driver.

## Use

```matlab
[sys,inter]=strychnine(spins)
```

- `spins` is a cell array of character strings naming the isotopes to retain; the source gives `{'1H','13C'}` as an example. The implementation checks that this is a cell array whose elements are character arrays, then selects matching isotopes. It does not define defaults or validate each isotope name against a fixed list.
- The isotope list contains 22 `1H`, 21 `13C`, and 2 `15N` sites. Request the isotopes needed for a calculation; the source applies an isotope-membership mask to the isotope and interaction data.
- `sys` and `inter` are Spinach input structures. The interaction data include scalar shifts, coordinates, and the available scalar coupling matrix.

## Data coverage and limits

The proton–proton couplings and isotropic shift values are attributed in the source to Berger and Braun, *200 and more NMR experiments: a practical course*. The code singles out the one-bond C18–H18b coupling as an exception, citing [J. Magn. Reson. (2014), DOI 10.1016/j.jmr.2014.02.003](https://doi.org/10.1016/j.jmr.2014.02.003); its stored value is 131.3. Coordinates are for the major conformer cited as [DOI 10.1039/C0CC04114A](https://doi.org/10.1039/C0CC04114A).

- Carbon–carbon couplings are not supplied. The source explicitly limits 13C use to natural-abundance simulations.
- The 15N sites have shifts and coordinates only; no couplings are provided for them.
- CSA tensors are absent, so this data set does not provide CSA contributions for relaxation calculations.
- The source gives numeric shift, coupling and coordinate arrays without declaring their units. Do not infer units from the values alone. The full parameter assignments remain in `etc/molecules/strychnine.m`.

## Source links

- [Spinach Wiki: strychnine.m](https://spindynamics.org/wiki/index.php?title=strychnine.m)
- [DOI 10.1016/j.jmr.2014.02.003](https://doi.org/10.1016/j.jmr.2014.02.003) — cited one-bond C18–H18b coupling.
- [DOI 10.1039/C0CC04114A](https://doi.org/10.1039/C0CC04114A) — cited major-conformer coordinates.
