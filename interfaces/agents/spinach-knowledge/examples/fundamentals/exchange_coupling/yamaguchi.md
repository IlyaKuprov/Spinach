# examples/fundamentals/exchange_coupling/yamaguchi.m

- MATLAB implementation: [examples/fundamentals/exchange_coupling/yamaguchi.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/exchange_coupling/yamaguchi.m)

- Signature: `yamaguchi()`

## Question and model

This example estimates exchange coupling using the Yamaguchi broken-symmetry DFT treatment for the bistrityl biradical with an alkynyl linker attributed in the source to Olav Schiemann. The script does not set up or run the DFT calculations; it consumes properties parsed from Gaussian log files.

## Inputs and calculation

From MATLAB's current folder, it passes `biradical_singlet.log` and `biradical_triplet.log` to `gparse`, storing the resulting property structures as `props_sing` and `props_trip`. It passes those structures to `brokensymm` as `J=brokensymm(props_sing,props_trip)`. The script displays `J/1e9` with the label `Exchange coupling:` and the unit `GHz`.

## Output and limitations

No fixed numeric value is present in the script: the displayed result depends on the two input logs and the `brokensymm` implementation. The example contains no explicit input-file checks, DFT method or geometry setup, numerical acceptance test, or error handling. The precise equation, sign convention, and calculation details therefore are not specified by this script and should not be inferred from its output label alone. Source: [examples/fundamentals/exchange_coupling/yamaguchi.m](../../../../../../examples/fundamentals/exchange_coupling/yamaguchi.m).