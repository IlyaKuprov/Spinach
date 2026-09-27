# examples/nmr_paramag/dft_density.m

- Signature: `dft_density()`

## Purpose

Simulation of the PCS field for the Eu(III) complex of 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. The example imports the electron probability density and DFT hyperfine and susceptibility data, then compares the distributed-density solution with point-model and HFC-derived PCS. The distributed model is based on the [Kuprov-equation paper](http://dx.doi.org/10.1039/C4CP03106G). The source notes that one outlier arises from an isotropic hyperfine contact shift for a nucleus; the point model and Kuprov equation used here do not include contact shifts.

## Physical / mathematical content

The calculation normalizes the imported three-dimensional electron probability density, obtains the susceptibility tensor from the DFT data, and solves the PCS field with `kpcs`. It also computes point-model PCS with `ppcs` and HFC-derived PCS with `hfc2pcs` for comparison.

## Numerical / algorithmic content

The DFT-derived PCS values are generated for the parsed nuclei using the corresponding isotope labels (¹H, ¹³C, or ¹⁴N). The example plots HFC-derived and point-model PCS against the distributed solution and displays the distributed PCS field.

## Implementation structure

- Load and normalize the electron probability density and read the DFT data.
- Derive the susceptibility tensor and compute point-model and distributed PCS.
- Calculate PCS from DFT hyperfine tensors and compare all models; plot the distributed field.
