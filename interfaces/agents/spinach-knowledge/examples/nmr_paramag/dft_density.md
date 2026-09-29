# examples/nmr_paramag/dft_density.m

- Signature: `dft_density()`
- Source: [dft_density.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/dft_density.m)

## Purpose and model

Compares point-centre and distributed-density PCS for the Eu(III) complex 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. This is a molecular-complex example, not a carbonic-anhydrase example. The source cites the [paper describing the delocalised-model equation](https://doi.org/10.1039/C4CP03106G).

## Inputs and calculations

- Loads electron probability density, grid extent, coordinates, and spacing from `tetra_py_probden.mat`; reads HFC and susceptibility data from `tetra_py_dft_run.log` with `gparse`. The source says HFCs are in Gauss.
- Normalises the probability density by its three-dimensional trapezoidal integral using `dx^3`.
- Converts `props.chi` to a rank-2 tensor representation, then obtains `chi` with `sphten2mat`.
- Computes point-model PCS with `ppcs(xyz,[0 0 0],chi)`; the point-centre coordinate is explicitly `[0 0 0]`, but its coordinate unit is not stated.
- Computes the distributed result with `kpcs(probden,chi,ext,xyz,'fft')`.
- For each DFT atom except the final atom, maps H, C, and N to `1H`, `13C`, and `14N`, then obtains HFC-derived PCS with `hfc2pcs`.

## Outputs and limits

Plots DFT-HFC PCS against distributed PCS and point-model PCS against distributed PCS; those plot axes are labelled in ppm. It transforms the distributed 3-D PCS grid for display and plots it with the molecular coordinates. The source gives no magnetic-field or temperature value and no units for the coordinates or susceptibility tensor. It notes one outlier caused by a nonzero isotropic HFC contact shift; neither the point model nor the Kuprov-equation calculation includes contact shifts. No spectrum is simulated.
