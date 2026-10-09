# examples/kinetics/nonlinear/spinless_sink_network.m

Demonstrates reversible `A+Q -> S`, with a proton-bearing reactant, a spin-free quencher, and a spin-free sink. Empty matching traces out the proton in capture and creates unpolarised spin on reverse arrival. Concentrations and magnetic order propagate together through the production state-dependent generator. The example plots populations and surviving proton polarisation, and reports the error in the conserved atom-equivalent concentration `A+Q+2*S`.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
