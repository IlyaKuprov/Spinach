# examples/esr_sol_swept/sle_esr_sat_b1.m

Demonstrates microwave saturation of an X-band nitroxide through the driven stochastic Liouville steady state. The four fields B1=[0.01,0.1,1,5] G and magnetic tensors follow Figure 1 of Gupta, Dzikovski, and Freed, DOI 10.1007/s00723-021-01390-7. Isotropic rotational diffusion is 1e8 per second, and the model uses electron R1=2e5 and R2=1e5 per second with IME thermalisation at 298 K.

Plots field derivatives of minus the imaginary transverse expectation over 3300–3420 G. The offset-to-field conversion centres the isotropic-g resonance at 3360 G. No independent per-panel amplitude normalisation is applied. This is a uniform-orientation SLE example, not a general orienting-potential MOMD fit. The no-pure-dephasing Lindblad model and thermal amplitude are explicit modelling choices.
