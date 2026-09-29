# examples/nmr_diffusion/flow_test_2.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/flow_test_2.m)

This example propagates a 2D phantom under diffusion and flow, without spin interactions. It loads `R1` from `phantom_a.mat` as the initial state, and assigns a uniform two-component flow field. The source specifies periodic boundary conditions and estimates calculation time as minutes.

The domain dimensions are `[0.02 0.02]` on a `[108 90]` grid, with derivative setting `{'period',7}`. Both flow components are `0.2`; the diffusion tensor is diagonal, with `dxx=dyy=5e-5` and cross-components `dxy=dyx=0`. The source does not state units for the domain, flow, or diffusion values. It configures a ghost spin with empty Zeeman and coupling matrices, then creates and inflates the transport generator using `v2fplanck(spin_system,parameters)` and `inflate`.

The trajectory call is `evolution(spin_system,F,[],R1(:),5e-4,200,'trajectory')`. Each trajectory column is reshaped to `108-by-90` and shown with `imagesc`; each frame is followed by `drawnow` and a `0.025`-second pause.