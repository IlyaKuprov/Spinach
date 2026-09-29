# experiments/imaging/qxspen_kernel.m

Source: [MATLAB on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/qxspen_kernel.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=qxspen_kernel.m)

- Signature: `[K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)`

## Purpose and parameters

This function constructs the QxSPEN distortion kernel and its derivatives; it is a numerical kernel, not an imaging pulse sequence. `FOVy` is the Y field of view in mm; `NSR` is the reconstruction grid-point count (finite integer at least 2); `Nyacq` is the acquired-point count (finite positive integer); `gam` is linear-phase slope in rad/mm; and `del` is constant phase in radians. The source does not document physical meanings or units for `alp` and `bet`; it describes them as uncertain parameters and directs readers to Ke Dai, so no units are inferred here.

## Kernel and output axes

The acquired coordinate vector spans `-FOVy/2` to `FOVy/2` with `Nyacq` points; reconstruction coordinates span the same endpoints with `NSR` points. For each reconstruction point, the code integrates over its cell using a 100-point `yint` grid and `trapz`, with a sinc term `sinc(alp*(yint-yacq))`, the quadratic phase `exp(1i*bet*(yint-yacq).^2)`, and additional quadratic, linear, and constant phase factors. Thus `K` has `Nyacq` rows (acquired Y locations) and `NSR` columns (reconstruction Y locations); `dK_dgam` and `dK_ddel` have the same shape. The source computes them as `1i*yacq.*K` and `1i*K`, respectively, then scales all three outputs by the same `norm(K,2)` factor from the unscaled kernel.

The code explicitly notes that the fixed integration point count is a dangerous numerical-integration stage and should be replaced by an analytical expression; it also marks adaptive accuracy control as TODO. The fixed count is part of this implementation, not evidence of a validated error bound. No run result, numeric worked example, or DOI is present in the source or existing page.
