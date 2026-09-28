# kernel/optimcon/penalty.m

- Signature: `[pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)`

For a real waveform `wf` with channels in rows, returns a penalty, a gradient shaped like `wf`, and a Hessian indexed by time point then channel. Let `N=size(wf,2)`:

- `none`: zero penalty, gradient and Hessian.
- `NS`: sum of `wf.^2` divided by `N`; gradient `2*wf/N`, Hessian `2*I/N`.
- `DNS`: sum of `(wf*D').^2` divided by `N`, where `D=fdmat(N,5,2,'wall')`; requires at least five time points.
- `SNS`: sum of squared deviations outside `[fb,cb]`, divided by `N`.
- `SNSA`: squared excess amplitude above `cb`, divided by `N`; rows must be paired `[Xa Ya Xb Yb ...]` and bounds scalar.

`fb` and `cb` are real scalars or arrays shaped like `wf`, with `cb>=fb`.

Source: https://spindynamics.org/wiki/index.php?title=penalty.m