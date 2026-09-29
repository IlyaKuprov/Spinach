# kernel/optimcon/penalty.m

Signature: `[pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)`


For a real numeric waveform `wf`, rows are controls and columns are time samples. The scalar `pen_term` and requested derivative outputs are initialised to zero: `pen_grad` has the shape of `wf`, while `pen_hess` is a `numel(wf)`-by-`numel(wf)` matrix in MATLAB column-major waveform order. It calculates only outputs requested by `nargout`. There is no default penalty type.

Let `N=size(wf,2)`. The implemented choices are:

- `none`: leaves the penalty, gradient and Hessian at zero.
- `NS`: sum of squared waveform entries divided by `N`; gradient is `2*wf/N` and Hessian is `2*I/N`.
- `DNS`: forms `D=fdmat(N,5,2,'wall')` and penalises the squared entries of `wf*D'`, divided by `N`. Its gradient is `2*(wf*D')*D/N` and its Hessian is `2*kron(D'*D,speye(size(wf,1)))/N`.
- `SNS`: penalises squared deviations only where entries are strictly above `cb` or strictly below `fb`, divided by `N`. The gradient is twice the active deviation divided by `N`; the Hessian is diagonal with `2/N` for active entries and zero elsewhere. Entries exactly on either bound are inactive.
- `SNSA`: interprets rows in interleaved Cartesian pairs `[Xa Ya Xb Yb ...]`, computes each pair amplitude, and penalises only amplitude strictly above scalar `cb`, divided by `N`. It maps the amplitude gradient back to Cartesian waveform coordinates and assembles the Hessian in paired-coordinate blocks. Inactive amplitudes are set to 1 for the polar-to-Cartesian derivative calculation to avoid the polar singularity. `fb` is checked but does not enter this case's penalty.

Both bounds must be real numeric scalars or arrays the same size as `wf`; all entries must satisfy `cb>=fb`. For `SNSA`, both bounds must instead be scalars and `wf` must have an even number of rows. The DNS case requires at least five time samples. `type` must be a character string; an unrecognised value raises an error.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/penalty.m)
[Spinach Wiki](https://spindynamics.org/wiki/index.php?title=penalty.m)
