# kernel/includes/redfield_integral_async.m

- Signature: `job_id=brw_compute_kernel(spin_system,w,job_id,upper_lim)`

## Purpose

This is the asynchronous parallel include for Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluation, called from the `relaxation.m` theory blocks. It follows the notation of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2010.12.004); its numerical quadrature method is superseded here by the faster auxiliary-matrix method described in [the second cited paper](http://dx.doi.org/10.1063/1.4928978).

## Theory parameters

- `rlx_onshell`: true selects the back-rotated kernel, which reduces to Redfield theory at zero shift; false selects the Nakajima-Zwanzig resolvent kernel.
- `rlx_shift`: the Laplace evaluation point, in Hz. Redfield theory is the on-shell form at zero shift.

## Algorithm

The include queues asynchronous integral jobs for significant spherical-tensor and correlation-function terms, then gathers their sparse contributions into the relaxation superoperator. The worker evaluates each contribution with the auxiliary-matrix integral routine `expmint`.

## Source documentation

https://spindynamics.org/wiki/index.php?title=redfield_integral_async.m
