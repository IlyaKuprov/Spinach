# etc/textbook/rlx_nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_nqi.m`
- Signature: `[r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)`
- Total lines: 88

## Purpose

Redfield theory expressions for quadrupolar relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- I -nuclear spin quantum number
- omega -nuclear Zeeman frequency, rad/s
- C_q -quadrupolar coupling constant,
- e^2*q*Q/h in Hz
- eta_q -quadrupolar tensor asymmetry
- tau_c -rotational correlation time, seconds

## Outputs

- r1 -longitudinal relaxation rate, Hz
- r2 -transverse relaxation rate, Hz
- t1 -longitudinal relaxation time, seconds
- t2 -transverse relaxation time, seconds

## Implementation structure

- Redfield theory expressions for quadrupolar relaxation rates,
- isotropic tumbling in liquid phase. Syntax:
- [r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)
- I -nuclear spin quantum number
- omega -nuclear Zeeman frequency, rad/s
- C_q -quadrupolar coupling constant,
- e^2*q*Q/h in Hz
- eta_q -quadrupolar tensor asymmetry
- tau_c -rotational correlation time, seconds
- r1 -longitudinal relaxation rate, Hz
- r2 -transverse relaxation rate, Hz
- t1 -longitudinal relaxation time, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `eeqq2nqi()`, `blinv()`, `spden()`, `isscalar()`.
