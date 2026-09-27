# examples/nmr_overtone/cpmas_valine_accum.m

- Signature: `cpmas_valine_accum()`

## Purpose

Shows the MAS cross-polarisation accumulation profile between protons and the 14N overtone transition in N-acetylvaline, using the Fokker–Planck formalism. The source estimates hours of computation and credits Ilya Kuprov, M. Carravetta, and M. Concistre.

## Physical / mathematical content

The source cites the valine quadrupolar-tensor data ([DOI](http://dx.doi.org/10.1039/c4cp03994g)). It specifies 14N and 1H at 14.10220742 T; the 14N quadrupolar tensor is set by 3.21 MHz and η=0.27, with the listed Zeeman eigenvalues and Euler angles. Damping relaxation is used with diagonal retention, zero equilibrium, and rate 2000; the sphten-liouv basis is unapproximated.

## Numerical / algorithmic content

The calculation uses the 800-point spherical powder grid `rep_2ang_800pts_sph`, rank 9, a MAS rate of −19.840 kHz, and a 70–105 kHz spectrum window sampled at 256 points. It steps the RF contact duration from 10 to 100 μs in ten increments and computes each spectrum using `singlerot` with `@overtone_cp`.

## Implementation structure

The function builds the spin system and basis, sets the MAS and RF preparation/detection operators, loops over contact durations, and plots each spectrum in a panel.
