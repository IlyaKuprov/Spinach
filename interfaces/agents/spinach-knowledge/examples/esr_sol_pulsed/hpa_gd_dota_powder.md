# examples/esr_sol_pulsed/hpa_gd_dota_powder.m

- Signature: `hpa_gd_dota_powder()`

## Purpose

Simulates a powder-averaged W-band pulsed ESR spectrum of a Gd(III) DOTA complex, using an ideal pulse and a numerical second-order rotating-frame transformation.

## Model and calculation

The model is a single `E8` spin at 9.40 T, with an isotropic g value of 1.9918 and an axial ZFS tensor. The script uses a Zeeman Liouville basis and averages the lab-frame acquisition over a 12,800-point spherical powder grid. The acquisition spans 6 GHz around a 1.5 GHz offset with 4,096 points, zero-filled to 16,384; the apodised signal is Fourier transformed and plotted. The source estimates a run time of minutes.
