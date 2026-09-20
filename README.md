# Pod Silencer Model Reproduction

## Overview

This folder contains a MATLAB implementation used to reproduce the
analytical results of the selected pod-silencer study.

## Method

The implementation uses the analytical formulation reported in the
reference study to determine the acoustic response of the silencer.

The characteristic equation is solved numerically to obtain the complex
axial wavenumber. The resulting wavenumber is then used in transfer-matrix
calculations to determine transmission loss.

## Numerical Method

- Symbolic formulation of the characteristic matrix
- Bessel and Neumann functions for radial acoustic fields
- Newton-Raphson iteration for complex axial wavenumber
- Frequency-by-frequency root continuation
- Transfer-matrix formulation
- Transmission-loss calculation

## Parameters Investigated

The implementation evaluates the effect of:

- Pod radius
- Absorptive-section geometry
- Frequency

## Results

The script generates transmission-loss curves corresponding to the
configurations investigated in the reference study.

## Reference

[Add the full citation of the paper here.]
