# Design and Analysis of Novel Sound Attenuators for HVAC Duct Systems

## Overview

Room-to-room crosstalk through HVAC duct systems can provide an unintended path for speech and other sounds to propagate between spaces. This project focuses on the design and analysis of sound attenuators for reducing acoustic crosstalk through HVAC duct systems.

The project uses established analytical and transfer-matrix methods to evaluate the acoustic performance of different crosstalk silencer geometries. Published analytical models are first reproduced and validated against the results reported in the literature. The validated methods are then used to investigate and compare different geometrical configurations.

---

## Objectives

- Study the mechanism of acoustic crosstalk transmission through HVAC duct systems.
- Review existing analytical models for duct acoustics and crosstalk silencers.
- Reproduce selected published analytical models using MATLAB.
- Validate the implemented models by comparing reproduced results with published results.
- Apply established analytical and transfer-matrix methods to evaluate different crosstalk silencer geometries.
- Investigate the influence of geometry and frequency on transmission loss.
- Study the effect of acoustic modes on the performance of the silencers.
- Evaluate pressure-drop characteristics alongside acoustic performance for selected geometries.

---

## Methodology

The project is divided into two main stages:

### Stage 1: Literature Reproduction and Validation

Analytical models reported in published literature are implemented in MATLAB and their results are reproduced.

The reproduction process involves:

1. Studying the theoretical formulation presented in the reference paper.
2. Deriving and implementing the required acoustic equations.
3. Implementing the characteristic equations and boundary conditions.
4. Numerically solving for the required acoustic wavenumbers and modal quantities.
5. Constructing the relevant transfer matrices.
6. Calculating transmission loss over the required frequency range.
7. Comparing the reproduced results with the corresponding published results.

### Stage 2: Evaluation of Crosstalk Silencer Geometries

After establishing and validating the analytical methods, the same approaches are used to evaluate different HVAC crosstalk silencer geometries.

The geometries are compared based on their predicted acoustic performance, with particular focus on transmission loss and frequency-dependent attenuation.

Pressure-drop characteristics will subsequently be considered alongside acoustic performance when evaluating selected configurations.

---

# Literature Reproduction

The analytical models used in this project are based on previously published research. The purpose of the reproduction work is to understand, implement, and validate the analytical formulations before applying them to the evaluation of new geometrical configurations.

The original research papers are **not redistributed in this repository**. The corresponding references are provided in the README files within each paper folder.

---

## Reproduced Paper 1: Pod Silencers

### Overview

The first reproduction implements the analytical model presented in the reference study for the analysis and design of pod silencers.

The MATLAB implementation reproduces the acoustic formulation and transmission-loss calculations for the configurations investigated in the paper.

### Analytical Approach

The implementation includes:

- Formulation of the acoustic characteristic equation.
- Radial acoustic field representation using Bessel and Neumann functions.
- Numerical determination of complex axial wavenumbers.
- Newton-Raphson iteration for solving the characteristic equation.
- Frequency-by-frequency tracking of the required roots.
- Transfer-matrix formulation of the acoustic elements.
- Calculation of the overall transmission loss.

### Parameters Investigated

The implementation is used to investigate the effect of:

- Pod geometry
- Pod radius
- Silencer configuration
- Frequency

### Validation

The calculated transmission-loss results are compared with the corresponding results reported in the reference paper.

This reproduction serves as a validation step before applying the analytical framework to further geometrical configurations.


## Reproduced Paper 2: Comparison of Analytical Methods
### Overview

The second reproduction implements and compares analytical approaches for acoustic propagation and transmission loss in lined circular ducts.

The study is used to understand the different formulations and establish an analytical framework suitable for further evaluation of HVAC silencer configurations.

### Analytical Approach

The MATLAB implementation includes:

Forward and backward propagating acoustic modes.
Mean-flow effects.
Complex axial and radial wavenumbers.
Acoustic impedance modelling.
Characteristic equations for the acoustic modes.
Numerical solution using Newton-Raphson iteration.
Transfer-matrix formulations.
Comparison of multiple analytical approaches.
Parameters Investigated

### Parameters 
The implementation investigates the effects of parameters including:

Mach number
Duct geometry
Lining characteristics
Annular cavity dimensions
Geometrical ratios
Frequency
Validation and Comparison

The predicted transmission-loss characteristics obtained using the different analytical approaches are compared with the corresponding results reported in the reference study.

The comparison is used to understand the applicability and behaviour of the different formulations before applying the analytical methods to the present crosstalk-silencer problem.
