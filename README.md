# Worst-Case Phase Analysis of a Passive Two-Port Network

## Introduction

This repository contains the MATLAB code developed for the dissertation project
**“Worst-Case Phase Analysis of a Passive Two-Port Network”**.

The project studies phase-related behaviour in a passive two-port electrical
network. In the classical one-port case, phase angle and power factor are
obtained directly from a single complex response. In the two-port case, the
frequency response becomes a complex matrix, so a single scalar phase angle is
no longer directly available.

To study this problem, the dissertation defines a worst-case phase-related
scalar quantity from a frequency-dependent matrix and computes it numerically.
The code in this repository supports the numerical part of the project.

## Role of the Software in the Project

The repository contains two main MATLAB routines:

- `plot_right_half_normalized_numerical_range.m`  
  This script is used to approximate and plot the right-half part of the
  normalised numerical range in the complex plane.

- `solve_min_normalized_real.m`  
  This function is used to compute the worst-case scalar quantity by solving
  an optimisation problem involving the real part of a normalised matrix
  expression.

Together, these files support the geometric and scalar numerical results
presented in the dissertation.

## Contextual Overview

For a given complex \(2 \times 2\) matrix \(A\), the project considers scalar
quantities of the form

\[
z=\frac{x^*Ax}{\|x\|_2\|Ax\|_2},
\]

where \(x\) is a nonzero complex direction vector. The normalised numerical
range is approximated by sampling many directions \(x\) and plotting the
resulting complex values.

The main scalar quantity computed in the project is the worst-case value

\[
\eta^*=\min_{x\neq 0}\frac{\operatorname{Re}(x^*Ax)}{\|x\|_2\|Ax\|_2},
\]

together with the corresponding angle

\[
\theta_{\max}=\arccos(\eta^*).
\]

These quantities are used in the dissertation to describe the phase-related
behaviour of a passive two-port network across frequency.

## Software Requirements

The code was written for MATLAB and uses only standard MATLAB functionality.

Recommended environment:
- MATLAB R2024b or a similar recent version

No special third-party libraries are required.

## Files Included

### `plot_right_half_normalized_numerical_range.m`
This script:
- asks the user to enter a \(2 \times 2\) complex matrix \(A\),
- generates many random complex direction vectors,
- computes sampled values of the normalised quantity,
- keeps only the right-half part of the sampled set,
- plots the sampled points together with the unit circle and an outer boundary.

### `solve_min_normalized_real.m`
This function:
- accepts a \(2 \times 2\) complex matrix \(A\),
- solves the optimisation problem for the worst-case scalar quantity,
- returns the minimum value, an associated direction vector, and iteration
  history.

The routine uses a Dinkelbach-type outer iteration together with coarse angular
search and local refinement.

## How to Run the Code

### 1. Open MATLAB
Start MATLAB and make sure the current working directory is the folder
containing the repository files.

### 2. Plot the right-half normalised numerical range
Run the plotting script:

```matlab
plot_right_half_normalized_numerical_range
