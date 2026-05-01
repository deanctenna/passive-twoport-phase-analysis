# Worst-Case Phase Analysis of a Passive Two-Port Network

## Introduction

This repository contains the MATLAB code developed for the dissertation project
“Worst-Case Phase Analysis of a Passive Two-Port Network”.

The project studies phase-related behaviour in a passive two-port electrical
network. In the one-port case, phase angle and power factor are obtained
directly from a single complex response. In the two-port case, the frequency
response becomes a complex matrix, so a single scalar phase angle is no longer
directly available.

The code in this repository supports the numerical part of the project.

## Files Included

### `plot_right_half_normalized_numerical_range.m`
This script approximates and plots the right-half part of the normalised
numerical range in the complex plane.

### `solve_min_normalized_real.m`
This function computes the worst-case scalar quantity numerically and returns
the corresponding direction and iteration history.

## Software Requirements

- MATLAB
- A recent version such as MATLAB R2024b is recommended

No special third-party libraries are required.

## How to Run the Code

Open MATLAB and set the working directory to the folder containing these files.

To plot the right-half normalised numerical range, run:

```matlab
plot_right_half_normalized_numerical_range
