# Passive Two-Port Phase Analysis

This repository contains the MATLAB code developed for the dissertation project
“Worst-Case Phase Analysis of a Passive Two-Port Network”.

The code is used to:
- plot the right-half normalised numerical range,
- compute the worst-case scalar quantity,
- compute the corresponding worst-case phase response.

## Files
- `plot_right_half_normalized_numerical_range.m`
- `solve_min_normalized_real.m`

## How to use
Run the MATLAB scripts from the MATLAB working directory. The optimisation routine
computes the worst-case quantity, while the plotting routine visualises the
normalised numerical range.

## Notes
The current implementation is designed for the two-port case and may be extended
in future work to larger multiport networks.
