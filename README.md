# Airfoil Inverse Design via MGM Method

A MATLAB implementation of the Modified Garabedian-McFadden (MGM) inverse design method for aerodynamic shape optimization.

## Project Overview

This project solves the **inverse design problem** in aerodynamics: determining an optimal airfoil shape that produces a desired surface pressure distribution. The algorithm automatically morphs a baseline geometry (NACA 0012) to match a target pressure profile.

## Methodology

- **Algorithm**: Modified Garabedian-McFadden (MGM) inverse method
- **Flow Solver**: Integrated with **XFOIL** for high-fidelity inviscid flow analysis
- **Process**:
  1. XFOIL analyzes current airfoil and computes pressure distribution (`C_p`)
  2. MGM calculates geometric modifications (`Δz`) from `C_p` error
  3. Airfoil shape is updated iteratively until convergence

## Results

Successfully morphed NACA 0012 airfoil toward target design:
- Final geometry shows clear convergence to target `C_p` distribution
- Demonstrates effective computational aerodynamic optimization

## Usage

1. Ensure MATLAB and XFOIL are installed
2. Run `main_code.m`
3. Review generated output files:
   - `airfoil_final.dat` (optimized geometry)
   - `cp_alfa_aoa2.txt` (final pressure distribution)

## Repository Structure
- main_code.m # Main algorithm
- targetCp.csv # Target pressure distribution
- Naca0012_CorrectedFormat.dat # Baseline airfoil
- xfoil_input.txt # XFOIL configuration
- airfoil_final.dat # Optimized airfoil (output)
