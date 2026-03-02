# A Sieve-Accelerated Quadrature Method for Exact Privacy Accounting in the 2020 U.S. Decennial Census

This repository accompanies the project/paper:

**“A Sieve-Accelerated Quadrature Method for Exact Privacy Accounting in the 2020 U.S. Decennial Census.”**

It provides code and precomputed results for computing **exact privacy accounting** (e.g., ε–δ curves and trade-off curves) for compositions of mechanisms used in the **2020 Census DHC File**, with an emphasis on efficient, high-precision computation.

## Repository structure

- **`Code/`**  
  Core implementation and scripts for generating privacy accounting outputs and plots.

- **`privacy budget allocation/`**  
  Records the privacy budget allocations **\(\rho_i\)** used in the 2020 Census DHC File.

- **`results/`**  
  Precomputed outputs and figures.
  - **`results/epsilon_delta_curve/`**  
    ε–δ curves characterizing privacy levels for different compositions/paths.
  - **`results/trade_off_curve/`**  
    Trade-off curves organized using the same rationale.

## Key result files (ε–δ curves)

The folder `results/epsilon_delta_curve/` contains PDFs that summarize privacy accounting results:

- **`epsilon_delta_curve_all_path_max.pdf`**  
  Characterizes the **overall privacy level** in the 2020 Census DHC File.

- **`epsilon_delta_curve_path_main_to_main.pdf`**  
  The privacy level of the composed mechanism **\([\widetilde{\mathbf{M}}_0, \widetilde{\mathbf{M}}_0]\)**  
  (corresponding to the enlarged panel in the lower-right corner of **Figure 4**).

- **`epsilon_delta_curve_path_sensitivity.pdf`**  
  Characterizes the **percentage difference** between the privacy level of  
  **\([\widetilde{\mathbf{M}}_0, \widetilde{\mathbf{M}}_0]\)** and the **overall privacy level**.

## Trade-off curves

We organize trade-off curves in **`results/trade_off_curve/`** following a similar rationale to the ε–δ curve outputs.

## Getting started

### Requirements
This repository is intended to be run with **Python 3**.  
(If you add a `requirements.txt` later, list dependencies there.)

### Quick run
Browse the scripts in `Code/` and run the relevant entry points to reproduce figures.
For example, if a script is provided for trade-off curves:
```bash
python Code/trade_off_curve.py
