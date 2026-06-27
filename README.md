# A Sieve-Accelerated Quadrature Method for Exact Privacy Accounting in the 2020 U.S. Decennial Census

This repository accompanies the project/paper:

**“A Sieve-Accelerated Quadrature Method for Exact Privacy Accounting in the 2020 U.S. Decennial Census.”**

It provides code and precomputed results for computing exact privacy accounting (e.g., epsilon-delta curves and trade-off curves) for compositions of mechanisms used in the 2020 Census DHC File, with an emphasis on efficient, high-precision computation.

## Repository structure

- **`code/`**
  Core implementation and scripts for generating privacy accounting outputs and plots.

- **`privacy budget allocation/`**
  Records the privacy budget allocations \rho_i used in the 2020 Census DHC File.

- **`results/`**
  Precomputed curve CSV outputs and summary figures.
  - **`results/epsilon_delta_curve/`**
    Epsilon-delta curves characterizing privacy levels for different compositions/paths.
    - **`results/epsilon_delta_curve/eps_delta_results/`**
      Raw epsilon-delta curve CSV outputs for all computed path pairs.
  - **`results/trade_off_curve/`**
    Trade-off curves organized using the same rationale.
    - **`results/trade_off_curve/trade_off_results/`**
      Raw trade-off curve CSV outputs for all computed path pairs.

## Key result files 

### Epsilon-delta curves

The folder `results/epsilon_delta_curve/` contains PDFs that summarize privacy accounting results:

- **`epsilon_delta_curve_all_path_max.pdf`**  
  Characterizes the overall privacy level in the 2020 Census DHC File (corresponding to the right panel of Figure 2).

- **`epsilon_delta_curve_path_main_to_main.pdf`**  
  Characterizes the privacy level of the composed mechanism [\tilde{M}_0, \tilde{M}_0] (corresponding to the enlarged panel in the lower-right corner of Figure 4).

- **`epsilon_delta_curve_path_sensitivity.pdf`**  
  Characterizes the percentage difference between the privacy level of  
  [\tilde{M}_0, \tilde{M}_0] and the overall privacy level (corresponding to the right panel of Figure 5).

The folder **`results/epsilon_delta_curve/eps_delta_results/`** contains the corresponding precomputed epsilon-delta curve CSV outputs.

### Trade-off curves

We organize trade-off curves in **`results/trade_off_curve/`** following a similar rationale to the epsilon-delta curve outputs. The folder **`results/trade_off_curve/trade_off_results/`** contains the corresponding precomputed trade-off curve CSV outputs.

