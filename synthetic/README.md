# Model-Free Bounds for Exotic Options — MATLAB Code

MATLAB code for the synthetic numerical experiments in Sections 8 and 9 of the dissertation *Model-Free Bounds for Exotic Options*.

## Requirements

- MATLAB
- CVX
- SDPT3

CVX must be installed and available on the MATLAB path. Run the scripts from this folder.

## Run order

1. `run_synthetic_experiment_cvx.m`  
   Baseline two-marginal experiment, numerical residual checks, and primal-dual validation at `K = 1`. Produces the two-marginal bound and normalised-width figures.

2. `run_resolution_study_cvx.m`  
   Support-resolution sensitivity study for the two-marginal problem.

3. `run_three_marginal_comparison_cvx.m`  
   Compares two- and three-marginal bounds for the main synthetic intermediate marginal. Produces the bound-comparison and interval-width figures.

4. `run_middle_marginal_sensitivity_cvx.m`  
   Varies the intermediate marginal through the parameter `lambda` and compares the resulting interval widths.

## Files

- `build_synthetic_marginals.m` — baseline first and final synthetic marginals.
- `check_convex_order.m` — checks equality of means and discrete convex order.
- `mot_bounds_cvx.m` — two-marginal primal MOT lower and upper bounds.
- `mot_dual_cvx.m` — finite two-marginal dual problem used for validation.
- `synthetic_binomial_pair.m` — constructs the marginals used in the support-resolution study.
- `build_three_marginals.m` — constructs the main three-marginal synthetic specification.
- `mot_bounds_3marginal_cvx.m` — three-marginal primal MOT lower and upper bounds.
- `build_middle_marginal_lambda.m` — constructs the intermediate marginal used in the sensitivity experiment.
- `run_synthetic_experiment_cvx.m` — runs the baseline two-marginal experiment.
- `run_resolution_study_cvx.m` — runs the support-resolution study.
- `run_three_marginal_comparison_cvx.m` — runs the main two- versus three-marginal comparison.
- `run_middle_marginal_sensitivity_cvx.m` — runs the intermediate-marginal sensitivity experiment.

The scripts write numerical results to CSV files and save only the figures used in the dissertation as PDF files.
