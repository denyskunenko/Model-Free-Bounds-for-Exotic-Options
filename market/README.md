# Market-data numerical experiment

MATLAB code for the SPY 13/30/44-day market-data experiment in the dissertation *Model-Free Bounds for Exotic Options*.

## Requirements

- MATLAB
- CVX
- SDPT3

Keep `SPY_options_2026-08-12.csv` in this folder. Run the scripts from this folder.

## Run order

1. `run_market_diagnostics.m` — loads the SPY option data, matches calls and puts, estimates forwards and writes the data diagnostics.
2. `run_market_marginal_recovery.m` — constructs the three fitted market-implied marginals and saves `market_marginals.mat`. It also creates the recovered-marginals figure used in the dissertation.
3. `run_market_mot_comparison.m` — runs the main 17-point two- versus three-marginal comparison across the strike grid and creates the bounds figure used in the dissertation.
4. `run_market_mot_resolution_check.m` — repeats the comparison at `K = 1` on the 17-, 21- and 23-point supports and creates the support-resolution figure used in the dissertation.
5. `run_dual_self_test_synthetic.m` — checks the explicit two- and three-marginal dual implementations against the synthetic `K = 1` benchmark.
6. `run_market_optimality_validation.m` — performs the full-grid primal-dual objective validation for the main 17-point market experiment.
7. `run_market_optimality_resolution_check.m` — performs the primal-dual objective validation at `K = 1` for the 17-, 21- and 23-point supports.

`00_RUN_ORDER_SPY_MAIN.m` contains the same run sequence in MATLAB form.

## Run files

- `run_market_diagnostics.m`: market-data and forward diagnostics.
- `run_market_marginal_recovery.m`: fitted marginal recovery.
- `run_market_mot_comparison.m`: main market MOT comparison.
- `run_market_mot_resolution_check.m`: support-resolution experiment.
- `run_dual_self_test_synthetic.m`: synthetic primal-dual self-test.
- `run_market_optimality_validation.m`: main full-grid optimality validation.
- `run_market_optimality_resolution_check.m`: support-level optimality validation.

The remaining `.m` files are helper functions called by these run scripts.
