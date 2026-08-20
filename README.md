# Model-Free Bounds for Exotic Options

This repository contains the MATLAB code used for the numerical experiments in my dissertation *Model-Free Bounds for Exotic Options*. 
The implementation studies finite-dimensional martingale optimal transport bounds for a forward-start option.

The repository is divided into two parts:

- `synthetic/` contains the synthetic two- and three-marginal experiments, including numerical validation, support-resolution checks, and sensitivity to the intermediate marginal.
- `market/` contains the SPY market-data experiment, including marginal recovery from option quotes, two- versus three-marginal MOT bounds, support-resolution checks, and primal-dual optimality diagnostics.

Each folder contains its own `README.md` with the required run order and a short description of the scripts.

The code was developed in MATLAB using CVX with the SDPT3 solver. The market code was developed with help of LLMs (Copilot and Chat GPT).
