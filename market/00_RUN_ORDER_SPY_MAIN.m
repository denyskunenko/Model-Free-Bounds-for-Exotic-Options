% RUN ORDER: MAIN SPY 13/30/44 EXPERIMENT -- UPDATED VALIDATION VERSION
% External requirement: CVX with SDPT3.
% The option CSV is included in this folder.
%
% IMPORTANT UPDATE:
% Feasibility residuals do not certify objective optimality. The two new
% validation scripts solve explicit finite duals and report primal-dual
% objective gaps. Run them before finalising Sections 24--25.

%% 1. Raw-data / forward diagnostics
clear all
cvx_clear
cvx_solver sdpt3
run_market_diagnostics

%% 2. Recover refined finite marginals
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_market_marginal_recovery

%% 3. Main 17-point primal comparison and figures
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_market_mot_comparison


%% OPTIONAL BUT RECOMMENDED: quick self-test of the new dual formulations
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_dual_self_test_synthetic

%% 4. REQUIRED: explicit primal-dual validation over the full K grid
% This is the key new run for the nesting/optimality issue.
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_market_optimality_validation

%% 5. Original K=1 support-resolution comparison
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_market_mot_resolution_check

%% 6. REQUIRED: K=1 primal-dual validation on 17/21/23 supports
% clear all
% cvx_clear
% cvx_solver sdpt3
% run_market_optimality_resolution_check
