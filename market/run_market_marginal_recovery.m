clear;
clc;
close all;

% ===============================================================
% REAL MARKET DATA - STAGE 2
%
% Recover three finite, arbitrage-consistent marginals from the SPY
% option data diagnosed in Stage 1.
%
% This script still does NOT calculate the forward-start MOT bounds.
% It prepares and validates mu1, mu2 and mu3 for that final stage.
% ===============================================================

cvx_clear;
% Solver selection is deliberately not changed inside this script.
% Set SDPT3 once in the MATLAB Command Window before the first run:
%     cvx_solver sdpt3

cfg = market_data_config();
cfg2 = market_marginal_config();

%% Output directory

if ~exist(cfg2.outputDir,'dir')
    mkdir(cfg2.outputDir);
end

%% Load the same raw option data used in Stage 1

T = load_spy_option_data(cfg);

%% Construct approximate normalized call observations

[obs,marketInfo] = ...
    build_normalized_market_observations(T,cfg,cfg2);

fprintf('\nStage 2 market observations\n');
fprintf('Total retained OTM observations: %d\n',height(obs));
fprintf('Normalized strike range        : %.4f to %.4f\n', ...
        min(obs.NormalizedStrike), ...
        max(obs.NormalizedStrike));

%% Joint marginal fit

supportNorm = cfg2.supportNormalized;

fit = fit_market_marginals_cvx( ...
    obs,supportNorm,cfg2);

P = fit.probabilities;

%% Monetary martingale scale
%
% M_T = S0 * S_T/F_T
%
% Every marginal therefore has common mean S0.

supportMoney = cfg.spot*supportNorm;

meanMoney = P*supportMoney;

%% ---------------------------------------------------------------
% Convex-order verification using the existing helper
% ---------------------------------------------------------------

cx12 = check_convex_order( ...
    supportMoney,P(1,:)', ...
    supportMoney,P(2,:)', ...
    cfg2.convexOrderTolerance);

cx23 = check_convex_order( ...
    supportMoney,P(2,:)', ...
    supportMoney,P(3,:)', ...
    cfg2.convexOrderTolerance);

%% ---------------------------------------------------------------
% Console diagnostics
% ---------------------------------------------------------------

fprintf('\nCVX fit status: %s\n',fit.status);
fprintf('CVX objective : %.6f\n',fit.objective);

fprintf('\nProbability diagnostics\n');

for t = 1:3
    fprintf(['Maturity %d: mass = %.12f, normalized mean = %.12f, ' ...
             'money mean = %.8f, min p = %.3e\n'], ...
             t,fit.mass(t),fit.means(t),meanMoney(t), ...
             fit.minProbability(t));
end

fprintf('\nConvex-order diagnostics\n');
fprintf('mu1 <=cx mu2: max violation = %.3e\n', ...
        fit.maxConvexOrderViolation12);

fprintf('mu2 <=cx mu3: max violation = %.3e\n', ...
        fit.maxConvexOrderViolation23);

fprintf('Existing helper check mu1 <=cx mu2: %d\n', ...
        cx12.inConvexOrder);

fprintf('Existing helper check mu2 <=cx mu3: %d\n', ...
        cx23.inConvexOrder);

fprintf('\nMarket fit diagnostics\n');
fprintf('Normalized RMSE          : %.6e\n', ...
        fit.rmseNormalized);

fprintf('Dollar RMSE              : %.6f\n', ...
        fit.rmseDollar);

fprintf('Median absolute $ error  : %.6f\n', ...
        fit.medianAbsDollarError);

fprintf('Maximum absolute $ error : %.6f\n', ...
        max(abs(fit.fittedObservations.DollarError)));

fprintf('Fraction inside bid/ask  : %.2f%%\n', ...
        100*fit.fractionInsideBidAsk);

%% ---------------------------------------------------------------
% Save observation-level diagnostics
% ---------------------------------------------------------------

writetable(fit.fittedObservations, ...
    fullfile(cfg2.outputDir, ...
    'market_option_fit_diagnostics.csv'));

writetable(marketInfo, ...
    fullfile(cfg2.outputDir, ...
    'market_maturity_inputs.csv'));

%% ---------------------------------------------------------------
% Save marginal probabilities
% ---------------------------------------------------------------

marginalTable = table( ...
    supportNorm, ...
    supportMoney, ...
    P(1,:)', ...
    P(2,:)', ...
    P(3,:)', ...
    'VariableNames',{ ...
    'NormalizedSupport', ...
    'MartingaleSupport', ...
    'Probability_T1', ...
    'Probability_T2', ...
    'Probability_T3'});

writetable(marginalTable, ...
    fullfile(cfg2.outputDir, ...
    'market_marginals.csv'));

%% ---------------------------------------------------------------
% Save MAT file for the final MOT stage
% ---------------------------------------------------------------

marketMarginals = struct();

marketMarginals.quoteDate = cfg.quoteDate;
marketMarginals.expiries = cfg.expiries;
marketMarginals.daysToExpiry = cfg.daysToExpiry;

marketMarginals.spot = cfg.spot;

marketMarginals.supportNormalized = supportNorm;
marketMarginals.supportMartingale = supportMoney;

marketMarginals.probabilities = P;

marketMarginals.forwardEstimates = ...
    marketInfo.ForwardEstimate;

marketMarginals.discountFactors = ...
    marketInfo.DiscountFactor;

marketMarginals.treasuryRates = ...
    marketInfo.TreasuryRate;

marketMarginals.fitStatus = fit.status;
marketMarginals.fitRMSEDollar = fit.rmseDollar;
marketMarginals.fitFractionInsideBidAsk = ...
    fit.fractionInsideBidAsk;

save(fullfile(cfg2.outputDir, ...
    'market_marginals.mat'), ...
    'marketMarginals');

%% ===============================================================
% FIGURE 3: recovered marginal probabilities
% ===============================================================

figure('Color','w');

plot(supportMoney,P(1,:),'-o', ...
    'LineWidth',1.1,'MarkerSize',4, ...
    'DisplayName','T_1');

hold on;

plot(supportMoney,P(2,:),'-s', ...
    'LineWidth',1.1,'MarkerSize',4, ...
    'DisplayName','T_2');

plot(supportMoney,P(3,:),'-^', ...
    'LineWidth',1.1,'MarkerSize',4, ...
    'DisplayName','T_3');

xlabel('Martingale-scaled SPY level');
ylabel('Probability mass');
grid on;
legend('Location','best');

title('Recovered discrete market-implied marginals');

exportgraphics(gcf, ...
    fullfile(cfg2.outputDir, ...
    'market_marginal_probabilities.pdf'), ...
    'ContentType','vector');

%% Final message

fprintf('\nStage 2 complete.\n');
fprintf(['If the probability, convex-order and market-fit diagnostics ' ...
         'are satisfactory,\n']);
fprintf(['market_marginals.mat is the input for the final real-data ' ...
         'two- versus three-marginal MOT comparison.\n']);
