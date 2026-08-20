function cfg = market_data_config()
% MARKET_DATA_CONFIG
% Central configuration for the first real-market-data diagnostics stage.
%
% This stage does NOT yet construct risk-neutral marginals. It only:
%   - loads the raw SPY quotes;
%   - matches calls and puts;
%   - constructs bid/ask midpoints and spread diagnostics;
%   - approximates short discount factors;
%   - estimates a near-ATM synthetic forward for each expiry;
%   - identifies a preliminary OTM quote set.
%
% SPY options are American-style. Therefore the parity-based forward is
% treated as a diagnostic approximation rather than an exact identity.

    %% Files

    cfg.optionCsv = 'SPY_options_2026-08-12.csv';

    %% Market identifiers

    cfg.symbol = "SPY";
    cfg.quoteDate = datetime(2026,8,12);

    %% Spot input

    cfg.spot = 772.49;

    %% Selected expiries

    cfg.expiries = [ ...
        datetime(2026,8,25);
        datetime(2026,9,11);
        datetime(2026,9,25)];

    cfg.daysToExpiry = [13;30;44];

    %% Treasury inputs
    %
    % Rates are decimal annual rates, not percentages.
    % 1m = 3.78%, 2m = 3.80%, 3m = 3.87%.

    cfg.treasuryTenorDays = [30;60;90];

    cfg.treasuryRates = [ ...
        0.0378;
        0.0380;
        0.0387];

    %% Forward-estimation filters
    %
    % Use near-the-money strikes only. This reduces dependence on deep-ITM
    % American-option quotes, where early-exercise effects are more relevant.

    cfg.parityMoneynessWindow = 0.05;   % +/-5% around spot
    cfg.maxParityRelativeSpread = 0.25; % 25%
    cfg.minParityCandidates = 3;

    %% Preliminary OTM quote filter
    %
    % This is deliberately not the final marginal-construction filter.
    % It is only used to reveal which quotes look usable.

    cfg.maxOtmRelativeSpread = 0.50;    % 50%
    cfg.requirePositiveBid = true;

    %% Output folder

    cfg.outputDir = 'market_diagnostics_output';
end
