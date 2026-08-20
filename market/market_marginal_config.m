function cfg2 = market_marginal_config()
% MARKET_MARGINAL_CONFIG
% Configuration for Stage 2: recovery of discrete market-implied marginals.
%
% The marginal is fitted in normalized units
%
%       X_T = S_T / F_T,
%
% so that each fitted marginal has mean 1.
%
% For later MOT work we also save the monetary martingale scale
%
%       M_T = S0 * X_T,
%
% whose common mean is S0.

    %% Common finite support in normalized units
    %
    % Observed retained strikes lie roughly between 0.70 and 1.20 times
    % the estimated forward. The extra points provide finite tail support.
    %
    % The same support is used for all three marginals so the convex-order
    % comparison is especially transparent.

    cfg2.supportNormalized = [ ...
        0.00;
        (0.60:0.01:1.30)';
        1.40];

    %% Weighted least-squares floor
    %
    % Observed bid/ask half-spreads can be extremely small in the tails.
    % Without a floor, those tiny quotes would receive excessive weight.
    %
    % This value is in normalized call-price units.

    cfg2.fitScaleFloor = 5e-5;

    %% Numerical tolerances used only for reporting

    cfg2.probabilityTolerance = 1e-8;
    cfg2.convexOrderTolerance = 1e-7;

    %% Output folder

    cfg2.outputDir = 'market_marginal_output';
end
