function cfg3 = market_mot_config()
% MARKET_MOT_CONFIG
% Configuration for the final real-market two- versus three-marginal MOT
% comparison.
%
% The refined market calibration is intentionally compressed before the
% three-period MOT problem is solved. This keeps the joint state space
% manageable while preserving the fitted option-price information.

    %% Main strike-multiplier grid for the empirical comparison

    cfg3.Kgrid = (0.80:0.05:1.20)';

    %% Main compressed support
    %
    % 17 support points:
    %   0,
    %   0.60, 0.65, ..., 1.30,
    %   1.40
    %
    % The interior range contains all economically relevant fitted mass,
    % while the two external points retain finite tail support.

    cfg3.mainSupportNormalized = [ ...
        0.00;
        (0.60:0.05:1.30)';
        1.40];

    %% Progressive supports used only for the K=1 robustness check
    %
    % Main:   17 points
    % Medium: 21 points
    % Finer:  23 points
    %
    % These are intentionally smaller than the previous 31-point check,
    % which can be numerically too large for SDPT3 in the three-marginal
    % CVX formulation.

    cfg3.checkSupport21 = [ ...
        0.00;
        linspace(0.60,1.30,19)';
        1.40];

    cfg3.checkSupport23 = [ ...
        0.00;
        linspace(0.60,1.30,21)';
        1.40];

    %% Numerical reporting tolerances

    cfg3.displayTolerance = 1e-7;
    cfg3.nestingTolerance = 1e-6;

    % This threshold only flags a computed sign reversal. It is NOT an
    % optimality tolerance. Use run_market_optimality_validation.m to
    % assess whether a sub-cent difference is resolved by primal-dual gaps.

    %% Output folder

    cfg3.outputDir = 'market_mot_output';
end
