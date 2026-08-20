clear;
clc;
close all;

% ===============================================================
% FINAL REAL-MARKET EXPERIMENT
%
% Compare:
%
%   two marginals:   T1 and T3
%
% against
%
%   three marginals: T1, T2 and T3
%
% using the refined SPY market-implied marginals recovered in Stage 2.
%
% IMPORTANT:
% The MOT variable is
%
%       M_t = S0 * S_t/F_t,
%
% so the original payoff (S3-K*S1)^+ must be transformed before the
% existing MOT solvers are called.
% ===============================================================

cvx_clear;

cfg3 = market_mot_config();

if ~exist(cfg3.outputDir,'dir')
    mkdir(cfg3.outputDir);
end

%% ---------------------------------------------------------------
% Load refined Stage-2 marginals
% ---------------------------------------------------------------

marketMarginals = load_market_marginals();

S0 = marketMarginals.spot;

fineSupport = marketMarginals.supportNormalized(:);
fineProb = marketMarginals.probabilities;

F = marketMarginals.forwardEstimates(:);
D = marketMarginals.discountFactors(:);

F1 = F(1);
F3 = F(3);
D3 = D(3);

%% ---------------------------------------------------------------
% Compress marginals onto the main MOT support
% ---------------------------------------------------------------

compression = compress_market_marginals_cvx( ...
    fineSupport, ...
    fineProb, ...
    cfg3.mainSupportNormalized);

supportNorm = compression.support;
P = compression.probabilities;

% IMPORTANT NUMERICAL SCALING:
% Solve the MOT problem directly on the normalized state X_t = S_t/F_t.
% These marginals all have mean 1 and the support is O(1), which is much
% better conditioned for CVX/SDPT3 than using the monetary scale S0*X_t.
supportMOT = supportNorm;

p1 = P(1,:)';
p2 = P(2,:)';
p3 = P(3,:)';

fprintf('Main MOT support points: %d\n',numel(supportNorm));

fprintf('Compression CVX status : %s\n', ...
        compression.status);

fprintf('Compression call RMSE  : %.3e\n', ...
        compression.rmseCall);

fprintf('Compression max error  : %.3e\n', ...
        compression.maxAbsCallError);

fprintf('Compressed means       : %.12f  %.12f  %.12f\n\n', ...
        compression.means(1), ...
        compression.means(2), ...
        compression.means(3));

%% ---------------------------------------------------------------
% Verify compressed marginals using existing helper
% ---------------------------------------------------------------

cx12 = check_convex_order( ...
    supportMOT,p1, ...
    supportMOT,p2);

cx23 = check_convex_order( ...
    supportMOT,p2, ...
    supportMOT,p3);

if ~cx12.inConvexOrder
    error('Compressed mu1 and mu2 fail convex order.');
end

if ~cx23.inConvexOrder
    error('Compressed mu2 and mu3 fail convex order.');
end

%% ---------------------------------------------------------------
% Save compressed marginals
% ---------------------------------------------------------------

compressedMarginals = struct();

compressedMarginals.supportNormalized = supportNorm;
compressedMarginals.supportMartingale = S0*supportNorm;
compressedMarginals.supportMOTNormalized = supportMOT;
compressedMarginals.probabilities = P;
compressedMarginals.spot = S0;
compressedMarginals.forwardEstimates = F;
compressedMarginals.discountFactors = D;

save(fullfile(cfg3.outputDir, ...
    'compressed_market_marginals.mat'), ...
    'compressedMarginals');

marginalTable = table( ...
    supportNorm, ...
    S0*supportNorm, ...
    p1,p2,p3, ...
    'VariableNames',{ ...
    'NormalizedSupport','MartingaleSupport', ...
    'Probability_T1','Probability_T2','Probability_T3'});

writetable(marginalTable, ...
    fullfile(cfg3.outputDir, ...
    'compressed_market_marginals.csv'));

%% ---------------------------------------------------------------
% Payoff transformation
% ---------------------------------------------------------------
%
% M_t = S0*S_t/F_t
%
% therefore
%
% S_t = (F_t/S0) M_t.
%
% Thus
%
% (S3-K*S1)^+
%
% = (F3/S0) ...
%   (M3 - K*(F1/F3)*M1)^+.
%
% The time-zero option value is then discounted by D3.
%
% Hence:
%
%   K_MOT = K * F1/F3
%
%   Price_0 = D3*(F3/S0)*MOT_value.

Kconversion = F1/F3;
priceScale = D3*F3;

fprintf('Payoff transformation\n');
fprintf('F1                 = %.8f\n',F1);
fprintf('F3                 = %.8f\n',F3);
fprintf('F1/F3              = %.8f\n',Kconversion);
fprintf('D3                 = %.8f\n',D3);
fprintf('Dollar price scale = %.8f\n\n',priceScale);

%% ---------------------------------------------------------------
% Storage
% ---------------------------------------------------------------

K = cfg3.Kgrid;
qN = numel(K);

Kmot = zeros(qN,1);

Lower2 = zeros(qN,1);
Upper2 = zeros(qN,1);
Width2 = zeros(qN,1);

Lower3 = zeros(qN,1);
Upper3 = zeros(qN,1);
Width3 = zeros(qN,1);

Width2Spot = zeros(qN,1);
Width3Spot = zeros(qN,1);

LowerTightening = zeros(qN,1);
UpperTightening = zeros(qN,1);
WidthReduction = zeros(qN,1);
RelativeWidthReduction = nan(qN,1);

MaxFeasResidual2 = zeros(qN,1);
MaxFeasResidual3 = zeros(qN,1);

StatusLower2 = strings(qN,1);
StatusUpper2 = strings(qN,1);
StatusLower3 = strings(qN,1);
StatusUpper3 = strings(qN,1);

%% ---------------------------------------------------------------
% Solve the market comparison
% ---------------------------------------------------------------

for q = 1:qN

    Kmot(q) = K(q)*Kconversion;

    %% Two-marginal problem

    r2 = mot_bounds_cvx( ...
        supportMOT,p1, ...
        supportMOT,p3, ...
        Kmot(q));

    %% Three-marginal problem

    r3 = mot_bounds_3marginal_cvx( ...
        supportMOT,p1, ...
        supportMOT,p2, ...
        supportMOT,p3, ...
        Kmot(q));

    %% Convert MOT expectations to time-zero SPY option prices

    Lower2(q) = priceScale*r2.lower;
    Upper2(q) = priceScale*r2.upper;
    Width2(q) = Upper2(q)-Lower2(q);

    Lower3(q) = priceScale*r3.lower;
    Upper3(q) = priceScale*r3.upper;
    Width3(q) = Upper3(q)-Lower3(q);

    Width2Spot(q) = Width2(q)/S0;
    Width3Spot(q) = Width3(q)/S0;

    LowerTightening(q) = Lower3(q)-Lower2(q);
    UpperTightening(q) = Upper2(q)-Upper3(q);

    WidthReduction(q) = Width2(q)-Width3(q);

    if Width2(q) > cfg3.displayTolerance
        RelativeWidthReduction(q) = ...
            WidthReduction(q)/Width2(q);
    end

    %% Nesting checks

    if Lower3(q) < Lower2(q)-cfg3.nestingTolerance

        warning(['Computed three-marginal lower bound is below the computed ' ...
                 'two-marginal lower bound at K=%.2f. This cannot be ' ...
                 'explained from feasibility residuals alone; run ' ...
                 'run_market_optimality_validation.m.'],K(q));
    end

    if Upper3(q) > Upper2(q)+cfg3.nestingTolerance

        warning(['Computed three-marginal upper bound is above the computed ' ...
                 'two-marginal upper bound at K=%.2f. This cannot be ' ...
                 'explained from feasibility residuals alone; run ' ...
                 'run_market_optimality_validation.m.'],K(q));
    end

    %% Feasibility diagnostics

    MaxFeasResidual2(q) = max([ ...
        r2.diagnostics.lower.rowResidual, ...
        r2.diagnostics.lower.colResidual, ...
        r2.diagnostics.lower.martingaleResidual, ...
        r2.diagnostics.lower.massResidual, ...
        r2.diagnostics.upper.rowResidual, ...
        r2.diagnostics.upper.colResidual, ...
        r2.diagnostics.upper.martingaleResidual, ...
        r2.diagnostics.upper.massResidual]);

    MaxFeasResidual3(q) = max([ ...
        r3.diagnostics.lower.mu1Residual, ...
        r3.diagnostics.lower.mu2Residual, ...
        r3.diagnostics.lower.mu3Residual, ...
        r3.diagnostics.lower.martingale12Residual, ...
        r3.diagnostics.lower.martingale23Residual, ...
        r3.diagnostics.lower.massResidual, ...
        r3.diagnostics.upper.mu1Residual, ...
        r3.diagnostics.upper.mu2Residual, ...
        r3.diagnostics.upper.mu3Residual, ...
        r3.diagnostics.upper.martingale12Residual, ...
        r3.diagnostics.upper.martingale23Residual, ...
        r3.diagnostics.upper.massResidual]);

    StatusLower2(q) = string(r2.statusLower);
    StatusUpper2(q) = string(r2.statusUpper);

    StatusLower3(q) = string(r3.statusLower);
    StatusUpper3(q) = string(r3.statusUpper);

    fprintf('Completed K = %.2f\n',K(q));
end

%% ---------------------------------------------------------------
% Clean tiny display values
% ---------------------------------------------------------------

tol = cfg3.displayTolerance;

Lower2(abs(Lower2) < tol) = 0;
Upper2(abs(Upper2) < tol) = 0;
Width2(abs(Width2) < tol) = 0;

Lower3(abs(Lower3) < tol) = 0;
Upper3(abs(Upper3) < tol) = 0;
Width3(abs(Width3) < tol) = 0;

LowerTightening(abs(LowerTightening) < tol) = 0;
UpperTightening(abs(UpperTightening) < tol) = 0;
WidthReduction(abs(WidthReduction) < tol) = 0;

%% ---------------------------------------------------------------
% Results table
% ---------------------------------------------------------------

results = table( ...
    K,Kmot, ...
    Lower2,Upper2,Width2,Width2Spot, ...
    Lower3,Upper3,Width3,Width3Spot, ...
    LowerTightening,UpperTightening, ...
    WidthReduction,RelativeWidthReduction, ...
    MaxFeasResidual2,MaxFeasResidual3, ...
    StatusLower2,StatusUpper2, ...
    StatusLower3,StatusUpper3);

disp(results);

writetable(results, ...
    fullfile(cfg3.outputDir, ...
    'market_two_vs_three_results.csv'));

%% ---------------------------------------------------------------
% Benchmark at K = 1
% ---------------------------------------------------------------

[~,idx1] = min(abs(K-1));

fprintf('\nFINAL MARKET COMPARISON AT K = 1\n');

fprintf('Original K                 = %.4f\n',K(idx1));
fprintf('MOT strike multiplier      = %.8f\n',Kmot(idx1));

fprintf('\nTwo-marginal interval\n');
fprintf('Lower                      = %.8f\n',Lower2(idx1));
fprintf('Upper                      = %.8f\n',Upper2(idx1));
fprintf('Width                      = %.8f\n',Width2(idx1));

fprintf('\nThree-marginal interval\n');
fprintf('Lower                      = %.8f\n',Lower3(idx1));
fprintf('Upper                      = %.8f\n',Upper3(idx1));
fprintf('Width                      = %.8f\n',Width3(idx1));

fprintf('\nTightening\n');
fprintf('Lower-bound increase       = %.8f\n', ...
        LowerTightening(idx1));

fprintf('Upper-bound decrease       = %.8f\n', ...
        UpperTightening(idx1));

fprintf('Absolute width reduction   = %.8f\n', ...
        WidthReduction(idx1));

fprintf('Raw computed relative reduction = %.2f%%\n', ...
        100*RelativeWidthReduction(idx1));

fprintf('\nMaximum feasibility residuals\n');
fprintf('Two marginal               = %.3e\n', ...
        max(MaxFeasResidual2));

fprintf('Three marginal             = %.3e\n', ...
        max(MaxFeasResidual3));

fprintf(['NOTE: feasibility residuals measure constraint satisfaction only.\n' ...
         'They do not certify objective optimality. Run the explicit\n' ...
         'primal-dual validation script before interpreting tiny nesting\n' ...
         'differences or percentage width reductions.\n']);


fprintf('\nNesting diagnostics (dollar prices)\n');
fprintf('Minimum lower tightening  = %.6e\n',min(LowerTightening));
fprintf('Minimum upper tightening  = %.6e\n',min(UpperTightening));

%% ===============================================================
% FIGURE 1: market two- versus three-marginal bounds
% ===============================================================

figure('Color','w');

plot(K,Lower2,'-o', ...
    'LineWidth',1.2,'MarkerSize',4, ...
    'DisplayName','Two-marginal lower');

hold on;

plot(K,Upper2,'-s', ...
    'LineWidth',1.2,'MarkerSize',4, ...
    'DisplayName','Two-marginal upper');

plot(K,Lower3,'--o', ...
    'LineWidth',1.3,'MarkerSize',4, ...
    'DisplayName','Three-marginal lower');

plot(K,Upper3,'--s', ...
    'LineWidth',1.3,'MarkerSize',4, ...
    'DisplayName','Three-marginal upper');

xlabel('Strike multiplier K');
ylabel('Time-zero price');
grid on;
legend('Location','best');

title('SPY market-implied two- versus three-marginal bounds');

exportgraphics(gcf, ...
    fullfile(cfg3.outputDir, ...
    'market_two_vs_three_bounds.pdf'), ...
    'ContentType','vector');

fprintf('\nFinal market comparison complete.\n');
fprintf(['For the dissertation, interpret these as empirical bounds ' ...
         'derived from approximate SPY market-implied marginals.\n']);
