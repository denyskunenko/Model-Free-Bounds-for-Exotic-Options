clear;
clc;
close all;

% ===============================================================
% OBJECTIVE-LEVEL VALIDATION: MAIN 17-POINT 13/30/44 EXPERIMENT
%
% This script supplements the feasibility checks in
% run_market_mot_comparison.m with explicit finite primal-dual checks.
%
% IMPORTANT:
% Small marginal/martingale residuals only show approximate feasibility.
% They do NOT by themselves show that the objective value is close to the
% optimum. Here each two- and three-marginal lower/upper primal problem is
% therefore paired with its explicit finite dual.
%
% Run this script before making claims about sub-cent nesting differences
% or the exact percentage tightening in the market experiment.
% ===============================================================

cvx_clear;
cvx_solver sdpt3;

cfg3 = market_mot_config();

if ~exist(cfg3.outputDir,'dir')
    mkdir(cfg3.outputDir);
end

%% Load and compress the refined marginals

marketMarginals = load_market_marginals();

S0 = marketMarginals.spot;
fineSupport = marketMarginals.supportNormalized(:);
fineProb = marketMarginals.probabilities;
F = marketMarginals.forwardEstimates(:);
beta = marketMarginals.discountFactors(:);

F1 = F(1);
F3 = F(3);
beta3 = beta(3);

Kconversion = F1/F3;
priceScale = beta3*F3;

compression = compress_market_marginals_cvx( ...
    fineSupport,fineProb,cfg3.mainSupportNormalized);

support = compression.support(:);
P = compression.probabilities;

p1 = P(1,:)';
p2 = P(2,:)';
p3 = P(3,:)';

fprintf('\nOBJECTIVE-LEVEL VALIDATION\n');
fprintf('Support points           : %d\n',numel(support));
fprintf('Compression status       : %s\n',compression.status);
fprintf('Compression RMSE         : %.3e\n',compression.rmseCall);
fprintf('Dollar price scale       : %.10f\n',priceScale);
fprintf('F1/F3                    : %.10f\n\n',Kconversion);

%% Storage

K = cfg3.Kgrid(:);
qN = numel(K);
Kmot = K*Kconversion;

% Primal and dual objective values in dollars.
L2P = nan(qN,1); L2D = nan(qN,1);
U2P = nan(qN,1); U2D = nan(qN,1);
L3P = nan(qN,1); L3D = nan(qN,1);
U3P = nan(qN,1); U3D = nan(qN,1);

GapL2 = nan(qN,1); GapU2 = nan(qN,1);
GapL3 = nan(qN,1); GapU3 = nan(qN,1);
MaxGap2 = nan(qN,1); MaxGap3 = nan(qN,1);

RawLowerTightening = nan(qN,1);
RawUpperTightening = nan(qN,1);
RawWidthReduction = nan(qN,1);
RawRelativeWidthReductionPct = nan(qN,1);

% Approximate primal-dual brackets for the exact finite-LP quantities.
LowerTighteningBracketLo = nan(qN,1);
LowerTighteningBracketHi = nan(qN,1);
UpperTighteningBracketLo = nan(qN,1);
UpperTighteningBracketHi = nan(qN,1);
WidthReductionBracketLo = nan(qN,1);
WidthReductionBracketHi = nan(qN,1);

MaxPrimalFeasResidual2 = nan(qN,1);
MaxPrimalFeasResidual3 = nan(qN,1);
MaxDualHedgeViolation2 = nan(qN,1);
MaxDualHedgeViolation3 = nan(qN,1);

PrimalSlvTol2 = nan(qN,1);
DualSlvTol2 = nan(qN,1);
PrimalSlvTol3 = nan(qN,1);
DualSlvTol3 = nan(qN,1);

StatusL2P = strings(qN,1); StatusU2P = strings(qN,1);
StatusL2D = strings(qN,1); StatusU2D = strings(qN,1);
StatusL3P = strings(qN,1); StatusU3P = strings(qN,1);
StatusL3D = strings(qN,1); StatusU3D = strings(qN,1);
RunStatus = strings(qN,1);

%% Full strike-grid validation

for q = 1:qN

    fprintf('\n------------------------------------------------\n');
    fprintf('Validating K = %.2f (K_MOT = %.10f)\n',K(q),Kmot(q));
    fprintf('------------------------------------------------\n');

    try
        %% Two marginal: primal and explicit dual
        p2res = mot_bounds_cvx(support,p1,support,p3,Kmot(q));
        d2res = mot_dual_cvx(support,p1,support,p3,Kmot(q));

        %% Three marginal: primal and explicit dual
        p3res = mot_bounds_3marginal_cvx( ...
            support,p1,support,p2,support,p3,Kmot(q));
        d3res = mot_dual_3marginal_cvx( ...
            support,p1,support,p2,support,p3,Kmot(q));

        %% Dollar objective values
        L2P(q) = priceScale*p2res.lower;
        L2D(q) = priceScale*d2res.lower;
        U2P(q) = priceScale*p2res.upper;
        U2D(q) = priceScale*d2res.upper;

        L3P(q) = priceScale*p3res.lower;
        L3D(q) = priceScale*d3res.lower;
        U3P(q) = priceScale*p3res.upper;
        U3D(q) = priceScale*d3res.upper;

        %% Signed primal-dual gaps in the theoretically correct direction
        % lower: primal(min) - dual(max) >= 0
        % upper: dual(min) - primal(max) >= 0
        GapL2(q) = L2P(q)-L2D(q);
        GapU2(q) = U2D(q)-U2P(q);
        GapL3(q) = L3P(q)-L3D(q);
        GapU3(q) = U3D(q)-U3P(q);

        MaxGap2(q) = max(abs([GapL2(q),GapU2(q)]));
        MaxGap3(q) = max(abs([GapL3(q),GapU3(q)]));

        %% Raw differences from primal solver outputs
        RawLowerTightening(q) = L3P(q)-L2P(q);
        RawUpperTightening(q) = U2P(q)-U3P(q);

        width2P = U2P(q)-L2P(q);
        width3P = U3P(q)-L3P(q);
        RawWidthReduction(q) = width2P-width3P;

        if abs(width2P) > cfg3.displayTolerance
            RawRelativeWidthReductionPct(q) = ...
                100*RawWidthReduction(q)/width2P;
        end

        %% Primal-dual brackets
        % Lower optimum L is bracketed approximately by [L_D,L_P].
        % Upper optimum U is bracketed approximately by [U_P,U_D].
        LowerTighteningBracketLo(q) = L3D(q)-L2P(q);
        LowerTighteningBracketHi(q) = L3P(q)-L2D(q);

        UpperTighteningBracketLo(q) = U2P(q)-U3D(q);
        UpperTighteningBracketHi(q) = U2D(q)-U3P(q);

        W2lo = U2P(q)-L2P(q);
        W2hi = U2D(q)-L2D(q);
        W3lo = U3P(q)-L3P(q);
        W3hi = U3D(q)-L3D(q);

        WidthReductionBracketLo(q) = W2lo-W3hi;
        WidthReductionBracketHi(q) = W2hi-W3lo;

        %% Feasibility diagnostics
        MaxPrimalFeasResidual2(q) = max([ ...
            p2res.diagnostics.lower.rowResidual, ...
            p2res.diagnostics.lower.colResidual, ...
            p2res.diagnostics.lower.martingaleResidual, ...
            p2res.diagnostics.lower.massResidual, ...
            p2res.diagnostics.upper.rowResidual, ...
            p2res.diagnostics.upper.colResidual, ...
            p2res.diagnostics.upper.martingaleResidual, ...
            p2res.diagnostics.upper.massResidual]);

        MaxPrimalFeasResidual3(q) = max([ ...
            p3res.diagnostics.lower.mu1Residual, ...
            p3res.diagnostics.lower.mu2Residual, ...
            p3res.diagnostics.lower.mu3Residual, ...
            p3res.diagnostics.lower.martingale12Residual, ...
            p3res.diagnostics.lower.martingale23Residual, ...
            p3res.diagnostics.lower.massResidual, ...
            p3res.diagnostics.upper.mu1Residual, ...
            p3res.diagnostics.upper.mu2Residual, ...
            p3res.diagnostics.upper.mu3Residual, ...
            p3res.diagnostics.upper.martingale12Residual, ...
            p3res.diagnostics.upper.martingale23Residual, ...
            p3res.diagnostics.upper.massResidual]);

        MaxDualHedgeViolation2(q) = max([0, ...
            d2res.diagnostics.lowerMaxHedgeViolation, ...
            d2res.diagnostics.upperMaxHedgeViolation]);

        MaxDualHedgeViolation3(q) = max([0, ...
            d3res.diagnostics.lowerMaxHedgeViolation, ...
            d3res.diagnostics.upperMaxHedgeViolation]);

        %% CVX-reported solver tolerances (supplementary only)
        PrimalSlvTol2(q) = max([p2res.slvTolLower,p2res.slvTolUpper]);
        DualSlvTol2(q) = max([d2res.slvTolLower,d2res.slvTolUpper]);
        PrimalSlvTol3(q) = max([p3res.slvTolLower,p3res.slvTolUpper]);
        DualSlvTol3(q) = max([d3res.slvTolLower,d3res.slvTolUpper]);

        %% Status strings
        StatusL2P(q) = string(p2res.statusLower);
        StatusU2P(q) = string(p2res.statusUpper);
        StatusL2D(q) = string(d2res.statusLower);
        StatusU2D(q) = string(d2res.statusUpper);
        StatusL3P(q) = string(p3res.statusLower);
        StatusU3P(q) = string(p3res.statusUpper);
        StatusL3D(q) = string(d3res.statusLower);
        StatusU3D(q) = string(d3res.statusUpper);

        RunStatus(q) = "Completed";

        fprintf('2M max primal-dual gap ($): %.6e\n',MaxGap2(q));
        fprintf('3M max primal-dual gap ($): %.6e\n',MaxGap3(q));
        fprintf('Raw width reduction ($)   : %.6e\n',RawWidthReduction(q));
        fprintf('Width-reduction bracket   : [%.6e, %.6e]\n', ...
            WidthReductionBracketLo(q),WidthReductionBracketHi(q));

    catch ME
        RunStatus(q) = "Failed";
        fprintf('\nValidation failed at K = %.2f:\n%s\n',K(q),ME.message);
        cvx_clear;
        cvx_solver sdpt3;
    end
end

%% Results table

validation = table( ...
    K,Kmot, ...
    L2P,L2D,GapL2,U2P,U2D,GapU2,MaxGap2, ...
    L3P,L3D,GapL3,U3P,U3D,GapU3,MaxGap3, ...
    RawLowerTightening,RawUpperTightening, ...
    RawWidthReduction,RawRelativeWidthReductionPct, ...
    LowerTighteningBracketLo,LowerTighteningBracketHi, ...
    UpperTighteningBracketLo,UpperTighteningBracketHi, ...
    WidthReductionBracketLo,WidthReductionBracketHi, ...
    MaxPrimalFeasResidual2,MaxPrimalFeasResidual3, ...
    MaxDualHedgeViolation2,MaxDualHedgeViolation3, ...
    PrimalSlvTol2,DualSlvTol2,PrimalSlvTol3,DualSlvTol3, ...
    StatusL2P,StatusU2P,StatusL2D,StatusU2D, ...
    StatusL3P,StatusU3P,StatusL3D,StatusU3D,RunStatus);

writetable(validation,fullfile(cfg3.outputDir, ...
    'market_optimality_validation_main17.csv'));

save(fullfile(cfg3.outputDir,'market_optimality_validation_main17.mat'), ...
    'validation','compression','priceScale','Kconversion');

fprintf('\n\nFULL OBJECTIVE-LEVEL VALIDATION TABLE\n');
disp(validation);

%% K=1 summary

[~,idx1] = min(abs(K-1));

fprintf('\n============================================================\n');
fprintf('K = 1 OBJECTIVE-LEVEL SUMMARY\n');
fprintf('============================================================\n');
fprintf('2M lower primal / dual ($): %.10f / %.10f\n',L2P(idx1),L2D(idx1));
fprintf('2M lower gap ($)          : %.6e\n',GapL2(idx1));
fprintf('2M upper primal / dual ($): %.10f / %.10f\n',U2P(idx1),U2D(idx1));
fprintf('2M upper gap ($)          : %.6e\n',GapU2(idx1));
fprintf('3M lower primal / dual ($): %.10f / %.10f\n',L3P(idx1),L3D(idx1));
fprintf('3M lower gap ($)          : %.6e\n',GapL3(idx1));
fprintf('3M upper primal / dual ($): %.10f / %.10f\n',U3P(idx1),U3D(idx1));
fprintf('3M upper gap ($)          : %.6e\n',GapU3(idx1));

fprintf('\nRaw primal width reduction ($): %.10f\n',RawWidthReduction(idx1));
fprintf('Raw relative reduction (%%)    : %.8f\n', ...
    RawRelativeWidthReductionPct(idx1));
fprintf('Approx width-reduction bracket : [%.10f, %.10f]\n', ...
    WidthReductionBracketLo(idx1),WidthReductionBracketHi(idx1));

fprintf('\nLargest gaps over full K grid\n');
validGap2 = MaxGap2(~isnan(MaxGap2));
validGap3 = MaxGap3(~isnan(MaxGap3));
if isempty(validGap2), validGap2 = NaN; end
if isempty(validGap3), validGap3 = NaN; end
fprintf('2M max |gap| ($): %.6e\n',max(validGap2));
fprintf('3M max |gap| ($): %.6e\n',max(validGap3));

fprintf('\nIMPORTANT INTERPRETATION\n');
fprintf(['Feasibility residuals and cvx_slvtol are supplementary diagnostics.\n' ...
         'For the dissertation, compare the explicit primal-dual gaps and\n' ...
         'the width-reduction bracket with the sub-cent differences being\n' ...
         'discussed. If the bracket contains zero, the sign/magnitude of\n' ...
         'that tiny tightening is not numerically resolved by these runs.\n']);

fprintf('\nObjective-level validation complete.\n');
