clear;
clc;
close all;

% ===============================================================
% OBJECTIVE-LEVEL K=1 VALIDATION ACROSS 17/21/23-POINT SUPPORTS
%
% This is the primal-dual counterpart of run_market_mot_resolution_check.
% It is intended to determine whether the tiny computed width reductions
% are larger than the numerical objective gaps.
% ===============================================================

cvx_clear;
cvx_solver sdpt3;

cfg3 = market_mot_config();

if ~exist(cfg3.outputDir,'dir')
    mkdir(cfg3.outputDir);
end

marketMarginals = load_market_marginals();

fineSupport = marketMarginals.supportNormalized(:);
fineProb = marketMarginals.probabilities;
F = marketMarginals.forwardEstimates(:);
beta = marketMarginals.discountFactors(:);

F1 = F(1);
F3 = F(3);
beta3 = beta(3);

K = 1.0;
Kmot = K*F1/F3;
priceScale = beta3*F3;

supportCells = { ...
    cfg3.mainSupportNormalized, ...
    cfg3.checkSupport21, ...
    cfg3.checkSupport23};
labels = ["Main17","Medium21","Finer23"];

nCase = numel(supportCells);

NumSupport = zeros(nCase,1);
CompressionRMSE = nan(nCase,1);

L2P = nan(nCase,1); L2D = nan(nCase,1);
U2P = nan(nCase,1); U2D = nan(nCase,1);
L3P = nan(nCase,1); L3D = nan(nCase,1);
U3P = nan(nCase,1); U3D = nan(nCase,1);

GapL2 = nan(nCase,1); GapU2 = nan(nCase,1);
GapL3 = nan(nCase,1); GapU3 = nan(nCase,1);
MaxGap2 = nan(nCase,1); MaxGap3 = nan(nCase,1);

RawWidth2 = nan(nCase,1);
RawWidth3 = nan(nCase,1);
RawWidthReduction = nan(nCase,1);
RawRelativeReductionPct = nan(nCase,1);
WidthReductionBracketLo = nan(nCase,1);
WidthReductionBracketHi = nan(nCase,1);

MaxPrimalFeasResidual2 = nan(nCase,1);
MaxPrimalFeasResidual3 = nan(nCase,1);
MaxDualHedgeViolation2 = nan(nCase,1);
MaxDualHedgeViolation3 = nan(nCase,1);

Status = strings(nCase,1);

for r = 1:nCase

    fprintf('\n================================================\n');
    fprintf('Objective validation: %s support\n',labels(r));
    fprintf('================================================\n');

    try
        comp = compress_market_marginals_cvx( ...
            fineSupport,fineProb,supportCells{r});

        support = comp.support(:);
        P = comp.probabilities;
        p1 = P(1,:)';
        p2 = P(2,:)';
        p3 = P(3,:)';

        NumSupport(r) = numel(support);
        CompressionRMSE(r) = comp.rmseCall;

        p2res = mot_bounds_cvx(support,p1,support,p3,Kmot);
        d2res = mot_dual_cvx(support,p1,support,p3,Kmot);

        p3res = mot_bounds_3marginal_cvx( ...
            support,p1,support,p2,support,p3,Kmot);
        d3res = mot_dual_3marginal_cvx( ...
            support,p1,support,p2,support,p3,Kmot);

        L2P(r)=priceScale*p2res.lower; L2D(r)=priceScale*d2res.lower;
        U2P(r)=priceScale*p2res.upper; U2D(r)=priceScale*d2res.upper;
        L3P(r)=priceScale*p3res.lower; L3D(r)=priceScale*d3res.lower;
        U3P(r)=priceScale*p3res.upper; U3D(r)=priceScale*d3res.upper;

        GapL2(r)=L2P(r)-L2D(r);
        GapU2(r)=U2D(r)-U2P(r);
        GapL3(r)=L3P(r)-L3D(r);
        GapU3(r)=U3D(r)-U3P(r);

        MaxGap2(r)=max(abs([GapL2(r),GapU2(r)]));
        MaxGap3(r)=max(abs([GapL3(r),GapU3(r)]));

        RawWidth2(r)=U2P(r)-L2P(r);
        RawWidth3(r)=U3P(r)-L3P(r);
        RawWidthReduction(r)=RawWidth2(r)-RawWidth3(r);
        RawRelativeReductionPct(r)=100*RawWidthReduction(r)/RawWidth2(r);

        W2lo = U2P(r)-L2P(r);
        W2hi = U2D(r)-L2D(r);
        W3lo = U3P(r)-L3P(r);
        W3hi = U3D(r)-L3D(r);

        WidthReductionBracketLo(r)=W2lo-W3hi;
        WidthReductionBracketHi(r)=W2hi-W3lo;

        MaxPrimalFeasResidual2(r) = max([ ...
            p2res.diagnostics.lower.rowResidual, ...
            p2res.diagnostics.lower.colResidual, ...
            p2res.diagnostics.lower.martingaleResidual, ...
            p2res.diagnostics.lower.massResidual, ...
            p2res.diagnostics.upper.rowResidual, ...
            p2res.diagnostics.upper.colResidual, ...
            p2res.diagnostics.upper.martingaleResidual, ...
            p2res.diagnostics.upper.massResidual]);

        MaxPrimalFeasResidual3(r) = max([ ...
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

        MaxDualHedgeViolation2(r)=max([0, ...
            d2res.diagnostics.lowerMaxHedgeViolation, ...
            d2res.diagnostics.upperMaxHedgeViolation]);
        MaxDualHedgeViolation3(r)=max([0, ...
            d3res.diagnostics.lowerMaxHedgeViolation, ...
            d3res.diagnostics.upperMaxHedgeViolation]);

        Status(r)="Completed";

        fprintf('Raw width reduction ($): %.8f\n',RawWidthReduction(r));
        fprintf('2M max primal-dual gap : %.3e\n',MaxGap2(r));
        fprintf('3M max primal-dual gap : %.3e\n',MaxGap3(r));
        fprintf('Reduction bracket       : [%.8f, %.8f]\n', ...
            WidthReductionBracketLo(r),WidthReductionBracketHi(r));

    catch ME
        Status(r)="Failed";
        fprintf('\n%s validation failed:\n%s\n',labels(r),ME.message);
        cvx_clear;
        cvx_solver sdpt3;
    end
end

summary = table(labels',NumSupport,CompressionRMSE, ...
    L2P,L2D,GapL2,U2P,U2D,GapU2,MaxGap2, ...
    L3P,L3D,GapL3,U3P,U3D,GapU3,MaxGap3, ...
    RawWidth2,RawWidth3,RawWidthReduction,RawRelativeReductionPct, ...
    WidthReductionBracketLo,WidthReductionBracketHi, ...
    MaxPrimalFeasResidual2,MaxPrimalFeasResidual3, ...
    MaxDualHedgeViolation2,MaxDualHedgeViolation3,Status, ...
    'VariableNames',{ ...
    'Grid','NumSupport','CompressionRMSE', ...
    'L2Primal','L2Dual','GapL2','U2Primal','U2Dual','GapU2','MaxGap2', ...
    'L3Primal','L3Dual','GapL3','U3Primal','U3Dual','GapU3','MaxGap3', ...
    'RawWidth2','RawWidth3','RawWidthReduction','RawRelativeReductionPct', ...
    'WidthReductionBracketLo','WidthReductionBracketHi', ...
    'MaxPrimalFeasResidual2','MaxPrimalFeasResidual3', ...
    'MaxDualHedgeViolation2','MaxDualHedgeViolation3','Status'});

fprintf('\n\nK=1 OBJECTIVE-LEVEL SUPPORT SUMMARY\n');
disp(summary);

writetable(summary,fullfile(cfg3.outputDir, ...
    'market_optimality_validation_supports_K1.csv'));

save(fullfile(cfg3.outputDir, ...
    'market_optimality_validation_supports_K1.mat'), ...
    'summary','priceScale','Kmot');

fprintf('\nSupport-level objective validation complete.\n');
