clear;
clc;
close all;

% ===============================================================
% K=1 SUPPORT-COMPRESSION ROBUSTNESS CHECK
%
% Progressive check using 17, 21 and 23 support points.
%
% The previous 31-point three-marginal CVX problem was numerically too
% large for SDPT3 on this setup. These grids still provide a meaningful
% refinement while keeping the joint LP manageable.
% ===============================================================

cvx_clear;

cfg3 = market_mot_config();

if ~exist(cfg3.outputDir,'dir')
    mkdir(cfg3.outputDir);
end

marketMarginals = load_market_marginals();

S0 = marketMarginals.spot;

fineSupport = marketMarginals.supportNormalized(:);
fineProb = marketMarginals.probabilities;

F = marketMarginals.forwardEstimates(:);
D = marketMarginals.discountFactors(:);

F1 = F(1);
F3 = F(3);
D3 = D(3);

K = 1.0;
Kmot = K*F1/F3;

% Correct normalized-state dollar conversion:
priceScale = D3*F3;

supportCells = { ...
    cfg3.mainSupportNormalized, ...
    cfg3.checkSupport21, ...
    cfg3.checkSupport23};

labels = ["Main17","Medium21","Finer23"];

nCase = numel(supportCells);

NumSupport = zeros(nCase,1);

CompressionRMSE = nan(nCase,1);
CompressionMaxError = nan(nCase,1);

Lower2 = nan(nCase,1);
Upper2 = nan(nCase,1);
Width2 = nan(nCase,1);

Lower3 = nan(nCase,1);
Upper3 = nan(nCase,1);
Width3 = nan(nCase,1);

AbsoluteReduction = nan(nCase,1);
RelativeReductionPercent = nan(nCase,1);

MaxFeasResidual2 = nan(nCase,1);
MaxFeasResidual3 = nan(nCase,1);

Status = strings(nCase,1);

for r = 1:nCase

    fprintf('\n================================================\n');
    fprintf('Running %s support check...\n',labels(r));
    fprintf('================================================\n');

    try

        comp = compress_market_marginals_cvx( ...
            fineSupport, ...
            fineProb, ...
            supportCells{r});

        supportMOT = comp.support;
        P = comp.probabilities;

        p1 = P(1,:)';
        p2 = P(2,:)';
        p3 = P(3,:)';

        NumSupport(r) = numel(supportMOT);

        CompressionRMSE(r) = comp.rmseCall;
        CompressionMaxError(r) = comp.maxAbsCallError;

        fprintf('Support points        : %d\n',NumSupport(r));
        fprintf('Compression RMSE      : %.3e\n',CompressionRMSE(r));
        fprintf('Compression max error : %.3e\n',CompressionMaxError(r));

        %% Two-marginal problem

        r2 = mot_bounds_cvx( ...
            supportMOT,p1, ...
            supportMOT,p3, ...
            Kmot);

        %% Three-marginal problem

        r3 = mot_bounds_3marginal_cvx( ...
            supportMOT,p1, ...
            supportMOT,p2, ...
            supportMOT,p3, ...
            Kmot);

        %% Convert normalized MOT values to time-zero dollar prices

        Lower2(r) = priceScale*r2.lower;
        Upper2(r) = priceScale*r2.upper;
        Width2(r) = Upper2(r)-Lower2(r);

        Lower3(r) = priceScale*r3.lower;
        Upper3(r) = priceScale*r3.upper;
        Width3(r) = Upper3(r)-Lower3(r);

        AbsoluteReduction(r) = Width2(r)-Width3(r);

        RelativeReductionPercent(r) = ...
            100*AbsoluteReduction(r)/Width2(r);

        %% Numerical residuals

        MaxFeasResidual2(r) = max([ ...
            r2.diagnostics.lower.rowResidual, ...
            r2.diagnostics.lower.colResidual, ...
            r2.diagnostics.lower.martingaleResidual, ...
            r2.diagnostics.lower.massResidual, ...
            r2.diagnostics.upper.rowResidual, ...
            r2.diagnostics.upper.colResidual, ...
            r2.diagnostics.upper.martingaleResidual, ...
            r2.diagnostics.upper.massResidual]);

        MaxFeasResidual3(r) = max([ ...
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

        Status(r) = "Completed";

        fprintf('Two-marginal width    : %.8f\n',Width2(r));
        fprintf('Three-marginal width  : %.8f\n',Width3(r));
        fprintf('Raw computed reduction: %.4f%%\n', ...
                RelativeReductionPercent(r));

    catch ME

        Status(r) = "Failed";

        fprintf('\n%s grid failed:\n%s\n',labels(r),ME.message);
        fprintf('Continuing to any remaining grid checks.\n');

        % Reset any stale CVX model before continuing.
        cvx_clear;
    end
end

%% Summary table

summary = table( ...
    labels', ...
    NumSupport, ...
    CompressionRMSE, ...
    CompressionMaxError, ...
    Lower2,Upper2,Width2, ...
    Lower3,Upper3,Width3, ...
    AbsoluteReduction, ...
    RelativeReductionPercent, ...
    MaxFeasResidual2, ...
    MaxFeasResidual3, ...
    Status, ...
    'VariableNames',{ ...
    'Grid','NumSupport', ...
    'CompressionRMSE','CompressionMaxError', ...
    'Lower2','Upper2','Width2', ...
    'Lower3','Upper3','Width3', ...
    'AbsoluteWidthReduction', ...
    'RelativeReductionPercent', ...
    'MaxFeasResidual2', ...
    'MaxFeasResidual3', ...
    'Status'});

fprintf('\n\nK=1 SUPPORT ROBUSTNESS SUMMARY\n');
disp(summary);

writetable(summary, ...
    fullfile(cfg3.outputDir, ...
    'market_support_resolution_K1.csv'));

%% Simple comparison plot for successful runs

ok = Status == "Completed";

if any(ok)

    figure('Color','w');

    plot(NumSupport(ok),Width2(ok),'-o', ...
        'LineWidth',1.2,'MarkerSize',5, ...
        'DisplayName','Two marginals');

    hold on;

    plot(NumSupport(ok),Width3(ok),'-s', ...
        'LineWidth',1.2,'MarkerSize',5, ...
        'DisplayName','Three marginals');

    xlabel('Number of support points');
    ylabel('K=1 interval width');
    grid on;
    legend('Location','best');

    title('Market MOT support-resolution check');

    exportgraphics(gcf, ...
        fullfile(cfg3.outputDir, ...
        'market_support_resolution_K1.pdf'), ...
        'ContentType','vector');
end

fprintf('\nK=1 progressive robustness check complete.\n');
