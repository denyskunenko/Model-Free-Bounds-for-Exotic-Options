clear;
clc;
close all;
cvx_clear;
cvx_solver sdpt3;

[x,a,y,c,z,b,s0] = build_three_marginals();

cx12 = check_convex_order(x,a,y,c);
cx23 = check_convex_order(y,c,z,b);

fprintf('\nThree-marginal input checks\n');
fprintf('Mean(mu1) = %.12f\n',a'*x);
fprintf('Mean(mu2) = %.12f\n',c'*y);
fprintf('Mean(mu3) = %.12f\n',b'*z);

fprintf('Max convex-order violation mu1 <=cx mu2: %.3e\n', ...
        cx12.maxViolation);

fprintf('Max convex-order violation mu2 <=cx mu3: %.3e\n\n', ...
        cx23.maxViolation);

Kgrid = (0.50:0.05:1.50)';
qN = numel(Kgrid);

L2 = zeros(qN,1);
U2 = zeros(qN,1);
W2 = zeros(qN,1);

L3 = zeros(qN,1);
U3 = zeros(qN,1);
W3 = zeros(qN,1);

lowerTightening = zeros(qN,1);
upperTightening = zeros(qN,1);
widthReduction = zeros(qN,1);
relativeReduction = nan(qN,1);

maxFeasResidual3 = zeros(qN,1);

statusLower3 = strings(qN,1);
statusUpper3 = strings(qN,1);

tolNesting = 1e-6;

%% Solve for every K


for q = 1:qN

    K = Kgrid(q);

    r2 = mot_bounds_cvx(x,a,z,b,K);

    r3 = mot_bounds_3marginal_cvx(x,a,y,c,z,b,K);

    L2(q) = r2.lower;
    U2(q) = r2.upper;
    W2(q) = r2.width;

    L3(q) = r3.lower;
    U3(q) = r3.upper;
    W3(q) = r3.width;

    lowerTightening(q) = L3(q)-L2(q);
    upperTightening(q) = U2(q)-U3(q);
    widthReduction(q) = W2(q)-W3(q);

    if W2(q) > 1e-7
        relativeReduction(q) = widthReduction(q)/W2(q);
    end

    if L3(q) < L2(q)-tolNesting
        warning('Lower-bound nesting failed numerically at K=%.2f.',K);
    end

    if U3(q) > U2(q)+tolNesting
        warning('Upper-bound nesting failed numerically at K=%.2f.',K);
    end

    maxFeasResidual3(q) = max([ ...
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

    statusLower3(q) = string(r3.statusLower);
    statusUpper3(q) = string(r3.statusUpper);
end

%% Clean noise

tolDisplay = 1e-7;

L2(abs(L2) < tolDisplay) = 0;
U2(abs(U2) < tolDisplay) = 0;
W2(abs(W2) < tolDisplay) = 0;

L3(abs(L3) < tolDisplay) = 0;
U3(abs(U3) < tolDisplay) = 0;
W3(abs(W3) < tolDisplay) = 0;

lowerTightening(abs(lowerTightening) < tolDisplay) = 0;
upperTightening(abs(upperTightening) < tolDisplay) = 0;
widthReduction(abs(widthReduction) < tolDisplay) = 0;

%% RESULTS

results = table( ...
    Kgrid, ...
    L2,U2,W2, ...
    L3,U3,W3, ...
    lowerTightening,upperTightening,widthReduction,relativeReduction, ...
    maxFeasResidual3,statusLower3,statusUpper3, ...
    'VariableNames',{ ...
    'K','Lower2','Upper2','Width2', ...
    'Lower3','Upper3','Width3', ...
    'LowerTightening','UpperTightening', ...
    'WidthReduction','RelativeWidthReduction', ...
    'MaxFeasResidual3','StatusLower3','StatusUpper3'});

disp(results);

writetable(results,'three_marginal_comparison_results.csv');

[~,idx1] = min(abs(Kgrid-1));

fprintf('\nComparison at K = 1\n');

fprintf('Two-marginal lower  = %.10f\n',L2(idx1));
fprintf('Three-marginal lower= %.10f\n',L3(idx1));

fprintf('Two-marginal upper  = %.10f\n',U2(idx1));
fprintf('Three-marginal upper= %.10f\n',U3(idx1));

fprintf('Two-marginal width  = %.10f\n',W2(idx1));
fprintf('Three-marginal width= %.10f\n',W3(idx1));

fprintf('Absolute width reduction = %.10f\n',widthReduction(idx1));

if ~isnan(relativeReduction(idx1))
    fprintf('Relative width reduction = %.2f%%\n', ...
            100*relativeReduction(idx1));
end

fprintf('Maximum three-marginal residual over all K = %.3e\n', ...
        max(maxFeasResidual3));

%% Figure 1: lower and upper bounds

figure('Color','w');

plot(Kgrid,L2,'-o','LineWidth',1.1,'MarkerSize',4);
hold on;
plot(Kgrid,U2,'-s','LineWidth',1.1,'MarkerSize',4);

plot(Kgrid,L3,'--o','LineWidth',1.3,'MarkerSize',4);
plot(Kgrid,U3,'--s','LineWidth',1.3,'MarkerSize',4);

xlabel('Strike multiplier K');
ylabel('Model-independent price');

legend( ...
    'Two-marginal lower', ...
    'Two-marginal upper', ...
    'Three-marginal lower', ...
    'Three-marginal upper', ...
    'Location','best');

grid on;
title('Two- versus three-marginal model-independent bounds');

exportgraphics(gcf, ...
    'two_vs_three_marginal_bounds.pdf', ...
    'ContentType','vector');


%% Figure 2: interval widths

figure('Color','w');

plot(Kgrid,W2/s0,'-o','LineWidth',1.2,'MarkerSize',4);
hold on;

plot(Kgrid,W3/s0,'-s','LineWidth',1.2,'MarkerSize',4);

xlabel('Strike multiplier K');
ylabel('(U-L)/s_0');

legend( ...
    'Two marginals', ...
    'Three marginals', ...
    'Location','best');

grid on;
title('Effect of the intermediate marginal on bound width');

exportgraphics(gcf, ...
    'two_vs_three_marginal_widths.pdf', ...
    'ContentType','vector');
