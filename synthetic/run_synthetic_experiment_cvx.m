clear;
clc;
close all;
cvx_clear;
cvx_solver sdpt3;

[x,a,z,b,s0] = build_synthetic_marginals();

%% Input checks

cx = check_convex_order(x,a,z,b);

fprintf('\nInput checks\n');
fprintf('Mean(mu1) = %.12f\n',cx.mean1);
fprintf('Mean(mu3) = %.12f\n',cx.mean3);
fprintf('Maximum convex-order violation = %.3e\n\n',cx.maxViolation);

%% Strike grid

Kgrid = (0.50:0.05:1.50)';
qN = numel(Kgrid);

%% Storage

L = zeros(qN,1);
U = zeros(qN,1);
W = zeros(qN,1);
Wspot = zeros(qN,1);

rowRes = zeros(qN,1);
colRes = zeros(qN,1);
martRes = zeros(qN,1);
massRes = zeros(qN,1);

statusLower = strings(qN,1);
statusUpper = strings(qN,1);

%% Solve all strikes

for q = 1:qN
    K = Kgrid(q);

    r = mot_bounds_cvx(x,a,z,b,K);

    L(q) = r.lower;
    U(q) = r.upper;
    W(q) = r.width;
    Wspot(q) = W(q)/s0;

    rowRes(q) = max( ...
        r.diagnostics.lower.rowResidual, ...
        r.diagnostics.upper.rowResidual);

    colRes(q) = max( ...
        r.diagnostics.lower.colResidual, ...
        r.diagnostics.upper.colResidual);

    martRes(q) = max( ...
        r.diagnostics.lower.martingaleResidual, ...
        r.diagnostics.upper.martingaleResidual);

    massRes(q) = max( ...
        r.diagnostics.lower.massResidual, ...
        r.diagnostics.upper.massResidual);

    statusLower(q) = string(r.statusLower);
    statusUpper(q) = string(r.statusUpper);
end

%% Clean noise

tolDisplay = 1e-7;

L(abs(L) < tolDisplay) = 0;
U(abs(U) < tolDisplay) = 0;
W(abs(W) < tolDisplay) = 0;

%% Results

results = table( ...
    Kgrid,L,U,W,Wspot, ...
    rowRes,colRes,martRes,massRes, ...
    statusLower,statusUpper, ...
    'VariableNames',{ ...
    'Kgrid','Lower','Upper','Width','WidthSpot', ...
    'MaxRowResidual','MaxColResidual', ...
    'MaxMartResidual','MaxMassResidual', ...
    'StatusLower','StatusUpper'});

disp(results);

writetable(results,'synthetic_bounds_results_cvx.csv');

%% Maximum residuals

fprintf('\nMaximum numerical residuals across the strike grid\n');
fprintf('Row marginal     : %.3e\n',max(rowRes));
fprintf('Column marginal  : %.3e\n',max(colRes));
fprintf('Martingale       : %.3e\n',max(martRes));
fprintf('Total mass       : %.3e\n',max(massRes));

%% Primal-dual validation at K = 1

K0 = 1.0;

primal0 = mot_bounds_cvx(x,a,z,b,K0);
dual0 = mot_dual_cvx(x,a,z,b,K0);

fprintf('\nBenchmark and primal-dual validation at K = 1\n');
fprintf('Primal lower = %.10f\n',primal0.lower);
fprintf('Dual lower   = %.10f\n',dual0.lower);
fprintf('Lower gap    = %.3e\n',abs(primal0.lower-dual0.lower));

fprintf('Primal upper = %.10f\n',primal0.upper);
fprintf('Dual upper   = %.10f\n',dual0.upper);
fprintf('Upper gap    = %.3e\n',abs(primal0.upper-dual0.upper));

fprintf('Width        = %.10f\n',primal0.width);

%% Figure 1: bounds

figure('Color','w');

plot(Kgrid,L,'-o','LineWidth',1.2,'MarkerSize',4);
hold on;
plot(Kgrid,U,'-s','LineWidth',1.2,'MarkerSize',4);

xlabel('Strike multiplier K');
ylabel('Model-independent price');
legend('Lower bound','Upper bound','Location','best');
grid on;
title('Two-marginal model-independent bounds');

exportgraphics(gcf,'synthetic_bounds_cvx.pdf','ContentType','vector');

%% Figure 2: normalised width

figure('Color','w');

plot(Kgrid,Wspot,'-o','LineWidth',1.2,'MarkerSize',4);

xlabel('Strike multiplier K');
ylabel('(U-L)/s_0');
grid on;
title('Normalised width of the two-marginal interval');

exportgraphics(gcf,'synthetic_width_cvx.pdf','ContentType','vector');
