clear;
clc;
close all;
cvx_clear;
cvx_solver sdpt3;

Ngrid = [4 9 16 25];

Kgrid = (0.50:0.05:1.50)';

sigma1 = 10;
noiseSize = 20;

allResults = table();
summaryRows = [];

figure('Color','w');
hold on;

for r = 1:numel(Ngrid)

    N = Ngrid(r);

    %% Construct marginals

    [x,a,z,b,s0] = synthetic_binomial_pair(N,sigma1,noiseSize);

    cx = check_convex_order(x,a,z,b);

    if ~cx.inConvexOrder
        error('Resolution N=%d failed convex order.',N);
    end

    qN = numel(Kgrid);

    L = zeros(qN,1);
    U = zeros(qN,1);
    W = zeros(qN,1);
    Wspot = zeros(qN,1);
    feas = zeros(qN,1);

    %% Solve strike values

    for q = 1:qN

        K = Kgrid(q);

        primal = mot_bounds_cvx(x,a,z,b,K);

        L(q) = primal.lower;
        U(q) = primal.upper;
        W(q) = primal.width;
        Wspot(q) = W(q)/s0;

        feas(q) = max([ ...
            primal.diagnostics.lower.rowResidual, ...
            primal.diagnostics.lower.colResidual, ...
            primal.diagnostics.lower.martingaleResidual, ...
            primal.diagnostics.lower.massResidual, ...
            primal.diagnostics.upper.rowResidual, ...
            primal.diagnostics.upper.colResidual, ...
            primal.diagnostics.upper.martingaleResidual, ...
            primal.diagnostics.upper.massResidual]);
    end

    %% Clean noise

    tolDisplay = 1e-7;
    L(abs(L) < tolDisplay) = 0;
    U(abs(U) < tolDisplay) = 0;
    W(abs(W) < tolDisplay) = 0;
    Wspot(abs(Wspot) < tolDisplay) = 0;

    %% Results

    m = numel(x);
    n = numel(z);

    T = table( ...
        repmat(N,qN,1), ...
        repmat(m,qN,1), ...
        repmat(n,qN,1), ...
        Kgrid,L,U,W,Wspot,feas, ...
        'VariableNames',{ ...
        'N','m','n','K','Lower','Upper','Width','WidthSpot', ...
        'MaxFeasResidual'});

    if isempty(allResults)
        allResults = T;
    else
        allResults = [allResults; T]; 
    end

    [~,idx1] = min(abs(Kgrid-1));

    summaryRows = [summaryRows; ...
        N,m,n,m*n,max(feas),W(idx1)]; 

    %% Plot width curve

    plot(Kgrid,Wspot,'LineWidth',1.2, ...
        'DisplayName',sprintf('N=%d (m=%d,n=%d)',N,m,n));
end


xlabel('Strike multiplier K');
ylabel('(U-L)/s_0');
grid on;
legend('Location','best');
title('Support-resolution sensitivity');

exportgraphics(gcf,'resolution_width_comparison_cvx.pdf','ContentType','vector');

writetable(allResults,'resolution_study_all_cvx.csv');

%% Summary

summary = array2table(summaryRows, ...
    'VariableNames',{ ...
    'N','m','n','NumVariables','MaxFeasResidual','WidthAtK1'});

disp(summary);

writetable(summary,'resolution_study_summary_cvx.csv');
