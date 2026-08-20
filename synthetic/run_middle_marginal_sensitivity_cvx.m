clear;
clc;
close all;
cvx_clear;
cvx_solver sdpt3;

[x,a,z,b,s0] = build_synthetic_marginals();

lambdaGrid = [0 0.25 0.50 0.75 1.00];
Kgrid = (0.50:0.05:1.50)';

qN = numel(Kgrid);
rN = numel(lambdaGrid);

%% Two-marginal reference curve

W2 = zeros(qN,1);

for q = 1:qN
    r2 = mot_bounds_cvx(x,a,z,b,Kgrid(q));
    W2(q) = r2.width;
end

%% Three-marginal width curves

W3 = zeros(qN,rN);

for r = 1:rN

    lambda = lambdaGrid(r);

    [y,c] = build_middle_marginal_lambda(lambda);

    fprintf('Running lambda = %.2f with %d intermediate support points\n', ...
            lambda,numel(y));

    for q = 1:qN
        r3 = mot_bounds_3marginal_cvx( ...
            x,a,y,c,z,b,Kgrid(q));

        W3(q,r) = r3.width;
    end
end

%% Summary at K = 1

[~,idx1] = min(abs(Kgrid-1));

WidthAtK1 = W3(idx1,:)';
ReductionAtK1 = W2(idx1)-WidthAtK1;

RelativeReductionAtK1 = ...
    100*ReductionAtK1/W2(idx1);

summary = table( ...
    lambdaGrid', ...
    WidthAtK1, ...
    ReductionAtK1, ...
    RelativeReductionAtK1, ...
    'VariableNames',{ ...
    'Lambda','WidthAtK1','AbsoluteReductionAtK1', ...
    'RelativeReductionPercentAtK1'});

disp(summary);

writetable(summary, ...
    'middle_marginal_sensitivity_summary.csv');

%% Results

fullTable = table(Kgrid,W2, ...
    'VariableNames',{'K','TwoMarginalWidth'});

for r = 1:rN
    name = sprintf('Width_lambda_%03d',round(100*lambdaGrid(r)));
    fullTable.(name) = W3(:,r);
end

writetable(fullTable, ...
    'middle_marginal_sensitivity_all.csv');

%% Plot

figure('Color','w');

plot(Kgrid,W2/s0,'--','LineWidth',1.5, ...
    'DisplayName','Two marginals');

hold on;

for r = 1:rN
    plot(Kgrid,W3(:,r)/s0,'LineWidth',1.1, ...
        'DisplayName',sprintf('\\lambda = %.2f',lambdaGrid(r)));
end

xlabel('Strike multiplier K');
ylabel('(U-L)/s_0');

grid on;
legend('Location','best');

title('Sensitivity to the intermediate marginal');

exportgraphics(gcf, ...
    'middle_marginal_sensitivity.pdf', ...
    'ContentType','vector');
