function [y,c] = build_middle_marginal_lambda(lambda)

    if lambda < 0 || lambda > 1
        error('lambda must lie in [0,1].');
    end

    [x,a,z,b,~] = build_synthetic_marginals();

    yNu = [70;100;130];
    cNu = [0.20;0.60;0.20];

    yAll = sort(unique([x;yNu]));

    c1 = zeros(size(yAll));
    c2 = zeros(size(yAll));

    for i = 1:numel(x)
        idx = find(abs(yAll-x(i)) < 1e-12,1);
        c1(idx) = a(i);
    end

    for j = 1:numel(yNu)
        idx = find(abs(yAll-yNu(j)) < 1e-12,1);
        c2(idx) = cNu(j);
    end

    cAll = (1-lambda)*c1 + lambda*c2;

    keep = cAll > 1e-12;

    y = yAll(keep);
    c = cAll(keep);

    cx12 = check_convex_order(x,a,y,c);
    cx23 = check_convex_order(y,c,z,b);

    if ~cx12.inConvexOrder
        error('mu1 is not below mu2(lambda) in convex order.');
    end

    if ~cx23.inConvexOrder
        error('mu2(lambda) is not below mu3 in convex order.');
    end
end
