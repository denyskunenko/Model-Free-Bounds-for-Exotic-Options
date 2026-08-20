function out = check_convex_order(x,a,z,b,tol)
% CHECK_CONVEX_ORDER Check common mean and discrete convex order.
%
% For finite supports with a common mean, convex order can be checked
% through call-price inequalities at the union of the support points.

    if nargin < 5
        tol = 1e-10;
    end

    x = x(:);
    a = a(:);
    z = z(:);
    b = b(:);

    if numel(x) ~= numel(a) || numel(z) ~= numel(b)
        error('Support and probability vectors must have matching lengths.');
    end

    if any(a < -tol) || any(b < -tol)
        error('Probabilities must be non-negative.');
    end

    if abs(sum(a)-1) > 1e-8 || abs(sum(b)-1) > 1e-8
        error('Each marginal must sum to one.');
    end

    mean1 = a' * x;
    mean3 = b' * z;
    meanDifference = mean1 - mean3;

    strikeGrid = sort(unique([x; z]));

    call1 = zeros(size(strikeGrid));
    call3 = zeros(size(strikeGrid));

    for q = 1:numel(strikeGrid)
        k = strikeGrid(q);
        call1(q) = sum(a .* max(x-k,0));
        call3(q) = sum(b .* max(z-k,0));
    end

    violations = call1-call3;
    maxViolation = max(violations);

    out.mean1 = mean1;
    out.mean3 = mean3;
    out.meanDifference = meanDifference;
    out.strikeGrid = strikeGrid;
    out.call1 = call1;
    out.call3 = call3;
    out.maxViolation = maxViolation;
    out.commonMean = abs(meanDifference) <= 1e-8;
    out.inConvexOrder = out.commonMean && maxViolation <= 1e-8;
end
