function result = mot_bounds_cvx(x,a,z,b,K)

    x = x(:);
    a = a(:);
    z = z(:);
    b = b(:);

    m = numel(x);
    n = numel(z);
    cx = check_convex_order(x,a,z,b);

    if ~cx.commonMean
        error('Marginals do not have the same mean.');
    end

    if ~cx.inConvexOrder
        error('Marginals fail convex order.');
    end

    C = max(repmat(z.',m,1) - K*repmat(x,1,n),0);

    %% LOWER BOUND

    cvx_begin quiet
        variable piLower(m,n) nonnegative

        minimize( sum(sum(C .* piLower)) )

        subject to
            sum(piLower,2) == a;
            sum(piLower,1)' == b;
            piLower*z == a.*x;
    cvx_end

    statusLower = cvx_status;
    lower = cvx_optval;

    if ~status_ok(statusLower)
        error('Lower problem ended with CVX status: %s', statusLower);
    end

    %% UPPER BOUND

    cvx_begin quiet
        variable piUpper(m,n) nonnegative

        maximize( sum(sum(C .* piUpper)) )

        subject to
            sum(piUpper,2) == a;
            sum(piUpper,1)' == b;
            piUpper*z == a.*x;
    cvx_end

    statusUpper = cvx_status;
    upper = cvx_optval;

    if ~status_ok(statusUpper)
        error('Upper problem ended with CVX status: %s', statusUpper);
    end

    %% RESULTS

    width = upper-lower;
    if abs(width) < 1e-7
        width = 0;
    end

    result.lower = lower;
    result.upper = upper;
    result.width = width;

    result.piLower = piLower;
    result.piUpper = piUpper;
    result.payoffMatrix = C;

    result.statusLower = statusLower;
    result.statusUpper = statusUpper;

    result.mean1 = cx.mean1;
    result.mean3 = cx.mean3;
    result.convexOrderViolation = cx.maxViolation;

    result.diagnostics.lower = diagnostics(piLower,x,a,z,b);
    result.diagnostics.upper = diagnostics(piUpper,x,a,z,b);
end


function d = diagnostics(pi,x,a,z,b)
    d.rowResidual = max(abs(sum(pi,2)-a));
    d.colResidual = max(abs(sum(pi,1)'-b));
    d.martingaleResidual = max(abs(pi*z-a.*x));
    d.massResidual = abs(sum(pi(:))-1);
    d.minProbability = min(pi(:));
end


function tf = status_ok(status)
    tf = strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved');
end
