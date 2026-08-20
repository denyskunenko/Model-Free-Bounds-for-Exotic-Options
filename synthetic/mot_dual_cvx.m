function result = mot_dual_cvx(x,a,z,b,K)
    x = x(:);
    a = a(:);
    z = z(:);
    b = b(:);

    m = numel(x);
    n = numel(z);

    C = max(repmat(z.',m,1) - K*repmat(x,1,n),0);
    D = repmat(z.',m,1) - repmat(x,1,n);

    %% LOWER DUAL

    cvx_begin quiet
        variables uLower(m) vLower(n) deltaLower(m)

        maximize( a'*uLower + b'*vLower )

        subject to
            uLower(1) == 0;

            uLower*ones(1,n) ...
            + ones(m,1)*vLower' ...
            + (deltaLower*ones(1,n)).*D ...
            <= C;
    cvx_end

    statusLower = cvx_status;
    lower = cvx_optval;

    if ~(strcmp(statusLower,'Solved') || strcmp(statusLower,'Inaccurate/Solved'))
        error('Lower dual problem ended with CVX status: %s', statusLower);
    end

    %% UPPER DUAL

    cvx_begin quiet
        variables uUpper(m) vUpper(n) deltaUpper(m)

        minimize( a'*uUpper + b'*vUpper )

        subject to
            uUpper(1) == 0;

            uUpper*ones(1,n) ...
            + ones(m,1)*vUpper' ...
            + (deltaUpper*ones(1,n)).*D ...
            >= C;
    cvx_end

    statusUpper = cvx_status;
    upper = cvx_optval;

    if ~(strcmp(statusUpper,'Solved') || strcmp(statusUpper,'Inaccurate/Solved'))
        error('Upper dual problem ended with CVX status: %s', statusUpper);
    end

    %% RESULTS

    result.lower = lower;
    result.upper = upper;
    result.statusLower = statusLower;
    result.statusUpper = statusUpper;

    result.uLower = uLower;
    result.vLower = vLower;
    result.deltaLower = deltaLower;

    result.uUpper = uUpper;
    result.vUpper = vUpper;
    result.deltaUpper = deltaUpper;
end
