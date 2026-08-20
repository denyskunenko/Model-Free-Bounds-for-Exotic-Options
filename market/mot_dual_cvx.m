function result = mot_dual_cvx(x,a,z,b,K)
% MOT_DUAL_CVX
% Explicit finite two-marginal dual for primal-dual validation.
%
% Lower dual (subhedge):
%   maximise a'*u + b'*v
%   subject to
%       u_i + v_j + Delta_i (z_j-x_i) <= (z_j-K x_i)^+.
%
% Upper dual (superhedge):
%   minimise the same static cost subject to the reverse inequality.
%
% Comparing these objectives with MOT_BOUNDS_CVX gives an objective-level
% primal-dual gap. This is different from, and stronger than, merely
% checking primal feasibility residuals.

    x = x(:);
    a = a(:);
    z = z(:);
    b = b(:);

    m = numel(x);
    n = numel(z);

    C = max(repmat(z.',m,1) - K*repmat(x,1,n),0);
    Inc = repmat(z.',m,1) - repmat(x,1,n);

    %% LOWER DUAL

    cvx_begin quiet
        variables uLower(m) vLower(n) deltaLower(m)

        maximize( a'*uLower + b'*vLower )

        subject to
            % Gauge normalisation: adding a constant to u and subtracting
            % it from v leaves both hedge and objective unchanged.
            uLower(1) == 0;

            uLower*ones(1,n) ...
            + ones(m,1)*vLower' ...
            + (deltaLower*ones(1,n)).*Inc ...
            <= C;
    cvx_end

    statusLower = cvx_status;
    lower = cvx_optval;
    slvTolLower = cvx_slvtol;

    if ~status_ok(statusLower)
        error('Two-marginal lower dual ended with CVX status: %s', ...
              statusLower);
    end

    hedgeLower = ...
        uLower*ones(1,n) ...
        + ones(m,1)*vLower' ...
        + (deltaLower*ones(1,n)).*Inc;

    %% UPPER DUAL

    cvx_begin quiet
        variables uUpper(m) vUpper(n) deltaUpper(m)

        minimize( a'*uUpper + b'*vUpper )

        subject to
            uUpper(1) == 0;

            uUpper*ones(1,n) ...
            + ones(m,1)*vUpper' ...
            + (deltaUpper*ones(1,n)).*Inc ...
            >= C;
    cvx_end

    statusUpper = cvx_status;
    upper = cvx_optval;
    slvTolUpper = cvx_slvtol;

    if ~status_ok(statusUpper)
        error('Two-marginal upper dual ended with CVX status: %s', ...
              statusUpper);
    end

    hedgeUpper = ...
        uUpper*ones(1,n) ...
        + ones(m,1)*vUpper' ...
        + (deltaUpper*ones(1,n)).*Inc;

    %% Results

    result.lower = lower;
    result.upper = upper;

    result.statusLower = statusLower;
    result.statusUpper = statusUpper;
    result.slvTolLower = slvTolLower;
    result.slvTolUpper = slvTolUpper;

    result.uLower = uLower;
    result.vLower = vLower;
    result.deltaLower = deltaLower;

    result.uUpper = uUpper;
    result.vUpper = vUpper;
    result.deltaUpper = deltaUpper;

    % Positive values indicate a violation of the dual hedge inequality.
    result.diagnostics.lowerMaxHedgeViolation = ...
        max(max(hedgeLower-C));
    result.diagnostics.upperMaxHedgeViolation = ...
        max(max(C-hedgeUpper));
end


function tf = status_ok(status)
    tf = strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved');
end
