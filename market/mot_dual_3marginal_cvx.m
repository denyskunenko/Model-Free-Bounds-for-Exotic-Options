function result = mot_dual_3marginal_cvx(x,a,y,c,z,b,K)
% MOT_DUAL_3MARGINAL_CVX
% Explicit finite three-marginal dual for primal-dual validation.
%
% The semi-static hedge has the form
%
%   u_i + v_j + w_k
%   + Delta12_i (y_j-x_i)
%   + Delta23_ij (z_k-y_j).
%
% The lower dual maximises the static cost subject to this hedge lying
% below the payoff at every support triple. The upper dual minimises the
% static cost subject to the reverse inequality.
%
% This is the discrete counterpart of the three-period semi-static dual
% discussed in the dissertation. Its objective can be compared directly
% with MOT_BOUNDS_3MARGINAL_CVX to assess objective optimality.

    x = x(:);
    a = a(:);
    y = y(:);
    c = c(:);
    z = z(:);
    b = b(:);

    m = numel(x);
    n = numel(y);
    ell = numel(z);

    %% LOWER DUAL

    cvx_begin quiet
        variables uLower(m) vLower(n) wLower(ell)
        variables delta12Lower(m) delta23Lower(m,n)

        maximize( a'*uLower + c'*vLower + b'*wLower )

        subject to
            % Gauge normalisations. They remove redundant directions but
            % do not change the attainable hedge values or objective.
            uLower(1) == 0;
            vLower(1) == 0;
            delta12Lower(1) == 0;
            delta23Lower(1,1) == 0;

            for i = 1:m
                payoff = max(z - K*x(i),0);

                for j = 1:n
                    uLower(i) + vLower(j) + wLower ...
                        + delta12Lower(i)*(y(j)-x(i)) ...
                        + delta23Lower(i,j)*(z-y(j)) ...
                        <= payoff;
                end
            end
    cvx_end

    statusLower = cvx_status;
    lower = cvx_optval;
    slvTolLower = cvx_slvtol;

    if ~status_ok(statusLower)
        error('Three-marginal lower dual ended with CVX status: %s', ...
              statusLower);
    end

    lowerViolation = dual_violation_lower( ...
        x,y,z,K,uLower,vLower,wLower,delta12Lower,delta23Lower);

    %% UPPER DUAL

    cvx_begin quiet
        variables uUpper(m) vUpper(n) wUpper(ell)
        variables delta12Upper(m) delta23Upper(m,n)

        minimize( a'*uUpper + c'*vUpper + b'*wUpper )

        subject to
            uUpper(1) == 0;
            vUpper(1) == 0;
            delta12Upper(1) == 0;
            delta23Upper(1,1) == 0;

            for i = 1:m
                payoff = max(z - K*x(i),0);

                for j = 1:n
                    uUpper(i) + vUpper(j) + wUpper ...
                        + delta12Upper(i)*(y(j)-x(i)) ...
                        + delta23Upper(i,j)*(z-y(j)) ...
                        >= payoff;
                end
            end
    cvx_end

    statusUpper = cvx_status;
    upper = cvx_optval;
    slvTolUpper = cvx_slvtol;

    if ~status_ok(statusUpper)
        error('Three-marginal upper dual ended with CVX status: %s', ...
              statusUpper);
    end

    upperViolation = dual_violation_upper( ...
        x,y,z,K,uUpper,vUpper,wUpper,delta12Upper,delta23Upper);

    %% Results

    result.lower = lower;
    result.upper = upper;

    result.statusLower = statusLower;
    result.statusUpper = statusUpper;
    result.slvTolLower = slvTolLower;
    result.slvTolUpper = slvTolUpper;

    result.uLower = uLower;
    result.vLower = vLower;
    result.wLower = wLower;
    result.delta12Lower = delta12Lower;
    result.delta23Lower = delta23Lower;

    result.uUpper = uUpper;
    result.vUpper = vUpper;
    result.wUpper = wUpper;
    result.delta12Upper = delta12Upper;
    result.delta23Upper = delta23Upper;

    result.diagnostics.lowerMaxHedgeViolation = lowerViolation;
    result.diagnostics.upperMaxHedgeViolation = upperViolation;
end


function v = dual_violation_lower(x,y,z,K,u,v2,w,d12,d23)
    m = numel(x);
    n = numel(y);
    v = -Inf;

    for i = 1:m
        payoff = max(z-K*x(i),0);
        for j = 1:n
            hedge = u(i)+v2(j)+w ...
                + d12(i)*(y(j)-x(i)) ...
                + d23(i,j)*(z-y(j));
            v = max(v,max(hedge-payoff));
        end
    end
end


function v = dual_violation_upper(x,y,z,K,u,v2,w,d12,d23)
    m = numel(x);
    n = numel(y);
    v = -Inf;

    for i = 1:m
        payoff = max(z-K*x(i),0);
        for j = 1:n
            hedge = u(i)+v2(j)+w ...
                + d12(i)*(y(j)-x(i)) ...
                + d23(i,j)*(z-y(j));
            v = max(v,max(payoff-hedge));
        end
    end
end


function tf = status_ok(status)
    tf = strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved');
end
