function result = mot_bounds_3marginal_cvx(x,a,y,c,z,b,K)
    x = x(:);
    a = a(:);
    y = y(:);
    c = c(:);
    z = z(:);
    b = b(:);

    m = numel(x);
    n = numel(y);
    ell = numel(z);

    cx12 = check_convex_order(x,a,y,c);
    cx23 = check_convex_order(y,c,z,b);

    if ~cx12.inConvexOrder
        error('mu1 is not below mu2 in convex order.');
    end

    if ~cx23.inConvexOrder
        error('mu2 is not below mu3 in convex order.');
    end

    C = zeros(m*n,ell);

    for i = 1:m
        rows = (i-1)*n + (1:n);
        payoffRow = max(z.' - K*x(i),0);
        C(rows,:) = repmat(payoffRow,n,1);
    end

    %% LOWER BOUND

    cvx_begin quiet

        variable PLower(m*n,ell) nonnegative

        minimize( sum(sum(C .* PLower)) )

        subject to

            % First marginal: sum_{j,k} pi_{ijk} = a_i
            for i = 1:m
                rows = (i-1)*n + (1:n);
                sum(sum(PLower(rows,:))) == a(i);
            end

            % Intermediate marginal: sum_{i,k} pi_{ijk} = c_j
            for j = 1:n
                rows = j:n:(m*n);
                sum(sum(PLower(rows,:))) == c(j);
            end

            % Final marginal: sum_{i,j} pi_{ijk} = b_k
            sum(PLower,1)' == b;

            % First martingale step: E[S2 | S1=x_i] = x_i
            for i = 1:m
                rows = (i-1)*n + (1:n);
                y' * sum(PLower(rows,:),2) == x(i)*a(i);
            end

            % Second martingale step:
            % E[S3 | S1=x_i,S2=y_j] = y_j
            for i = 1:m
                for j = 1:n
                    r = (i-1)*n + j;
                    PLower(r,:)*z == y(j)*sum(PLower(r,:));
                end
            end

    cvx_end

    statusLower = cvx_status;
    lower = cvx_optval;

    if ~status_ok(statusLower)
        error('Three-marginal lower problem ended with CVX status: %s', ...
              statusLower);
    end

    %% UPPER BOUND

    cvx_begin quiet

        variable PUpper(m*n,ell) nonnegative

        maximize( sum(sum(C .* PUpper)) )

        subject to

            for i = 1:m
                rows = (i-1)*n + (1:n);
                sum(sum(PUpper(rows,:))) == a(i);
            end

            for j = 1:n
                rows = j:n:(m*n);
                sum(sum(PUpper(rows,:))) == c(j);
            end

            sum(PUpper,1)' == b;

            for i = 1:m
                rows = (i-1)*n + (1:n);
                y' * sum(PUpper(rows,:),2) == x(i)*a(i);
            end

            for i = 1:m
                for j = 1:n
                    r = (i-1)*n + j;
                    PUpper(r,:)*z == y(j)*sum(PUpper(r,:));
                end
            end

    cvx_end

    statusUpper = cvx_status;
    upper = cvx_optval;

    if ~status_ok(statusUpper)
        error('Three-marginal upper problem ended with CVX status: %s', ...
              statusUpper);
    end

    %% RESULTS

    width = upper-lower;

    if abs(width) < 1e-7
        width = 0;
    end

    piLower = permute(reshape(PLower,[n,m,ell]),[2 1 3]);
    piUpper = permute(reshape(PUpper,[n,m,ell]),[2 1 3]);

   
    pi13Lower = zeros(m,ell);
    pi13Upper = zeros(m,ell);

    for i = 1:m
        rows = (i-1)*n + (1:n);
        pi13Lower(i,:) = sum(PLower(rows,:),1);
        pi13Upper(i,:) = sum(PUpper(rows,:),1);
    end

    result.lower = lower;
    result.upper = upper;
    result.width = width;

    result.PLower = PLower;
    result.PUpper = PUpper;

    result.piLower = piLower;
    result.piUpper = piUpper;

    result.pi13Lower = pi13Lower;
    result.pi13Upper = pi13Upper;

    result.statusLower = statusLower;
    result.statusUpper = statusUpper;

    result.diagnostics.lower = diagnostics(PLower,x,a,y,c,z,b);
    result.diagnostics.upper = diagnostics(PUpper,x,a,y,c,z,b);
end


function d = diagnostics(P,x,a,y,c,z,b)

    m = numel(x);
    n = numel(y);

    %% First marginal residual

    r1 = 0;

    for i = 1:m
        rows = (i-1)*n + (1:n);
        r1 = max(r1,abs(sum(sum(P(rows,:)))-a(i)));
    end

    %% Intermediate marginal residual

    r2 = 0;

    for j = 1:n
        rows = j:n:(m*n);
        r2 = max(r2,abs(sum(sum(P(rows,:)))-c(j)));
    end

    %% Final marginal residual

    r3 = max(abs(sum(P,1)'-b));

    %% First martingale-step residual

    mart12 = 0;

    for i = 1:m
        rows = (i-1)*n + (1:n);
        err = abs(y' * sum(P(rows,:),2) - x(i)*a(i));
        mart12 = max(mart12,err);
    end

    %% Second martingale-step residual

    mart23 = 0;

    for i = 1:m
        for j = 1:n
            r = (i-1)*n + j;
            err = abs(P(r,:)*z - y(j)*sum(P(r,:)));
            mart23 = max(mart23,err);
        end
    end

    %% diagnostics

    d.mu1Residual = r1;
    d.mu2Residual = r2;
    d.mu3Residual = r3;
    d.martingale12Residual = mart12;
    d.martingale23Residual = mart23;
    d.massResidual = abs(sum(P(:))-1);
    d.minProbability = min(P(:));
end


function tf = status_ok(status)

    tf = strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved');

end
