function result = compress_market_marginals_cvx( ...
    fineSupport,fineProb,coarseSupport)
% COMPRESS_MARKET_MARGINALS_CVX
%
% Project the refined fitted marginals onto a smaller common support.
%
% The projection is NOT simple rounding. CVX chooses new probability masses
% so that the coarse marginal call curves remain close to the refined call
% curves while imposing:
%
%   - non-negative probabilities;
%   - unit total mass;
%   - normalized mean equal to one;
%   - mu1 <=cx mu2 <=cx mu3.
%
% The fitting targets are the refined call curves evaluated at every point
% of the refined support.

    fineSupport = fineSupport(:);
    coarseSupport = coarseSupport(:);

    if size(fineProb,1) ~= 3
        error('fineProb must have three rows, one for each maturity.');
    end

    if size(fineProb,2) ~= numel(fineSupport)
        error('fineProb column count must match fineSupport length.');
    end

    nT = 3;
    nFine = numel(fineSupport);
    nCoarse = numel(coarseSupport);

    %% Refined target call curves

    targetCalls = zeros(nT,nFine);

    for q = 1:nFine
        payoffFine = max(fineSupport-fineSupport(q),0);

        for t = 1:nT
            targetCalls(t,q) = fineProb(t,:)*payoffFine;
        end
    end

    %% CVX projection

    cvx_begin quiet

        variable Q(nT,nCoarse) nonnegative
        expression residual(nT*nFine)

        idx = 0;

        for t = 1:nT
            for q = 1:nFine

                idx = idx+1;

                payoffCoarse = ...
                    max(coarseSupport-fineSupport(q),0);

                residual(idx) = ...
                    Q(t,:)*payoffCoarse-targetCalls(t,q);
            end
        end

        minimize( norm(residual,2) )

        subject to

            % Probability mass
            Q*ones(nCoarse,1) == ones(nT,1);

            % Common normalized mean
            Q*coarseSupport == ones(nT,1);

            % Convex order on the common coarse support
            for g = 1:nCoarse

                payoffGrid = ...
                    max(coarseSupport-coarseSupport(g),0);

                Q(1,:)*payoffGrid ...
                    <= Q(2,:)*payoffGrid;

                Q(2,:)*payoffGrid ...
                    <= Q(3,:)*payoffGrid;
            end

    cvx_end

    if ~(strcmp(cvx_status,'Solved') || ...
         strcmp(cvx_status,'Inaccurate/Solved'))

        error('Compression CVX problem ended with status: %s', ...
              cvx_status);
    end

    %% Recalculate compressed call curves on refined strike grid

    compressedCalls = zeros(nT,nFine);

    for q = 1:nFine

        payoffCoarse = ...
            max(coarseSupport-fineSupport(q),0);

        for t = 1:nT
            compressedCalls(t,q) = ...
                Q(t,:)*payoffCoarse;
        end
    end

    errorCalls = compressedCalls-targetCalls;

    %% Diagnostics

    result.status = cvx_status;
    result.objective = cvx_optval;

    result.support = coarseSupport;
    result.probabilities = Q;

    result.targetStrikeGrid = fineSupport;
    result.targetCalls = targetCalls;
    result.compressedCalls = compressedCalls;
    result.callErrors = errorCalls;

    result.rmseCall = sqrt(mean(errorCalls(:).^2));
    result.maxAbsCallError = max(abs(errorCalls(:)));

    result.mass = sum(Q,2);
    result.means = Q*coarseSupport;

    result.maxConvexOrderViolation12 = ...
        max(compressedCalls(1,:)-compressedCalls(2,:));

    result.maxConvexOrderViolation23 = ...
        max(compressedCalls(2,:)-compressedCalls(3,:));
end
