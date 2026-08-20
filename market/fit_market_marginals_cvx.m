function result = fit_market_marginals_cvx(obs,support,cfg2)
% FIT_MARKET_MARGINALS_CVX
%
% Jointly fit three discrete probability distributions to the normalized
% market option observations.
%
% Let p_{t,j} be the probability that X_t takes support value x_j.
%
% The model call price at normalized strike k is
%
%       sum_j p_{t,j} max(x_j-k,0).
%
% Constraints:
%
%   1. p_{t,j} >= 0
%   2. sum_j p_{t,j} = 1
%   3. sum_j x_j p_{t,j} = 1
%   4. mu1 <=cx mu2 <=cx mu3
%
% The convex-order inequalities are imposed through call-price
% inequalities at every support point.
%
% The objective is weighted least squares against the approximate
% normalized call observations.

    support = support(:);

    nT = 3;
    nS = numel(support);
    nObs = height(obs);

    tObs = obs.MaturityIndex;
    kObs = obs.NormalizedStrike;
    yObs = obs.NormalizedCallMid;
    fitScale = obs.FitScale;

    %% -----------------------------------------------------------
    % CVX fit
    % ------------------------------------------------------------

    cvx_begin quiet

        variable Prob(nT,nS) nonnegative

        expression residual(nObs)

        for q = 1:nObs

            payoff = max(support-kObs(q),0);

            residual(q) = ...
                (Prob(tObs(q),:)*payoff-yObs(q))/fitScale(q);
        end

        % CVX recommends a norm formulation over an explicit
        % quadratic sum-of-squares objective when possible.
        % Minimising ||residual||_2 has the same minimiser as
        % minimising sum(residual.^2).
        minimize( norm(residual,2) )

        subject to

            % Unit mass at each maturity
            Prob*ones(nS,1) == ones(nT,1);

            % Common normalized mean equal to one
            Prob*support == ones(nT,1);

            % Convex order:
            % call(mu1,k) <= call(mu2,k) <= call(mu3,k)
            for g = 1:nS

                payoffGrid = max(support-support(g),0);

                Prob(1,:)*payoffGrid ...
                    <= Prob(2,:)*payoffGrid;

                Prob(2,:)*payoffGrid ...
                    <= Prob(3,:)*payoffGrid;
            end

    cvx_end

    if ~(strcmp(cvx_status,'Solved') || ...
         strcmp(cvx_status,'Inaccurate/Solved'))

        error('Market marginal fit ended with CVX status: %s', ...
              cvx_status);
    end

    %% -----------------------------------------------------------
    % Recalculate fitted observation values
    % ------------------------------------------------------------

    predictedNorm = zeros(nObs,1);

    for q = 1:nObs

        payoff = max(support-kObs(q),0);

        predictedNorm(q) = ...
            Prob(tObs(q),:)*payoff;
    end

    normalizedError = predictedNorm-yObs;

    predictedDollar = ...
        predictedNorm .* ...
        obs.DiscountFactor .* ...
        obs.ForwardEstimate;

    dollarError = ...
        predictedDollar-obs.SyntheticCallMid;

    insideBidAsk = ...
        predictedNorm >= obs.NormalizedCallBid & ...
        predictedNorm <= obs.NormalizedCallAsk;

    standardizedError = ...
        normalizedError ./ fitScale;

    %% -----------------------------------------------------------
    % Call curves on the support grid
    % ------------------------------------------------------------

    callGrid = zeros(nT,nS);

    for g = 1:nS

        payoffGrid = max(support-support(g),0);

        for t = 1:nT
            callGrid(t,g) = ...
                Prob(t,:)*payoffGrid;
        end
    end

    %% -----------------------------------------------------------
    % Diagnostics
    % ------------------------------------------------------------

    mass = sum(Prob,2);
    means = Prob*support;

    minProbability = min(Prob,[],2);

    maxConvexOrderViolation12 = ...
        max(callGrid(1,:)-callGrid(2,:));

    maxConvexOrderViolation23 = ...
        max(callGrid(2,:)-callGrid(3,:));

    rmseNorm = sqrt(mean(normalizedError.^2));
    rmseDollar = sqrt(mean(dollarError.^2));

    medianAbsDollar = median(abs(dollarError));

    fractionInsideBidAsk = mean(insideBidAsk);

    %% Output

    fittedObs = obs;

    fittedObs.PredictedNormalizedCall = ...
        predictedNorm;

    fittedObs.NormalizedError = ...
        normalizedError;

    fittedObs.StandardizedError = ...
        standardizedError;

    fittedObs.PredictedSyntheticCall = ...
        predictedDollar;

    fittedObs.DollarError = ...
        dollarError;

    fittedObs.InsideBidAsk = ...
        insideBidAsk;

    result.status = cvx_status;
    result.objective = cvx_optval;

    result.support = support;
    result.probabilities = Prob;
    result.callGrid = callGrid;

    result.fittedObservations = fittedObs;

    result.mass = mass;
    result.means = means;
    result.minProbability = minProbability;

    result.maxConvexOrderViolation12 = ...
        maxConvexOrderViolation12;

    result.maxConvexOrderViolation23 = ...
        maxConvexOrderViolation23;

    result.rmseNormalized = rmseNorm;
    result.rmseDollar = rmseDollar;
    result.medianAbsDollarError = medianAbsDollar;
    result.fractionInsideBidAsk = fractionInsideBidAsk;
end
