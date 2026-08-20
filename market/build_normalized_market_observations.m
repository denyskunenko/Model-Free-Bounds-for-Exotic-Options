function [obs,marketInfo] = ...
    build_normalized_market_observations(T,cfg,cfg2)
% BUILD_NORMALIZED_MARKET_OBSERVATIONS
%
% Convert the preliminary OTM quotes into approximate normalized
% European-style call observations.
%
% For K >= F:
%       use the observed OTM call.
%
% For K < F:
%       use the observed OTM put and the approximate parity conversion
%
%       C = P + D(F-K).
%
% Then normalize:
%
%       k = K/F,
%       c_norm(k) = C / (D F).
%
% Under an exact European model this would satisfy
%
%       c_norm(k) = E[(X_T-k)^+],    X_T = S_T/F_T,
%
% with E[X_T] = 1.
%
% SPY options are American-style, so these observations are treated as
% approximate targets rather than exact identities.

    nExp = numel(cfg.expiries);

    obs = table();

    Expiry = NaT(nExp,1);
    DaysToExpiry = zeros(nExp,1);
    TreasuryRate = zeros(nExp,1);
    DiscountFactor = zeros(nExp,1);
    ForwardEstimate = zeros(nExp,1);
    ForwardMAD = zeros(nExp,1);
    RetainedQuotes = zeros(nExp,1);

    for e = 1:nExp

        expiry = cfg.expiries(e);
        dte = cfg.daysToExpiry(e);

        [rate,D] = treasury_discount_factor(dte,cfg);

        chain = build_matched_chain(T,expiry);
        fwd = estimate_forward_spy(chain,D,cfg);

        otm = select_preliminary_otm_quotes( ...
            chain,fwd.forward,cfg);

        otm = otm(otm.Retain,:);

        N = height(otm);

        MaturityIndex = repmat(e,N,1);
        ExpiryObs = repmat(expiry,N,1);
        DaysObs = repmat(dte,N,1);
        DObs = repmat(D,N,1);
        FObs = repmat(fwd.forward,N,1);

        RawStrike = otm.Strike;
        NormalizedStrike = RawStrike/fwd.forward;

        SyntheticCallBid = zeros(N,1);
        SyntheticCallAsk = zeros(N,1);

        for q = 1:N

            if otm.Side(q) == "Put"

                parityTerm = D*(fwd.forward-RawStrike(q));

                SyntheticCallBid(q) = ...
                    otm.Bid(q) + parityTerm;

                SyntheticCallAsk(q) = ...
                    otm.Ask(q) + parityTerm;

            else

                SyntheticCallBid(q) = otm.Bid(q);
                SyntheticCallAsk(q) = otm.Ask(q);
            end
        end

        SyntheticCallMid = ...
            0.5*(SyntheticCallBid+SyntheticCallAsk);

        scaleDollar = D*fwd.forward;

        NormalizedCallBid = ...
            SyntheticCallBid/scaleDollar;

        NormalizedCallAsk = ...
            SyntheticCallAsk/scaleDollar;

        NormalizedCallMid = ...
            SyntheticCallMid/scaleDollar;

        HalfSpreadNormalized = ...
            0.5*(NormalizedCallAsk-NormalizedCallBid);

        FitScale = max( ...
            HalfSpreadNormalized, ...
            cfg2.fitScaleFloor*ones(N,1));

        Side = otm.Side;

        Te = table( ...
            MaturityIndex,ExpiryObs,DaysObs,DObs,FObs, ...
            RawStrike,NormalizedStrike,Side, ...
            SyntheticCallBid,SyntheticCallMid,SyntheticCallAsk, ...
            NormalizedCallBid,NormalizedCallMid,NormalizedCallAsk, ...
            HalfSpreadNormalized,FitScale, ...
            'VariableNames',{ ...
            'MaturityIndex','Expiry','DaysToExpiry', ...
            'DiscountFactor','ForwardEstimate', ...
            'RawStrike','NormalizedStrike','SourceSide', ...
            'SyntheticCallBid','SyntheticCallMid','SyntheticCallAsk', ...
            'NormalizedCallBid','NormalizedCallMid','NormalizedCallAsk', ...
            'HalfSpreadNormalized','FitScale'});

        if isempty(obs)
            obs = Te;
        else
            obs = [obs;Te]; %#ok<AGROW>
        end

        Expiry(e) = expiry;
        DaysToExpiry(e) = dte;
        TreasuryRate(e) = rate;
        DiscountFactor(e) = D;
        ForwardEstimate(e) = fwd.forward;
        ForwardMAD(e) = fwd.medianAbsoluteDeviation;
        RetainedQuotes(e) = N;
    end

    %% Check that every observation lies inside the chosen support range

    s = cfg2.supportNormalized;

    if min(obs.NormalizedStrike) <= min(s) || ...
       max(obs.NormalizedStrike) >= max(s)

        error(['Chosen normalized support does not contain all retained ' ...
               'market strikes strictly inside its endpoints.']);
    end

    marketInfo = table( ...
        (1:nExp)', ...
        Expiry,DaysToExpiry,TreasuryRate,DiscountFactor, ...
        ForwardEstimate,ForwardMAD,RetainedQuotes, ...
        'VariableNames',{ ...
        'MaturityIndex','Expiry','DaysToExpiry','TreasuryRate', ...
        'DiscountFactor','ForwardEstimate','ForwardMAD', ...
        'RetainedQuotes'});
end
