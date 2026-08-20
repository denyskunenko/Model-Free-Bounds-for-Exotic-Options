function out = estimate_forward_spy(chain,D,cfg)
% ESTIMATE_FORWARD_SPY
% Estimate an approximate forward from near-ATM call-put differences.
%
% For European options,
%
%       C - P = D(F-K)
%
% and therefore
%
%       F = K + (C-P)/D.
%
% SPY options are American-style, so this relation is used only as a
% practical near-ATM diagnostic. Deep-ITM quotes are excluded.
%
% The final estimate is the median of the strike-level implied forwards.

    nearATM = abs(chain.Strike/cfg.spot - 1) ...
              <= cfg.parityMoneynessWindow;

    positiveBid = ...
        chain.CallBid > 0 & ...
        chain.PutBid  > 0;

    spreadOK = ...
        chain.CallRelSpread <= cfg.maxParityRelativeSpread & ...
        chain.PutRelSpread  <= cfg.maxParityRelativeSpread;

    valid = ...
        chain.ValidCall & ...
        chain.ValidPut & ...
        nearATM & ...
        positiveBid & ...
        spreadOK;

    candidates = chain(valid,:);

    if height(candidates) < cfg.minParityCandidates
        error(['Only %d near-ATM parity candidates survived the filters. ' ...
               'The forward estimate is not sufficiently supported.'], ...
               height(candidates));
    end

    candidates.ParityForward = ...
        candidates.Strike ...
        + (candidates.CallMid-candidates.PutMid)/D;

    Fhat = median(candidates.ParityForward,'omitnan');

    madForward = median( ...
        abs(candidates.ParityForward-Fhat), ...
        'omitnan');

    candidates.ForwardDeviation = ...
        candidates.ParityForward-Fhat;

    out.forward = Fhat;
    out.medianAbsoluteDeviation = madForward;
    out.candidates = candidates;
    out.numberCandidates = height(candidates);
    out.minimumForward = min(candidates.ParityForward);
    out.maximumForward = max(candidates.ParityForward);
end
