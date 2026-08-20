function out = select_preliminary_otm_quotes(chain,Fhat,cfg)
% SELECT_PRELIMINARY_OTM_QUOTES
% Choose the OTM side of the option surface as a diagnostic.
%
% If K < Fhat, use the put.
% If K >= Fhat, use the call.
%
% Quotes are flagged Retain=true only when the selected OTM quote has a
% valid bid/ask pair, a positive bid (if requested), and a relative spread
% below the configured threshold.
%
% This is NOT yet the final no-arbitrage cleaning step.

    N = height(chain);

    Side = strings(N,1);
    Bid = nan(N,1);
    Ask = nan(N,1);
    Mid = nan(N,1);
    Spread = nan(N,1);
    RelSpread = nan(N,1);
    Vol = nan(N,1);

    for q = 1:N

        if chain.Strike(q) < Fhat

            Side(q) = "Put";
            Bid(q) = chain.PutBid(q);
            Ask(q) = chain.PutAsk(q);
            Mid(q) = chain.PutMid(q);
            Spread(q) = chain.PutSpread(q);
            RelSpread(q) = chain.PutRelSpread(q);
            Vol(q) = chain.PutVol(q);

        else

            Side(q) = "Call";
            Bid(q) = chain.CallBid(q);
            Ask(q) = chain.CallAsk(q);
            Mid(q) = chain.CallMid(q);
            Spread(q) = chain.CallSpread(q);
            RelSpread(q) = chain.CallRelSpread(q);
            Vol(q) = chain.CallVol(q);
        end
    end

    Retain = ...
        isfinite(Bid) & ...
        isfinite(Ask) & ...
        isfinite(Mid) & ...
        Bid >= 0 & ...
        Ask >= Bid & ...
        Mid > 0 & ...
        RelSpread <= cfg.maxOtmRelativeSpread;

    if cfg.requirePositiveBid
        Retain = Retain & Bid > 0;
    end

    Moneyness = chain.Strike/Fhat;

    out = table( ...
        chain.Strike, ...
        Moneyness, ...
        Side, ...
        Bid,Ask,Mid,Spread,RelSpread,Vol,Retain, ...
        'VariableNames',{ ...
        'Strike','StrikeOverForward','Side', ...
        'Bid','Ask','Mid','Spread','RelativeSpread','Vol','Retain'});
end
