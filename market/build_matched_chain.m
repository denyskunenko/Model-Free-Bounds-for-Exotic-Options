function chain = build_matched_chain(T,expiry)
% BUILD_MATCHED_CHAIN Create one row per strike containing both call and put.
%
% Output columns include:
%   Strike
%   CallBid, CallAsk, CallMid, CallSpread, CallRelSpread, CallVol
%   PutBid,  PutAsk,  PutMid,  PutSpread,  PutRelSpread,  PutVol

    G = T(T.expiration == expiry,:);

    calls = G(G.call_put == "Call",:);
    puts  = G(G.call_put == "Put",:);

    if isempty(calls) || isempty(puts)
        error('Missing call or put observations for expiry %s.', ...
              datestr(expiry,'yyyy-mm-dd'));
    end

    %% Avoid silently combining duplicate rows

    if numel(unique(calls.strike)) ~= height(calls)
        error('Duplicate call strikes found for expiry %s.', ...
              datestr(expiry,'yyyy-mm-dd'));
    end

    if numel(unique(puts.strike)) ~= height(puts)
        error('Duplicate put strikes found for expiry %s.', ...
              datestr(expiry,'yyyy-mm-dd'));
    end

    %% Match by common strike

    [K,ic,ip] = intersect(calls.strike,puts.strike,'sorted');

    if isempty(K)
        error('No matched call-put strikes for expiry %s.', ...
              datestr(expiry,'yyyy-mm-dd'));
    end

    chain = table( ...
        K, ...
        calls.bid(ic), calls.ask(ic), calls.vol(ic), ...
        puts.bid(ip),  puts.ask(ip),  puts.vol(ip), ...
        'VariableNames',{ ...
        'Strike', ...
        'CallBid','CallAsk','CallVol', ...
        'PutBid','PutAsk','PutVol'});

    %% Midpoints and spreads

    chain.CallMid = 0.5*(chain.CallBid + chain.CallAsk);
    chain.PutMid  = 0.5*(chain.PutBid  + chain.PutAsk);

    chain.CallSpread = chain.CallAsk-chain.CallBid;
    chain.PutSpread  = chain.PutAsk-chain.PutBid;

    chain.CallRelSpread = ...
        chain.CallSpread ./ max(chain.CallMid,1e-12);

    chain.PutRelSpread = ...
        chain.PutSpread ./ max(chain.PutMid,1e-12);

    %% Basic validity flags

    chain.ValidCall = ...
        chain.CallBid >= 0 & ...
        chain.CallAsk >= chain.CallBid & ...
        chain.CallMid > 0;

    chain.ValidPut = ...
        chain.PutBid >= 0 & ...
        chain.PutAsk >= chain.PutBid & ...
        chain.PutMid > 0;
end
