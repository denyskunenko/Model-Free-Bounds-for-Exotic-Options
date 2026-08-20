clear;
clc;
close all;

% ===============================================================
% REAL MARKET DATA - STAGE 1
%
% This script diagnoses the raw SPY option data before any attempt is
% made to construct risk-neutral marginals.
%
% It does NOT yet feed market data into the MOT solver.
% ===============================================================

cfg = market_data_config();

%% Output directory

if ~exist(cfg.outputDir,'dir')
    mkdir(cfg.outputDir);
end

%% Load raw data

T = load_spy_option_data(cfg);

nExp = numel(cfg.expiries);

Expiry = NaT(nExp,1);
DaysToExpiry = zeros(nExp,1);
TreasuryRate = zeros(nExp,1);
DiscountFactor = zeros(nExp,1);

MatchedStrikes = zeros(nExp,1);
ParityCandidates = zeros(nExp,1);
ForwardEstimate = zeros(nExp,1);
ForwardMAD = zeros(nExp,1);
ForwardMinusSpot = zeros(nExp,1);
RetainedOtmQuotes = zeros(nExp,1);

%% ===============================================================
% Process each maturity
% ===============================================================

for e = 1:nExp

    expiry = cfg.expiries(e);
    dte = cfg.daysToExpiry(e);

    Expiry(e) = expiry;
    DaysToExpiry(e) = dte;

    %% Discount factor

    [rate,D] = treasury_discount_factor(dte,cfg);

    TreasuryRate(e) = rate;
    DiscountFactor(e) = D;

    %% Match calls and puts

    chain = build_matched_chain(T,expiry);

    MatchedStrikes(e) = height(chain);

    %% Approximate synthetic forward

    fwd = estimate_forward_spy(chain,D,cfg);

    ForwardEstimate(e) = fwd.forward;
    ForwardMAD(e) = fwd.medianAbsoluteDeviation;
    ForwardMinusSpot(e) = fwd.forward-cfg.spot;
    ParityCandidates(e) = fwd.numberCandidates;

    %% Preliminary OTM quote selection

    otm = select_preliminary_otm_quotes( ...
        chain,fwd.forward,cfg);

    retained = otm(otm.Retain,:);

    RetainedOtmQuotes(e) = height(retained);

    %% Console output

    fprintf('Expiry %s\n',datestr(expiry,'yyyy-mm-dd'));
    fprintf('  Days to expiry          : %d\n',dte);
    fprintf('  Treasury rate proxy     : %.4f%%\n',100*rate);
    fprintf('  Discount factor         : %.8f\n',D);
    fprintf('  Matched strikes         : %d\n',height(chain));
    fprintf('  Parity candidates       : %d\n',fwd.numberCandidates);
    fprintf('  Forward estimate        : %.6f\n',fwd.forward);
    fprintf('  Forward - spot          : %.6f\n', ...
            fwd.forward-cfg.spot);
    fprintf('  Forward median abs dev  : %.6f\n', ...
            fwd.medianAbsoluteDeviation);
    fprintf('  Preliminary OTM retained: %d\n\n', ...
            height(retained));

    %% Save expiry-level CSVs

    expiryText = datestr(expiry,'yyyy-mm-dd');
    safeExpiry = strrep(expiryText,'-','_');

    writetable(chain, ...
        fullfile(cfg.outputDir, ...
        ['matched_chain_' safeExpiry '.csv']));

    writetable(fwd.candidates, ...
        fullfile(cfg.outputDir, ...
        ['forward_candidates_' safeExpiry '.csv']));

    writetable(otm, ...
        fullfile(cfg.outputDir, ...
        ['preliminary_otm_' safeExpiry '.csv']));

end

%% ===============================================================
% Summary table
% ===============================================================

summary = table( ...
    Expiry, ...
    DaysToExpiry, ...
    TreasuryRate, ...
    DiscountFactor, ...
    MatchedStrikes, ...
    ParityCandidates, ...
    ForwardEstimate, ...
    ForwardMAD, ...
    ForwardMinusSpot, ...
    RetainedOtmQuotes);

disp(summary);

writetable(summary, ...
    fullfile(cfg.outputDir,'market_diagnostics_summary.csv'));

%% Final note

fprintf('\nIMPORTANT\n');
fprintf(['These forward estimates are diagnostics based on near-ATM ' ...
         'call-put differences.\n']);
fprintf(['SPY options are American-style, so they should not be treated ' ...
         'as exact European parity forwards.\n']);
fprintf(['No risk-neutral marginal has been constructed yet. The next ' ...
         'stage will impose no-arbitrage shape conditions\n']);
fprintf(['and convert cleaned option-price curves into discrete ' ...
         'marginals for the MOT problem.\n']);
