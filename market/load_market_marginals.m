function marketMarginals = load_market_marginals()
% LOAD_MARKET_MARGINALS Locate and load the refined Stage-2 MAT file.

    candidates = { ...
        fullfile('market_marginal_output','market_marginals.mat'), ...
        'market_marginals.mat'};

    found = "";

    for q = 1:numel(candidates)
        if isfile(candidates{q})
            found = string(candidates{q});
            break;
        end
    end

    if strlength(found) == 0
        error(['Cannot find market_marginals.mat. Run the refined Stage-2 ' ...
               'script first, or copy market_marginals.mat into the current ' ...
               'MATLAB folder.']);
    end

    S = load(found);

    if ~isfield(S,'marketMarginals')
        error('%s does not contain the variable marketMarginals.',found);
    end

    marketMarginals = S.marketMarginals;

    required = { ...
        'spot', ...
        'supportNormalized', ...
        'probabilities', ...
        'forwardEstimates', ...
        'discountFactors'};

    for q = 1:numel(required)
        if ~isfield(marketMarginals,required{q})
            error('marketMarginals is missing field %s.',required{q});
        end
    end

    fprintf('\nLoaded refined market marginals from:\n  %s\n\n',found);
end
