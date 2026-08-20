function T = load_spy_option_data(cfg)
% LOAD_SPY_OPTION_DATA Load and validate the raw SPY option-chain CSV.

    if ~isfile(cfg.optionCsv)
        error(['Cannot find %s. Put the CSV in the same MATLAB folder ' ...
               'as the market-data code.'],cfg.optionCsv);
    end

    T = readtable(cfg.optionCsv,'TextType','string');

    required = [ ...
        "date","act_symbol","expiration","strike", ...
        "call_put","bid","ask","vol"];

    names = string(T.Properties.VariableNames);

    if ~all(ismember(required,names))
        missing = required(~ismember(required,names));
        error('CSV is missing required columns: %s', ...
              strjoin(missing,", "));
    end

    %% Convert dates

    if ~isdatetime(T.date)
        T.date = datetime(string(T.date),'InputFormat','yyyy-MM-dd');
    end

    if ~isdatetime(T.expiration)
        T.expiration = datetime(string(T.expiration), ...
                                'InputFormat','yyyy-MM-dd');
    end

    %% Convert text fields

    T.act_symbol = string(T.act_symbol);
    T.call_put = string(T.call_put);

    %% Convert numeric fields if necessary

    numericNames = ["strike","bid","ask","vol"];

    for q = 1:numel(numericNames)

        name = numericNames(q);

        if ~isnumeric(T.(name))
            T.(name) = str2double(string(T.(name)));
        end
    end

    %% Restrict to the intended experiment

    keep = ...
        T.act_symbol == cfg.symbol & ...
        T.date == cfg.quoteDate & ...
        ismember(T.expiration,cfg.expiries);

    T = T(keep,:);

    if isempty(T)
        error('No observations remain after applying symbol/date/expiry filters.');
    end

    %% Basic quote checks

    if any(~isfinite(T.strike)) || any(T.strike <= 0)
        error('Invalid strike values found in the raw CSV.');
    end

    if any(~isfinite(T.bid)) || any(~isfinite(T.ask))
        error('Missing or non-finite bid/ask values found in the raw CSV.');
    end

    if any(T.bid < 0)
        error('Negative bid values found in the raw CSV.');
    end

    if any(T.ask < T.bid)
        error('At least one option quote has ask < bid.');
    end

    %% Report

    fprintf('\nRaw option data loaded\n');
    fprintf('Rows retained       : %d\n',height(T));
    fprintf('Symbol              : %s\n',cfg.symbol);
    fprintf('Quote date          : %s\n',datestr(cfg.quoteDate,'yyyy-mm-dd'));
    fprintf('Selected expiries   : %d\n\n',numel(cfg.expiries));
end
