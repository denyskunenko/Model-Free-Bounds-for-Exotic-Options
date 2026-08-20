function [rate,D] = treasury_discount_factor(daysToExpiry,cfg)
% TREASURY_DISCOUNT_FACTOR
% Approximate a short zero-rate input from the observed Treasury curve.
%
% For maturities shorter than the first tenor, the 1-month yield is used.
% Between observed tenors, linear interpolation in maturity is used.
%
% The resulting yield is then treated as a continuously compounded
% short-rate proxy:
%
%       D(T) = exp(-r T),   T = days/365.
%
% This is an empirical approximation, not a bootstrapped zero curve.

    tenor = cfg.treasuryTenorDays(:);
    rates = cfg.treasuryRates(:);

    if daysToExpiry <= tenor(1)

        rate = rates(1);

    elseif daysToExpiry >= tenor(end)

        rate = interp1(tenor,rates,daysToExpiry,'linear','extrap');

    else

        rate = interp1(tenor,rates,daysToExpiry,'linear');
    end

    T = daysToExpiry/365;

    D = exp(-rate*T);
end
