function [x,a,z,b,s0] = synthetic_binomial_pair(N,sigma1,noiseSize)

    if nargin < 2 || isempty(sigma1)
        sigma1 = 10;
    end

    if nargin < 3 || isempty(noiseSize)
        noiseSize = 20;
    end

    s0 = 100;

    %% First marginal

    j = (0:N)';

    x = s0 + (sigma1/sqrt(N)) * (2*j-N);

    a = zeros(N+1,1);

    for q = 0:N
        a(q+1) = nchoosek(N,q)/(2^N);
    end

    %% Final marginal

    rawZ = [x-noiseSize;
            x+noiseSize];

    rawB = [0.5*a;
            0.5*a];

    %% Combine coinciding support values

    [z,~,idx] = unique(rawZ);
    b = accumarray(idx,rawB);

    [z,ord] = sort(z);
    b = b(ord);
end
