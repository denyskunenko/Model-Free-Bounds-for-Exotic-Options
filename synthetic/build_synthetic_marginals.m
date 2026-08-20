function [x,a,z,b,s0] = build_synthetic_marginals()
    s0 = 100;

    % First marginal
    x = [80;
         90;
         100;
         110;
         120];

    a = [0.10;
         0.20;
         0.40;
         0.20;
         0.10];

    % Final marginal induced by +/-20 mean-zero increment
    z = [60;
         70;
         80;
         90;
         100;
         110;
         120;
         130;
         140];

    b = [0.05;
         0.10;
         0.20;
         0.10;
         0.10;
         0.10;
         0.20;
         0.10;
         0.05];
end
