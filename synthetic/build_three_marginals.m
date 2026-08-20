function [x,a,y,c,z,b,s0] = build_three_marginals()
    [x,a,z,b,s0] = build_synthetic_marginals();

    y = [70;
         100;
         130];

    c = [0.20;
         0.60;
         0.20];

    cx12 = check_convex_order(x,a,y,c);
    cx23 = check_convex_order(y,c,z,b);

    if ~cx12.inConvexOrder
        error('mu1 is not below mu2 in convex order.');
    end

    if ~cx23.inConvexOrder
        error('mu2 is not below mu3 in convex order.');
    end
end
