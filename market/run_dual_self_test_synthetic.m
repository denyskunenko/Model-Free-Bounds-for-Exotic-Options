clear;
clc;

% ===============================================================
% QUICK SELF-TEST FOR THE NEW EXPLICIT DUAL FUNCTIONS
%
% Uses the dissertation's synthetic K=1 benchmark, where the previously
% validated primal values are approximately:
%   2M lower = 6.33333333, upper = 10.00000000
%   3M lower = 7.55555556, upper =  9.22222222
%
% This test is not a substitute for the market-data validation; it checks
% that the newly added explicit dual formulations match the known primal
% formulations on a small, well-behaved example.
% ===============================================================

cvx_clear;
cvx_solver sdpt3;

x = [80;90;100;110;120];
a = [0.10;0.20;0.40;0.20;0.10];

y = [70;100;130];
c = [0.20;0.60;0.20];

z = [60;70;80;90;100;110;120;130;140];
b = [0.05;0.10;0.20;0.10;0.10;0.10;0.20;0.10;0.05];

K = 1.0;

p2 = mot_bounds_cvx(x,a,z,b,K);
d2 = mot_dual_cvx(x,a,z,b,K);

p3 = mot_bounds_3marginal_cvx(x,a,y,c,z,b,K);
d3 = mot_dual_3marginal_cvx(x,a,y,c,z,b,K);

fprintf('\nSYNTHETIC DUAL SELF-TEST AT K=1\n');
fprintf('2M lower primal / dual : %.10f / %.10f\n',p2.lower,d2.lower);
fprintf('2M upper primal / dual : %.10f / %.10f\n',p2.upper,d2.upper);
fprintf('3M lower primal / dual : %.10f / %.10f\n',p3.lower,d3.lower);
fprintf('3M upper primal / dual : %.10f / %.10f\n',p3.upper,d3.upper);

fprintf('\nAbsolute primal-dual differences\n');
fprintf('2M lower: %.3e\n',abs(p2.lower-d2.lower));
fprintf('2M upper: %.3e\n',abs(p2.upper-d2.upper));
fprintf('3M lower: %.3e\n',abs(p3.lower-d3.lower));
fprintf('3M upper: %.3e\n',abs(p3.upper-d3.upper));

fprintf('\nExpected approximate primal benchmark values\n');
fprintf('2M: [6.33333333, 10.00000000]\n');
fprintf('3M: [7.55555556,  9.22222222]\n');
