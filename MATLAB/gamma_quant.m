function [lower, meanie, upper] = gamma_quant(meanie,vars)

% k is shape
% theta is scale
theta = vars ./ meanie;

% shape = 1/theta * mean
% shape = mean/theta
k = meanie ./ theta;

ci = gaminv([0.025,0.975], k, theta);
% x = gaminv(p,a,b) returns the icdf of the gamma distribution with 
% shape parameter a and the scale parameter b, evaluated at the values in p.

lower = ci(1);
upper = ci(2);
end