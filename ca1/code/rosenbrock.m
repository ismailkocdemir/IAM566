function [ f,gradf,hessf ] = rosenbrock(x)

% Description:
% Calculates Rosenbrock function and its gradient at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: Rosenbrock function at point x
%  gradx: gradient of Rosenbrock fuction at point x
%
% Usage:
% [f,gradf] = rosenbrock(x)

f     = 100*(x(1)^2 - x(2))^2 + (x(1)-1)^2;
gradf = [100*(2*(x(1)^2-x(2))*2*x(1)) + 2*(x(1)-1)...
         ;100*(-2*(x(1)^2-x(2)))];
hessf = [-400*x(2) + 1200*x(1)^2 + 2  -400*x(1) ...
		 ; -400*x(1)                 200];

end