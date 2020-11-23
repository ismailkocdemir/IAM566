function [f,gradf,hessf] = beale(x)
%
% Author: Ismail Hakki Kocdemir
%
% Description:
%  Calculates Beale function and its gradient and Hessian at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: Rosenbrock function at point x
%  gradf: gradient of Rosenbrock fuction at point x
%  hessf: hessian of Rosenbrock fuction at point x
%
% Usage:
% [f,gradf,hessf] = beale(x)

f = (1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1)+x(1)*x(2)^2)^2 + (2.625 - x(1)+x(1)*x(2)^3)^2;

gradf = [2*(x(2)^2 - 1)*(x(1)*x(2)^2 - x(1) + 9/4) + 2*(x(2)^3 - 1)*(x(1)*x(2)^3 - x(1) + 21/8) + 2*(x(2) - 1)*(x(1)*x(2) - x(1) + 3/2) ...
         ;2*x(1)*(x(1)*x(2) - x(1) + 3/2) + 4*x(1)*x(2)*(x(1)*x(2)^2 - x(1) + 9/4) + 6*x(1)*x(2)^2*(x(1)*x(2)^3 - x(1) + 21/8)];
hessf = [2*(x(2) - 1)^2 + 2*(x(2)^2 - 1)^2 + 2*(x(2)^3 - 1)^2,  6*x(2)^2*(x(1)*x(2)^3 - x(1) + 21/8) - 2*x(1) + 2*x(1)*x(2) + 2*x(1)*(x(2) - 1) + 4*x(2)*(x(1)*x(2)^2 - x(1) + 9/4) + 4*x(1)*x(2)*(x(2)^2 - 1) + 6*x(1)*x(2)^2*(x(2)^3 - 1) + 3 ...
    ; 6*x(2)^2*(x(1)*x(2)^3 - x(1) + 21/8) - 2*x(1) + 2*x(1)*x(2) + 2*x(1)*(x(2) - 1) + 4*x(2)*(x(1)*x(2)^2 - x(1) + 9/4) + 4*x(1)*x(2)*(x(2)^2 - 1) + 6*x(1)*x(2)^2*(x(2)^3 - 1) + 3,  8*x(1)^2*x(2)^2 + 18*x(1)^2*x(2)^4 + 4*x(1)*(x(1)*x(2)^2 - x(1) + 9/4) + 2*x(1)^2 + 12*x(1)*x(2)*(x(1)*x(2)^3 - x(1) + 21/8)];

end