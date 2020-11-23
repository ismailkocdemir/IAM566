function [ alpha ] = armijo(fhandle, x, p, alpha0, c, beta, amax)
% Description:
%  Obtains the step size of an iteration satisfied by Armijo condition
%
% Input:
%  fhandle: function handle
%  x: current point
%  alpha0: initial step size
%  c: parameter
%  beta: parameter
%  amax: maximum number of iterations
%
% Output:
%  alpha: The first step size that satisfies Armijo condition
%
% Usage:
%  armijo(fhandle, x, p, alpha0, c, beta, amax)

j=0;

[f,gradx]=feval(fhandle,x);

[fh,~]=feval(fhandle,x+alpha0*p);

while ((fh > f + c*alpha0*gradx'*p)&&(j<amax)) 
    
    alpha0 = alpha0 * beta ;
    
    a  =  x+ alpha0*p ;
    
    [fh,~]=feval(fhandle,a);
    
    j=j+1; 
    
end

alpha=alpha0;

end