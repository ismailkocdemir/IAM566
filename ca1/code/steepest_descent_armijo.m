function [X,Grad,it] = steepest_descent_armijo(fhandle,x0,tol,maxit,alpha0,c,beta,amax)

% Description:
% Obtains the iterations of a given function with given datas by using 
% steepest descent method with armijo condition
%
% Input:
%  fhandle: function handle
%  x0: initial point
%  tol: tolerance
%  maxit: maximum number of iterations
%  alpha0: initial step size
%  c: parameter
%  beta: parameter
%  amax: maximum number of iterations of Armijo function
%
% Output:
%  X: Iterations of fhandle up to convergence or the number of iterations
%  reaches nmax
%
% Usage:
%  steepest_descent_armijo(fhandle,x0,tol,maxit,alpha0,c,beta,amax)

it = 1;

% Calculate function values of initial point
[~,fgrad] = feval(fhandle,x0);

Grad(:,1)= norm(fgrad);

% Allocate initial point
x(:,1)=x0;

while( it < maxit &&  norm(fgrad) > tol )
  
  % Compute the search direction
  p = -fgrad/norm(fgrad);
  
  % Do the exact line search
  alpha= armijo(fhandle, x(:,it), p, alpha0, c, beta, amax);
  
  % Update the point
  x(:,it+1) = x(:, it) + alpha*p;
  
  % Compute gradient of function at current point for stopping criteria
  [~,fgrad] = feval(fhandle,x(:,it+1));
  
  Grad(:,it+1)=norm(fgrad);
    
  % Update the iteritaions
  it = it+1;
end
  X=x;
  it = it-1
  draw_trace_beale(X, fhandle)
end