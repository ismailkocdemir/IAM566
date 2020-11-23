function [X,Grad,it] = newton(fhandle,x0,tol,maxit)

% Description:
% Obtains the iterations of a given function with given datas by using 
% steepest descent method with exact line search
%
% Input:
%  fhandle: function handle
%  x0: initial point
%  tol: tolerance
%  maxit: maximum number of iterations
%
% Output:
%  X: Iterations of fhandle up to convergence or the number of iterations
%  reaches nmax
%
% Usage:
%  newton(fhandle,x0,tol,maxit)

it = 1;

% Calculate function values of initial point
[~,fgrad,hessian] = feval(fhandle,x0);

Grad(:,1)= norm(fgrad);

% Allocate initial point
x(:,1)=x0;

while( it < maxit && norm(fgrad) > tol )
  
  % Compute the search direction
  p = -1*inv(hessian)*fgrad;
  
  % Do the exact line search
  alpha = 1;
  
  % Update the point
  x(:,it+1)= x(:, it) + alpha*p;
  
  % Compute gradient of function at current point for stopping criteria
  [~,fgrad,hessian] = feval(fhandle,x(:,it+1));
  
  Grad(:,it+1)=norm(fgrad);
  
  % Update the iteritaions
  it = it+1;
end
  X=x;
  it = it-1
  draw_trace_beale(X, fhandle)
end