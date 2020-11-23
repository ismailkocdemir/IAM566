function [X,Df,NormDf,it] = BFGS_updated(fhandle, x0, tol, H0, maxit)
% Author:
%  Ismail Hakki Kocdemir
%
% Description:
%  Obtains the iterations of a given function with given initial point by using 
%  BFGS method with updating of the inverse
%
% Input:
%  fhandle: used to evaluate a function.
%  x0: initial point
%  tol: tolerance
%  H0 : Hessian matrix
%  c : vector
%
% Output:
%  X: Iterations of fhandle up to convergence 
%  Df: Gradients of f(x)
%  NormDf: Norm of Gradients of f(x)
%  it: Number of iterations
%
% Usage:
%  BFGS_updated(fhandle, x0, tol, H0, maxit)


% Allocate initial point
x(:,1)=x0;

H=H0;

[~,df(:,1),Q] = feval(fhandle,x0);

normdf(:,1)=norm(df(:,1));

i=1; 

%Inverse SR1 Update Algorithm
while ( i < maxit && normdf(:,i)> tol) %Stopping criteria

    % Compute the search direction
    p=-H*df(:, i); 
    
    % Compute step length
    alpha = -(df(:, i)'*p)/(p'*Q*p); 
    
    % Update the point
    x(:,i+1) = x(:,i) + alpha.*p;

    % Compute gradient and Hessian
    [~,df(:,i+1),Q] = feval(fhandle,x(:,i+1));  
    
    % Compute norm of gradient
    normdf(:,i+1)=norm(df(:,i+1)); 
    
    s=x(:,i+1)-x(:,i);
    
    y=df(:,i+1)-df(:,i);
    
    %Inverse update formula
    H = H + ((s-H*y)*(s-H*y)')/((s-H*y)'*y); 
    
    i=i+1;
end

X=x;           %Roots of f(x)
Df=df;         %Gradients of f(x) at roots 
NormDf=normdf; %Norm of Gradients of f(x)
it=i-1;        %Number of iterations
draw_trace_beale(X, fhandle)

end