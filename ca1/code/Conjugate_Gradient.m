function [X,R,NormR,it] = Conjugate_Gradient(A,b, x0, tol, maxit)

% Author: Ismail Hakki Kocdemir
%
% Description:
%  Obtains the iterations of a given function with given datas by using 
%  BFGS method with updating of the inverse
%
% Input:
%  A: Hilber matrix
%  b: right-side vector
%  x0: initial point
%  tol: tolerance
%  maxit : maximum iterations
%
% Output:
%  X: Iterations of fhandle up to convergence 
%  R: Residuals
%  NormR: Norms of residual
%  it: Number of iterations
%
% Usage:
%  Conjugate_Gradient(x0,tol,Q,c)


% Allocate initial points
x(:,1)=x0;
r(:,1) = b-A*x0;
norm_r(:,1)=norm(r(:,1));
d = r(:,1);
rr = r(:,1)'*r(:,1);
i=1; 

while ( i <= maxit && norm_r(:,i) > tol) %Stopping criteria
    
    % Compute step length
    Ad = A*d;
    alpha = rr/(d'*Ad); 
    % Update the point
    x(:,i+1) = x(:,i) + alpha*d;

    % Compute residual and d.
    r(:,i+1) = r(:,i) - alpha*Ad;  
    norm_r(:,i+1)=norm(r(:,i+1)); 

    beta = (r(:,i+1)'*r(:,i+1))/rr;
    d = r(:,i+1) + beta*d;
    rr = r(:,i+1)'*r(:,i+1);
    
    i=i+1;
end

X=x;           %Roots of f(x)
R=r;         %Residual at roots 
NormR=norm_r; %Norm of Residauls of f(x)
it=i-1;        %Number of iterations

end