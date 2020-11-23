function [X,F,G,H,it,status] = steihaugCG_TR(fun, x0, maxit, tol, tau1, tau2, Delta0)

% Description:
%  Obtains the iterations of a given function with given function
%  handler using trust region method and conjugate gradient for 
%  its subproblem.
%
% Input:
%  fun: function handler
%  x0: initial point
%  maxit: maximum number of iterations
%  tol: tolerance
%  tau1 : treshold for updating delta
%  tau1 : tresholf for updating delta
%  Delta0 : Initial radius
%
% Output:
%  X: Iterations of fhandle up to convergence 
%  F: Final function value 
%  G: Final Gradient
%  H: Final Hessian
%  it: Number of iterations
%  status : indicator of ending state.
%
% Usage:
%  steihaugCG_TR(fun, x0, maxit, tol, tau1, tau2, Delta0)


% Allocate initial point
x(:,1) = x0;
delta(:,1) = Delta0;
[F,G,H] = fun(x0);
status = 0;
model = @(p_x,p_y) F + G'*[p_x;p_y] + 0.5*([p_x p_y]*H*[p_x;p_y]);
i = 1; 

while ( i <= maxit && norm(G) > tol) %Stopping criteria
    
    % Compute p_k
    [p,iter,flag] = steihaugCG(H, G, delta(:,i), tol);
    %syms p_x p_y;
    model = @(p_x,p_y) F + G'*[p_x;p_y] + 0.5*([p_x p_y]*H*[p_x;p_y]);
    
    if mod(i-1, 5) == 0
        figure(i);
        fcontour(model, [-5,5]);
    end

    actual = ( F - fun( x(:,i) + p) );
    pred = ( model(0,0) - model(p(1), p(2)));
    ro =  actual/pred; 

    if ro < tau1
        delta(:,i+1) = delta(:,i) / 2;
        x(:,i+1) = x(:,i);
    else
        x(:,i+1) = x(:,i) + p;
        [F,G,H] = fun(x(:,i+1));
        if ro < tau2;
            delta(:,i+1) = delta(:,i);
        else
            delta(:,i+1) = 2*delta(:,i);
        end
    end

    i=i+1;
end

if i==maxit and norm(G) > tol
	status = 1;
end
X=x;           %Roots of f(x)
it=i-1;        %Number of iterations
figure(it);
fcontour(model, [-5,5]);
draw_trace_rosenbrock(X, fun)


end