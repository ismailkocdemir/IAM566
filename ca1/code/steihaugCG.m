function [p, iter, flag] = steihaugCG(B, g, radius, tol)
% Author:
%  Ismail Hakki Kocdemir
%
% Description:
%  Obtains the update step on a quadratic model with given 
%  hessian and gradient for minimization of the subproblem of trust region.
%
% Input:
%  B: Hessian matrix
%  g: Gradient
%  radius : Update radius
%  tol : tolerance
%
% Output:
%  X: Iterations of fhandle up to convergence 
%  Df: Gradients of f(x)
%  NormDf: Norm of Gradients of f(x)
%  it: Number of iterations
%
% Usage:
%  steihaugCG(B,g,radius,tol)


p = zeros([size(g, 1) 1]);
r(:,1) = g;
norm_r(:,1) = norm(r(:,1));
s = -g;

flag = 0;
i=1; 

%Inverse SR1 Update Algorithm
while ( tol*norm_r(:,1) < norm_r(:,i) ) %Stopping criteria
	sbs = s'*B*s;
	if sbs > 0
		alpha = (r(:,i)'*r(:,i))/sbs;
	else
		p1 = p(1);
		p2 = p(2);
		s1 = s(1);
		s2 = s(2);
		tau = solve('(p1+x*s1)^2 + (p2+x*s2)^2 == radius^2', 'x');
		tau = max(subs(tau));
		p = p + tau.x*s;		
		flag = -1;
		break;
	end

	if norm(p+alpha*s) < radius
		p = p + alpha*s;
	else
		p1 = p(1);
		p2 = p(2);
		s1 = s(1);
		s2 = s(2);
		
		tau = solve('(p1+x*s1)^2 + (p2+x*s2)^2 == radius^2', 'x');
		tau = max(subs(tau));
		p = p + tau*s;
		flag = 1;
		break
	end

	r(:,i+1) = r(:,i) + alpha.*B*s;
	norm_r(:,i+1) = norm(r(:,i+1));

	beta = (r(:,i+1)'*r(:,i+1))/(r(:,i)'*r(:,i));
    s = -r(:,i+1) + beta*s;

	i = i+1;
end

iter = i-1;
		
