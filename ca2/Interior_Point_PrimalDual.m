function [criteria , X, it] = Interior_Point_PrimalDual(A,b,C,X0,Y0,S0,Nmax, tol)
% Author:
%    Ismail Hakki Kocdemir
%
% Description:
%    This function implments the primal dual interior point algorithm.
%    Consider min  c'*X    Primal problem
%             s.t. Ax <= b 
%                   x >= 0
%            
%             max  b'*Y     Dual problem
%                  A'*Y <=c iff A'*Y + S  = c
%                                        S >= 0
%
% The input:
%   A: This is the Ax=b matrix. 
%   b: Vector. This is the right hand side of Ax=b.
%   c: Vector. This is from minimize  J(x) = c'x. 
%   X0: initial vector
%   S0: initial slack vector
%   Y0: initial vector
%   Nmax: Maximum number of iterations
%   tol: Tolerance value to stop the algorithm.
%
% Ouput:
%   criteria = c'*X
%   X:optimal value
%   it: number of iterations
%
% Usage:
%  [criteria , X, it] = Interior_Point_PrimalDual(A,b,C,X0,Y0,S0,Nmax, tol)

  % Diagonal Matrices
  X = diag(X0); 
  S = diag(S0); 

  % size of matrices
  n = size(X,1);
  m = size(Y0,1);
  e = ones(n,1);

  % Zeros matrix in DF
  H  = zeros(n,m);
  H1 = zeros(n,n);
  H2 = zeros(m,m);
  mu = 1;


  I=eye(n);
  sigma = 0.3 ;

  it = 0;
  for k = 1:Nmax 
    mu = (X0'*S0)/n;

    if mu < tol
      break
    end

    F1= zeros(n,1); %A'*Y0 + S0 - C ;
    F2= zeros(m,1); %A*X0 - b;
    F3= X*S*e - sigma*mu*e;
    
    % F(x,\lambda,s)
    F=[ F1;F2;F3];
    
    % Jacobian of F
    DF= [H1 A' I; A H2 H'; S H X];
    
    %Rc = A'*Y0 + S0 - C ;
    %Rb = A*X0 - b;
    %R = [Rc; Rb];

    % Search Direction
    P = - DF\F;
    
    dX = P(1:n); % 
    dY = P(n+1:n+m);
    dS = P(n+m+1:end);
        
    alpha_prim = 1;
    alpha_dual = 1;
    neg  = find(S0 + dS < 0);
    if length(neg)
      alpha_dual = min(1, min(S0(neg) ./ -dS(neg)));
    end

    neg  = find(X0 + dX < 0);
    if length(neg)
      alpha_prim = min(1, min(X0(neg) ./ -dX(neg)));
    end

    % Update all variables
    X0 = X0 + alpha_prim*dX ; 
    S0 = S0 + alpha_dual*dS ; 
    Y0 = Y0 + alpha_dual*dY; 

    X = diag(X0) ; 
    S = diag(S0) ;
    it = k;
  end
  % optimal point
  X=diag(X);
  % optimal value
  criteria = C'*X ; 

end

