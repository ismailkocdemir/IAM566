function [criteria , X, it] = Interior_Point_Mehrotra(A,b,C,X0,Y0,S0,Nmax, tol)
% Author:
%    Ismail Hakki Kocdemir
%
% Description:
%    This function implments Mehrotra's interior point algorithm.
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
%  [criteria , X, it] = Interior_Point_Mehrotra(A,b,C,X0,Y0,S0,Nmax, tol)

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
  eta = 0.995;

  I=eye(n);

  it = 0;
  for k = 1:Nmax 
    mu = (X0'*S0)/n;
    % Finish if tolerance value is reacehed    
    if mu < tol
      break
    end

    F1 = A'*Y0 + S0 - C ;
    F2 = A*X0 - b;
    F3 = X*S*e;
    
    % F(x,\lambda,s)
    F=[ F1;F2;F3];
    
    % Jacobian of F
    DF= [H1 A' I; A H2 H'; S H X];
    
    % Search Direction
    Paff = - DF\F;
    

    dX0_aff = Paff(1:n); % 
    dY0_aff = Paff(n+1:n+m);
    dS0_aff = Paff(n+m+1:end);
        
    alpha_aff_prim = 1;
    alpha_aff_dual = 1;

    neg  = find(S0 + dS0_aff < 0);
    if length(neg)
      alpha_aff_dual =  min(1, min(S0(neg) ./ -dS0_aff(neg)));
    end

    neg  = find(X0 + dX0_aff < 0);
    if length(neg)
      alpha_aff_prim =  min(1, min(X0(neg) ./ -dX0_aff(neg)));
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X0_aff = X0 + alpha_aff_prim*dX0_aff ; 
    S0_aff = S0 + alpha_aff_dual*dS0_aff ; 
    
    mu_aff = (X0_aff'*S0_aff)/n;
    sigma = power((mu_aff/mu),3);

    dX_aff = diag(dX0_aff);
    dS_aff = diag(dS0_aff);
    
    F1= A'*Y0 + S0 - C ;
    F2= A*X0 - b;
    F3= X*S*e + dX_aff*dS_aff*e - sigma*mu*e;
    F=[ F1;F2;F3];

    P = -DF\F;

    dX0 = P(1:n); %
    dY0 = P(n+1:n+m);
    dS0 = P(n+m+1:end);

    alpha_prim_max = 1;
    alpha_dual_max = 1;
    neg  = find(S0 + dS0 < 0);
    if length(neg)
      alpha_dual_max = min(1, min(S0(neg) ./ -dS0(neg)));
    end

    neg  = find(X0 + dX0 < 0);
    if length(neg)
      alpha_prim_max = min(1, min(X0(neg) ./ -dX0(neg)));
    end

    alpha_prim = min(1, eta*alpha_prim_max);
    alpha_dual = min(1, eta*alpha_dual_max);
        
    X0 = X0 + alpha_prim*dX0 ; 
    S0 = S0 + alpha_dual*dS0 ; 
    Y0 = Y0 + alpha_dual*dY0; 

    X = diag(X0) ; 
    S = diag(S0) ;
    it = k;

  end
  % optimal point
  X=diag(X);
  % optimal value
  criteria = C'*X ; 

end

