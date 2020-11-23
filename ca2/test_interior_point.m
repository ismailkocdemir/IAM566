% test_interior_point.m
%
% Author: Ismail Hskki Kocdemir
%
% Decription:
%  Test script for comparing the algorithms.
%
% The input: None
%
% Ouput: None
%
% Usage:
%  Uncomment 'Feasible input creation' and change x0,x1,s2,s3 to construct feasible initial values and run the script
%   to see the output of Primal-Dual, Mehrotra and linprog's interior point outputs.
%
%  Uncomment 'single infeasable input' to see the outputs for infeasible initial values.


A=[1, 1 ,1 ,0; 2,1,0,1];
C=[-1;1;0;0];
b=[40;60];

% Feasible input creation
x0 = 29;
x1 = 1;
X0 = [x0; x1; 40-(x0+x1); 60-(2*x0+x1)]

s2 = 0.34;
s3 = 0.34;
S0 = [s2+2*s3-1; s2+s3+1; s2; s3];

Y0 = linsolve(A', C-S0)
S0

% single infeasible input 
%X0 = [10;10; 20; 20]
%Y0 = [-10;-10]
%S0 = [1;1;10;10]

intial_value = C'*X0;
intial_value

tol = 1e-6;
Nmax=1e3;


%% CA2 Implementations %%
[criteria_pd , X_pd, it_pd] = Interior_Point_PrimalDual(A,b,C,X0,Y0,S0,Nmax,tol)
[criteria_meh , X_meh, it_meh] = Interior_Point_Mehrotra(A,b,C,X0,Y0,S0,Nmax,tol)


%% linprog implementation %%

lb = zeros(size(X0));

options = optimoptions('linprog','Algorithm','interior-point');
[X_linprog, criteria_linprog, flag, output] = linprog(C,A,b, [], [], lb, [], X0, options);
it_linprog = output.iterations;

criteria_linprog
X_linprog
it_linprog
