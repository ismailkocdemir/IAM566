x0 = [0;-1];
tol=1.0e-3;
tau1=1/4;
tau2=3/4;
Delta0 = 1;
maxit=1000;

[X,F,G,H,it,status] = steihaugCG_TR(@rosenbrock, x0, maxit, tol, tau1, tau2, Delta0);
status
X
it