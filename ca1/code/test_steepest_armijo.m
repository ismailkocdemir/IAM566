x=[3.1,0.4];

tol=1.0e-3;
maxit=1000;
alpha0 = 1; 
c = 1.0e-2; 
beta = 0.5; 
amax = 100; 

[X,Grad,it] = steepest_descent_armijo(@beale,x,tol,maxit, alpha0, c, beta, amax);
%[X,Grad,it] = steepest_descent_armijo(@rosenbrock,x,tol,maxit, alpha0, c, beta, amax);

x_star = [3.025; 0.474];

pred = X(:,end)
tol
it
diff = norm(X(:,end) - x_star)
normDf = Grad(:,end)

figure(2)
plot(0:it,Grad,'-o')
xlabel('iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('Steepest Descent with Armijo','fontsize',18)