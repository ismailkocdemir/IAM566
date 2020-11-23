x=[3.1,0.4];

tol=1.0e-3;
maxit=1000;

[X,Grad,it] = newton(@beale,x,tol,maxit);
%[X,Grad,it] = newton(@rosenbrock,x,tol,maxit);

figure(2)

x_star = [3.025; 0.474];
pred = X(:,end)
tol
it
diff = norm(X(:,end) - x_star);
diff
normDf = Grad(:,end)


plot([0:it],Grad,'-o');
xlabel('iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('Newton','fontsize',18)