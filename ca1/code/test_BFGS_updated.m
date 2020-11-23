x0=[3.1,0.4];

tol = 1.0e-4;
maxit = 1000;
H = [1 0 ; 0 1];

[X,Df,NormDf,it] = BFGS_updated(@beale,x0,tol,H,maxit);

x_star = [3.025; 0.474];

pred = X(:,end)
tol
it
diff = norm(X(:,end) - x_star)
normdf = normDf(:,end)

figure(2)
plot(0:it,NormDf,'-o');
xlabel('iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('BFGS updated','fontsize',18)
