tol=1.0e-6;
maxit=1000;

n = 12;
A = hilb(n);
x = zeros([n 1]);
b = ones([n 1]);

[X,R,NormR,it] = Conjugate_Gradient(A,b,x,tol,maxit);

last_res = NormR(:,end)
it
last_x = X(:, end)
condA = cond(A)

figure(2)
plot(0:it,NormR,'-o');
xlabel('iterations','fontsize',18)
ylabel('Norm of Residual','fontsize',18)
title('Conjugate Gradient','fontsize',18)
