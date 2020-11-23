function draw_trace(X, f)
n = size(X,2);
Z_res = zeros(1,n);

for i = 1:n
    [Z_res(1,i),~] = f(X(:,i));
end

syms rosenbrock_symsfun(x,y)
rosenbrock_symsfun(x,y) = 100*(x^2 - y)^2 + (x-1)^2;

xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));

figure;
ezsurf(rosenbrock_symsfun,[xmin, xmax, ymin, ymax]);
hold on;
scatter3(X(1, :), X(2, :), Z_res);
hold on
plot3(X(1, :), X(2, :), Z_res, '-x', 'LineWidth', 2);
hold off


end