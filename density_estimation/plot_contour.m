x = linspace(0,3);
y = linspace(0,3);
[X,Y] = meshgrid(x,y);
pdf = zeros(100,100);
for i = 1:100
    for j = 1:100
        pdf(j,i) = fun_pdf([x(:,i), y(:,j)], par);
    end
end
pdf_est = zeros(100,100);
for i = 1:100
    for j = 1:100
        pdf_est(j,i) = fun_pdf([x(:,i), y(:,j)], mean(est_total));
    end
end
subplot(1,2,1)
contour(X,Y,pdf,15)
xlabel('X1')
ylabel('X2')
title('observed pdf for gamma distribution')
subplot(1,2,2)
contour(X,Y,pdf_est,15)
xlabel('X1')
ylabel('X2')
title('estimated pdf for gamma distribution')
