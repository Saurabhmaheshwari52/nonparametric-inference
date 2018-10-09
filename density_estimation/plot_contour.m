m = 2;
x = linspace(0,m,100);
y = linspace(0,m,100);
[X,Y] = meshgrid(x,y);
pdf = zeros(100,100);
pdf_est = zeros(100,100);
val = zeros(1,n_dim);
for i = 1:100
    val(1) = x(i);
    for j = 1:100
        rng(i+j)
        val(4) = rand(1)*2;
        val(3) = rand(1)*2;
        val(2) = y(j);
        pdf(j,i) = fun_pdf(par, val, n_comp, n_dim, is_corr);
        pdf_est(j,i) = fun_pdf(mean(est_total)', val, n_comp, n_dim, is_corr);
    end
end
subplot(1,2,1)
contour(X,Y,pdf,45)
xlabel('X1')
ylabel('X2')
title('observed pdf for gm distribution')
subplot(1,2,2)
contour(X,Y,pdf_est,45)
xlabel('X1')
ylabel('X2')
title('estimated pdf for gm distribution')
