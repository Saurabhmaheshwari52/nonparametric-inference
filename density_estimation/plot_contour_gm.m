m = 2;
x = linspace(0,m,100);
y = linspace(0,m,100);
[X,Y] = meshgrid(x,y);
pdf_or = zeros(100,100);
pdf_mar = zeros(100,100);
val = zeros(1,2);
dim1 = 3;
dim2 = 4;
%----------------
par_loc_or = par_loc(:,[dim1, dim2]);
par_scalecorr_or = par_scalecorr([dim1, dim2],[dim1, dim2],:);
gm_or = gmdistribution(par_loc_or, par_scalecorr_or);
%---------------- 
par_est = mean(est_total);
par_loc_mar = par_est([dim1*2-1, dim1*2, dim2*2-1, dim2*2]);
par_loc_mar = reshape(par_loc_mar, 2,2);
par_scale_est = par_est((n_comp*n_dim+1):(end-num_corr));
par_scale_est = reshape(par_scale_est, 1, n_dim, n_comp);
corr_mat_est = zeros(n_dim,n_dim);
corr_mat_est(triu(ones(n_dim),1) == 1) = par_est((end-num_corr+1):end);
corr_mat_est = corr_mat_est+corr_mat_est';
for i = 1:size(par_scale,3)
    prod_mat = prods(par_scale_est(:,:,i), n_dim);
    par_scalecorr_est(:,:,i) = diag(par_scale_est(:,:,i))+corr_mat_est.*prod_mat;
end
par_scalecorr_mar = par_scalecorr_est([dim1, dim2],[dim1, dim2],:);
gm_mar = gmdistribution(par_loc_mar, par_scalecorr_mar);
%----------------
for i = 1:100
    for j = 1:100
        val = [x(i), y(j)];
        pdf_or(j,i) = pdf(gm_or, val);
        pdf_mar(j,i) = pdf(gm_mar, val);
    end
end
subplot(1,2,1)
contour(X,Y,pdf_or,45)
xlabel(strcat('X', num2str(dim1)))
ylabel(strcat('X', num2str(dim2)))
title('observed pdf for gm distribution')
subplot(1,2,2)
contour(X,Y,pdf_mar,45)
xlabel(strcat('X', num2str(dim1)))
ylabel(strcat('X', num2str(dim2)))
title('estimated pdf for gm distribution')
