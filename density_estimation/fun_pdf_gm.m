% For bimodal distributions
function [fun_pdf_val] = fun_pdf_gm(par,x,n_comp, n_dim, is_corr)
par_l = reshape(par(1:n_comp*n_dim), n_comp, n_dim);
if ~is_corr
    par_s = reshape(par(n_comp*n_dim+1:end)', 1, n_dim, n_comp);
    par_c = [];
    par_scalecorr = par_s;
else
    n_corr = n_dim*(n_dim-1)/2;
    par_s = reshape(par(n_comp*n_dim+1:(end-n_corr))', 1, n_dim, n_comp);
    par_c = par((end-n_corr+1):size(par,1));
    corr_mat = zeros(n_dim,n_dim);
    corr_mat(triu(ones(n_dim),1) == 1) = par_c;
    corr_mat = corr_mat+corr_mat';
    rng(1);
    for i = 1:size(par_s,3)
        prod_mat = prods(par_s(:,:,i), n_dim);
        par_scalecorr(:,:,i) = diag(par_s(:,:,i))+corr_mat.*prod_mat;
    end%     [~,b] = chol(par_scalecorr(:,:,i))
% 	par_scalecorr(:,:,i)
end
%par_l1 = par(1:n_c);
%par_l2 = par((n_c+1):n_c*n_d);
%par_s1 = reshape(par(n_c*n_d+1:n_c*n_d+n_c), 1, 1, n_c);
%par_s2 = reshape(par(n_c*n_d+n_c+1:end), 1, 1, n_c);
%gm1 = gmdistribution(par_l1, par_s1);
%gm2 = gmdistribution(par_l2, par_s2);
%fun_pdf_val = pdf(gm1, x(:,1)).*pdf(gm2, x(:,2));
%par_l
%par_scalecorr
gm1 = gmdistribution(par_l, par_scalecorr);
fun_pdf_val = pdf(gm1, x);
end
