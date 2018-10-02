% For bimodal distributions
function [fun_pdf_val] = fun_pdf_gm(par,x,n_c, n_d, is_corr)
par_l = reshape(par(1:n_c*n_d), n_c, n_d);
if ~is_corr
    par_s = reshape(par(n_c*n_d+1:end)', 1, n_d, n_c);
    par_c = [];
    par_scalecorr = par_s;
else
    par_s = reshape(par(n_c*n_d+1:(end-1))', 1, n_d, n_c);
    par_c = par(size(par,1));
    for i = 1:size(par_s,3)
    par_scalecorr(:,:,i) = diag(par_s(:,:,i))+flip(diag(prod(par_s(:,:,i))*par_c*ones(1,size(x,2))));
%     [~,b] = chol(par_scalecorr(:,:,i))
% 	par_scalecorr(:,:,i)
    end
end
%par_l1 = par(1:n_c);
%par_l2 = par((n_c+1):n_c*n_d);
%par_s1 = reshape(par(n_c*n_d+1:n_c*n_d+n_c), 1, 1, n_c);
%par_s2 = reshape(par(n_c*n_d+n_c+1:end), 1, 1, n_c);
%gm1 = gmdistribution(par_l1, par_s1);
%gm2 = gmdistribution(par_l2, par_s2);
%fun_pdf_val = pdf(gm1, x(:,1)).*pdf(gm2, x(:,2));
gm1 = gmdistribution(par_l, par_scalecorr);
fun_pdf_val = pdf(gm1, x);
end
