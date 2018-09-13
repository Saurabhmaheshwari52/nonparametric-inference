% For bimodal distributions
function [fun_pdf_val] = fun_pdf_gm(par, x, n_c, n_d)
% par_l = reshape(par(1:n_c*n_d), n_c, n_d);
% par_s = reshape(par(n_c*n_d+1:end)', 1, n_d, n_c);
par_l1 = par(1:n_c);
par_l2 = par((n_c+1):n_c*n_d);
par_s1 = reshape(par(n_c*n_d+1:n_c*n_d+n_c), 1, 1, n_c);
par_s2 = reshape(par(n_c*n_d+n_c+1:end), 1, 1, n_c);
gm1 = gmdistribution(par_l1, par_s1);
gm2 = gmdistribution(par_l2, par_s2);
fun_pdf_val = pdf(gm1, x(:,1)).*pdf(gm2, x(:,2));
end
