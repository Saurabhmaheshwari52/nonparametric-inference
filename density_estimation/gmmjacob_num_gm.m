function [G,diff_mat]=gmmjacob_num_gm(f,coef,x,param, n_comp, n_dim, is_corr)
diff_mat=numeric_jacob(@(param)f(param, x, n_comp, n_dim, is_corr), param);
G=coef'*diff_mat;
end