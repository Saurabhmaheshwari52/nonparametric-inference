function [G,diff_mat]=gmmjacob_num_gm(f,coef,x,param, n_comp, n_dim)
diff_mat=numeric_jacob(@(param)f(param, x, n_comp, n_dim), param);
G=coef'*diff_mat;
end