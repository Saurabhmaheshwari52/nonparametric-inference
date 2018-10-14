function [val,grad] = gmmobj(param,method,x,coef,moment,weight,fun_pdf, n_comp, n_dim, is_corr)
% n_comp
% n_dim
if ~strcmp(method, 'gm')
%     num_dim_x=size(x,2);
    diff=(fun_pdf(x,param)'*coef)'-moment;
else
%     size(coef'*fun_pdf(param, x, n_comp, n_dim))
%     size(moment)
    diff = coef'*fun_pdf(param, x, n_comp, n_dim, is_corr)-moment;
end
 
% Braess
% diff=((ones(size(x))*1/param.*(x<=param))'*coef)'-moment;
val=diff'*weight*diff;
%grad = 1;
grad=2*diff'*weight*gmmjacob_num_gm(fun_pdf, coef, x, param, n_comp, n_dim, is_corr);
end

