function [list] = fun_map_gm(x)
n_dim = size(x,2);
list = [];
prod_all = prod(x, 2);
for i = 1:n_dim
    list = [list, (x(:,i)-1).^2, x(:,i).^2./prod_all];
end
list = [list, prod_all, sum(x,2)+1];

    