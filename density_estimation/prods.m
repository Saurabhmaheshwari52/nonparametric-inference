function [prod_mat] = prods(mat, n_dim)
prod_mat = zeros(n_dim,n_dim);
for i = 1:(n_dim-1)
    for j = (i+1):n_dim
        prod_mat(i,j) = mat(i)*mat(j);
    end
end
prod_mat = prod_mat+prod_mat';
end



