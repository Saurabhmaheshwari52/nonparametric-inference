% Correlated normal
clear
method = 'gm';
% fast function test
% fun_map = @(x) [x(:,1)+x(:,2)+x(:,3)+x(:,4)+1,...
%     (x(:,1)-1).^2,(x(:,2)-1).^2,(x(:,3)-1).^2,(x(:,4)-1).^2,... 
%     x(:,1).*x(:,2).*x(:,3).*x(:,4),...
%     x(:,1)./(x(:,2).*x(:,3).*x(:,4)), x(:,2)./(x(:,1).*x(:,3).*x(:,4)),...
%     x(:,3)./(x(:,2).*x(:,1).*x(:,4)), x(:,4)./(x(:,1).*x(:,2).*x(:,3))];
fun_map = @(x)fun_map_gm(x);
num_dim_x_full = 10;
num_corr = num_dim_x_full*(num_dim_x_full-1)/2;
% true density
par_loc = [ones(1, num_dim_x_full/2)*0.8,ones(1, num_dim_x_full/2) ;
    ones(1,num_dim_x_full/2)*1.2, ones(1,num_dim_x_full/2)*1.4];
par_scale= cat(3,ones(1, num_dim_x_full)*0.01,ones(1, num_dim_x_full)*0.02);
rng(24)
par_corr = rand(1, num_corr);
par = [par_loc(:); par_scale(:); par_corr(:)];
n_comp = size(par_loc,1);
n_dim = size(par_loc,2);
is_corr = true;
fun_pdf = @(par,x,n_comp,n_dim,is_corr) fun_pdf_gm(par,x,n_comp,n_dim,is_corr);

% generate data points
corr_mat = zeros(n_dim,n_dim);
corr_mat(triu(ones(n_dim),1) == 1) = par_corr;
corr_mat = corr_mat+corr_mat';
rng(1);
for i = 1:size(par_scale,3)
    prod_mat = prods(par_scale(:,:,i), n_dim);
    par_scalecorr(:,:,i) = diag(par_scale(:,:,i))+corr_mat.*prod_mat;
end
gmdist = gmdistribution(par_loc, par_scalecorr);
obs_x_pool_full = random(gmdist, 1e6);
obs_y_pool_full = fun_map(obs_x_pool_full);
% naive integration
fun_pdf_smp = @(x) 1/(2^num_dim_x_full);

% generate numerical sampling points
smp_x_full = [];
for i = 1:n_dim
    rng(i)
    smp_x_full = [smp_x_full, rand(1e6,1)*2];
end
smp_y_full = fun_map(smp_x_full);
if ~exist('pdf_aux','var')
    pdf_aux=ones(size(smp_x_full,1),1);
end

% estimation specification
num_sample=200;
num_obs = 250;
num_smp = 1e4;
num_smp_iter = 1e3;
num_sample_test=25;

% only consider observed y
num_dim_y_full = size(obs_y_pool_full,2);
indices_y=linspace(1,num_dim_y_full,num_dim_y_full);num_dim_y=length(indices_y); 
density_data;

% only consider selected moments
num_moment=2;
filter_moment_select= @(my)(sum(my,2)<=10);
density_moment;

% coefficient computation
density_coef;
% only consider estimated x
indices_x=linspace(1,num_dim_x_full,num_dim_x_full);num_dim_x=length(indices_x);
density_estimate_p;
