% Correlated normal
clear
method = 'gm';
% fast function test
fun_map = @(x) [x(:,1)+x(:,2)+1,(x(:,1)-1).^2,(x(:,2)-1).^2, x(:,1).*x(:,2)];
num_dim_y_full = 4;num_dim_x_full = 2;

% true density
par_loc = [[0.8,0.8];[1.2,1.2]];
par_scale= cat(3,[0.01,0.02],[0.01,0.02]);
par_corr = 0.5;
par = [par_loc(:); par_scale(:); par_corr(:)];
n_comp = size(par_loc,1);
n_dim = size(par_loc,2);
is_corr = true;
fun_pdf = @(par, x, n_comp, n_dim, is_corr) fun_pdf_gm(par,x,n_comp, n_dim, is_corr);

% generate data points
rng(1);
for i = 1:size(par_scale,3)
    par_scalecorr(:,:,i) = diag(par_scale(:,:,i))+flip(diag(prod(par_scale(:,:,i))*...
    par_corr*ones(1,num_dim_x_full)));
end
gmdist = gmdistribution(par_loc, par_scalecorr);
obs_x_pool_full = random(gmdist, 1e6);
obs_y_pool_full = fun_map(obs_x_pool_full);
% naive integration
fun_pdf_smp = @(x) 1/(2^num_dim_x_full);

% generate numerical sampling points
rng(2);
smp_x_full = [rand(1e6,1)*2,rand(1e6,1)*2];
smp_y_full = fun_map(smp_x_full);
if ~exist('pdf_aux','var')
    pdf_aux=ones(size(smp_x_full,1),1);
end

% estimation specification
num_sample=200;
num_obs = 250;
num_smp = 1e5;
num_smp_iter = 1e4;
num_sample_test=50;

% only consider observed y
indices_y=[1,2,3,4];num_dim_y=length(indices_y); 
density_data;

% only consider selected moments
num_moment=4;
filter_moment_select= @(my)(sum(my,2)<=10);
density_moment;

% coefficient computation
density_coef;

% only consider estimated x
indices_x=[1,2];num_dim_x=length(indices_x);
density_estimate_p;
