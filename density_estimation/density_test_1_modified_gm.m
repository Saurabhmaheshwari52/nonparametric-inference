clear

% fast function test
% fun_map = @(x) [exp(x(:,1)+x(:,2)+1),(x(:,1)-1).^2,(x(:,2)-1).^2,...
%     exp(x(:,2)+x(:,3)+1),(x(:,3)-1).^2];
fun_map = @(x) [(x(:,1)-1).^2, (x(:,2)-1).^2,...
    (x(:,1)+1).^2.*(x(:,2)+1).^2, (x(:,1).*x(:,2)),...
    x(:,1).^1.5+x(:,2).^1.5,...
    x(:,1)./(1e-3+x(:,2))];
num_dim_y_full = 6;num_dim_x_full = 2;
method = 'gm';
% true density
par_loc = [[1,1];[2,2]];
par_scale= cat(3,[0.01,0.02],[0.01,0.02]);
par_corr = [];
par = [par_loc(:); par_scale(:)];
n_comp = size(par_loc,1);
n_dim = size(par_loc,2);
fun_pdf = @(par, x, n_c, n_d)fun_pdf_gm(par, x, n_c, n_d);

% generate data points
rng(1);
gm_loc_1 = par_loc(:,1);
gm_loc_2 = par_loc(:,2);
gm_scale_1 = reshape(par_scale(:,:,1), 1, 1, n_comp);
gm_scale_2 = reshape(par_scale(:,:,2), 1, 1, n_comp);
gm_obs_1 = gmdistribution(gm_loc_1, gm_scale_1);
gm_obs_2 = gmdistribution(gm_loc_2, gm_scale_2);
obs_x_pool_full = [random(gm_obs_1, 1e7), random(gm_obs_2, 1e7)];
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
num_obs = 50000;
num_smp = 1e6;
num_smp_iter = 1e5;
num_sample_test=50;

% only consider observed y
indices_y=[1,2,3,4,5,6];num_dim_y=length(indices_y); 
density_data;

% only consider selected moments
num_moment=4;
filter_moment_select= @(my)(sum(my,2)<=num_moment);
density_moment;

% coefficient computation
density_coef;

% only consider estimated x
indices_x=[1,2];num_dim_x=length(indices_x);
density_estimate_p;
