% Correlated normal
clear
method = 'log';
% fast function test
fun_map = @(x) [exp(x(:,1)+x(:,2)+1),(x(:,1)-1).^2,(x(:,2)-1).^2];
num_dim_y_full = 3;num_dim_x_full = 2;

% true density
par_loc = ones(1,num_dim_x_full)*-.1;
par_scale= 0.2*ones(1,num_dim_x_full);
par_corr = 0.5;
par = [par_loc,par_scale,par_corr];
fun_pdf = @(x,par) copula_nataf(method, x,par(1:num_dim_x_full),par(num_dim_x_full+1:num_dim_x_full*2),par(num_dim_x_full*2+1:end));

% generate data points
rng(1);
obs_xmvn_pool_full = mvnrnd(zeros(1,num_dim_x_full),diag(ones(1,num_dim_x_full).^2)+flip(diag(par_corr*ones(1,num_dim_x_full))),1e6);
obs_x_pool_full = normcdf(obs_xmvn_pool_full);
for i =1:num_dim_x_full, obs_x_pool_full(:,i)=logninv(obs_x_pool_full(:,i),par_loc(i),par_scale(i)); end
obs_y_pool_full = fun_map(obs_x_pool_full);

% naive integration
fun_pdf_smp = @(x) 1/(2.5^num_dim_x_full);

% generate numerical sampling points
rng(2);
smp_x_full = [rand(1e6,1)*2.5,rand(1e6,1)*2.5];
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
indices_y=[1,2,3];num_dim_y=length(indices_y); 
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
