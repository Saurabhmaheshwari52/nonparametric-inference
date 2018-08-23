clear

% fast function test
fun_map = @(x) [exp(x(:,1)+x(:,2)+1),(x(:,1)-1).^2,(x(:,2)-1).^2,...
    exp(x(:,2)+x(:,3)+1),(x(:,3)-1).^2];
num_dim_y_full = 5;num_dim_x_full = 3;

% true density
par_loc = 0.5*ones(1,num_dim_x_full); 
par_scale= 1.35*ones(1,num_dim_x_full);
par_corr = [];
par = [par_loc,par_scale,par_corr];
fun_pdf = @(x,par) wblpdf(x(:,1),par(1),par(4)).*wblpdf(x(:,2),par(2),par(5)).*...
    wblpdf(x(:,3),par(3),par(6));

% generate data points
rng(1);
obs_x_pool_full = [wblrnd(par_loc(1),par_scale(1),1e6,1),wblrnd(par_loc(2),...
    par_scale(2),1e6,1),wblrnd(par_loc(3),par_scale(3),1e6,1)];
obs_y_pool_full = fun_map(obs_x_pool_full);

% naive integration
fun_pdf_smp = @(x) 1/(2^num_dim_x_full);

% generate numerical sampling points
rng(2);
smp_x_full = [rand(1e6,1)*2,rand(1e6,1)*2,rand(1e6,1)*2];
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
indices_y=[1,2,3,4,5];num_dim_y=length(indices_y); 
density_data;

% only consider selected moments
num_moment=4;
filter_moment_select= @(my)(sum(my,2)<=10);
density_moment;

% coefficient computation
density_coef;

% only consider estimated x
indices_x=[1,2,3];num_dim_x=length(indices_x);
moment_obs_y(:,1)=moment_pdfsmp_y;
density_estimate_p;
