a1 = load('est_total_normal_0.2_diag_moment_4.mat');
a1 = a1.est_total;
a2 = load('est_total_normal_0.3_diag_moment_4.mat');
a2 = a2.est_total;
a3 = load('est_total_normal_0.4_diag_moment_4.mat');
a3 = a3.est_total;
index = reshape(1:6,2,3).';
for ind = 1:6
    subplot(3,2,index(ind))
    plot_density(a1,a2,a3,ind,'normal')
end
% hl = legend('scale = 1.2','scale = 1.35','scale = 1.5');
newPosition = [0.85 0.7 0.2 0.2];
newUnits = 'normalized';
set(legend({'scale = 0.2','scale = 0.3','scale = 0.4'},'Orientation','Vertical'),...
    'Position', newPosition,'Units', newUnits);
suptitle('Normal Observed Distribution (Moment bound = 4) only diagonal in gmm matrix')
