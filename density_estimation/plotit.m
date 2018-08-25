a1 = load('est_total_lognormal_0.2.mat');
a1 = a1.est_total;
a2 = load('est_total_lognormal_0.3.mat');
a2 = a2.est_total;
a3 = load('est_total_lognormal_0.4.mat');
a3 = a3.est_total;

for ind = 1:6
    subplot(4,2,ind)
    plot_density(a1,a2,a3,ind,'lognormal')
end
% hl = legend('scale = 1.2','scale = 1.35','scale = 1.5');
newPosition = [0.4 0.1 0.2 0.2];
newUnits = 'normalized';
set(legend('scale = 0.2','scale = 0.3','scale = 0.4'),'Position', newPosition,'Units', newUnits);
