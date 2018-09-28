function [pdf_joint]=copula_nataf(method, x ,par_loc_x,par_scale_x,par_corr_x)
    norm_x = zeros(size(x));
    if strcmp(method, 'log')
        for i=1:size(x,2), norm_x(:,i)=cdf2normx(logncdf(x(:,i),par_loc_x(i),par_scale_x(i))); end
        pdf_joint = lognpdf(x(:,1),par_loc_x(1),par_scale_x(1)).*lognpdf(x(:,2),par_loc_x(2),par_scale_x(2))...
                    ./ prod(normpdf(norm_x),2) .* mvnpdf(norm_x,zeros(1,size(x,2)),diag(ones(1,size(x,2)))+flip(diag(ones(1,size(x,2))*par_corr_x)));
    end
    if strcmp(method, 'gam')
        for i=1:size(x,2), norm_x(:,i)=cdf2normx(gamcdf(x(:,i),par_loc_x(i),par_scale_x(i))); end
        pdf_joint = gampdf(x(:,1),par_loc_x(1),par_scale_x(1)).*gampdf(x(:,2),par_loc_x(2),par_scale_x(2))...
                    ./ prod(normpdf(norm_x),2) .* mvnpdf(norm_x,zeros(1,size(x,2)),diag(ones(1,size(x,2)))+flip(diag(ones(1,size(x,2))*par_corr_x)));
    end
end