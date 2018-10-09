% 3. Coefficient Computation
% integrate gx using monte carlo integration
weight=1./fun_pdf_smp(smp_x_full)/num_smp./pdf_aux(1:num_smp,:);
coefpure=zeros(num_smp,size(momentlist,1));
for i=1:size(momentlist,1)

    coefpure(:,i)=prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight;

end 

coefpure=coefpure.*(pdf_aux(1:num_smp,:)*ones(1,size(momentlist,1)));
%-----------------------------------------------------------------------------
tic
% generate integrated population moment
moment_pdfsmp_y = zeros(size(momentlist,1),1);
if ~exist('method')
    method = '';
end
if (strcmp(method, 'gm'))
    pdf_smp = fun_pdf(par,smp_x_full(1:num_smp,:),n_comp,n_dim,is_corr);
else
    n_comp = 0;
    n_dim = 0;
    is_corr = [];
    pdf_smp = fun_pdf(smp_x_full(1:num_smp,:),par);
end

% sample y pdf
for i=1:size(momentlist,1)
    moment_pdfsmp_y(i)=dot(prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight,pdf_smp.*pdf_aux(1:num_smp,:));
end
toc

%------------------------------------------------------------------------------
tic
 [sigmadata]=weight_obs(obs_pool_y(1:1e4,:),moment_obs_pool_y(moment_select),momentlist(moment_select,:));
toc

if (strcmp(method, 'gm'))
    sigmasimulation=weight_smp(smp_y(1:1e4,:),1./fun_pdf_smp(smp_x_full(1:1e4,:)).*fun_pdf_gm(par, smp_x_full(1:1e4,:), n_comp, n_dim, is_corr),moment_obs_pool_y(moment_select),momentlist(moment_select,:));
else
   sigmasimulation=weight_smp(smp_y(1:1e4,:),1./fun_pdf_smp(smp_x_full(1:1e4,:)).*fun_pdf(smp_x_full(1:1e4,:),par),moment_obs_pool_y(moment_select),momentlist(moment_select,:));
end
