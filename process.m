function [err,ttl,mis_rate,p_err]=process(filename,ext,sorting_column,sorting_method,num_of_class,row_bound,column_bound,prior_enable,truncate_class,dim,r_dim,divergence,transformed,m_iter,tol,solver)
    [means,covs,prior,data,row_bound,column_bound] = preprocessing(filename,ext,sorting_column,sorting_method,num_of_class,row_bound,column_bound,dim,prior_enable,truncate_class);
    [proj,transformer]=SLDA(means,covs,divergence,dim,r_dim,m_iter,tol,solver);
    [r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,proj]=criteria(means(:,1,1),means(:,1,2),covs(:,:,1),covs(:,:,2),proj,transformer);
    [err,ttl,mis_rate,p_err]=hypTest(data,row_bound,column_bound,proj,r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,prior);
end
