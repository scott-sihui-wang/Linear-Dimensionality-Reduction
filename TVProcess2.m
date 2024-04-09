function [err,ttl,mis_rate,p_err]=TVProcess2(filename,ext,sorting_column,sorting_method,row_bound,column_bound,dim,r_dim,m_iter,tol,solver)
    [means,covs,prior,data,row_bound,column_bound] = preprocessing(filename,ext,sorting_column,sorting_method,2,row_bound,column_bound,dim,true,false);
    proj=FID(means,covs,prior,"TV",dim,r_dim,m_iter,tol,solver);
    [r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,proj]=criteria(means(:,1,1),means(:,1,2),covs(:,:,1),covs(:,:,2),proj,0);
    [err,ttl,mis_rate,p_err]=hypTest(data,row_bound,column_bound,proj,r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,[1 1]);
end
