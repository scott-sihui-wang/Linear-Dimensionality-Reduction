function [proj,transformer]=SLDA(means,covs,method,n_dim,r_dim,transformed,m_iter,tol,solver)
    X=means(:,1,2)-means(:,1,1);
    if transformed
        A=eye(n_dim);
        transformer=covs(:,:,1)^(-1/2);
        B=transformer*covs(:,:,2)*transformer;
        X=transformer*X;
    else
        transformer=0;
        A=covs(:,:,1);
        B=covs(:,:,2);
    end
    manifold=grassmannfactory(n_dim,r_dim);
    problem.M=manifold;
    if method=="KL"
        problem.cost = @(M) -0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.5*trace((M'*A*M)\(M'*B*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
        problem.egrad = @(M) -(A*M)/(M'*A*M)+(B*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M);
    else
        if method=="rKL"
            problem.cost = @(M) -0.5*logdet(M'*B*M)+0.5*logdet(M'*A*M)-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*B*M)*(M'*X);
            problem.egrad = @(M) -(B*M)/(M'*B*M)+(A*M)/(M'*A*M)-(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);
        else
            if method=="SKL"
            problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
            problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);
            else
                if method=="H"
                    C=(A+B)/2;
                    problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
                    problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);
                end
            end
        end
    end
    options.maxiter=m_iter;
    options.tolgradnorm=tol;
    if solver==1
        [proj, xcost, info, options] = steepestdescent(problem);
    else
        if solver==2
            [proj, xcost, info, options] = trustregions(problem,[],options);
        end
    end
end