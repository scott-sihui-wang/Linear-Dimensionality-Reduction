function [M] = dim_reduction_2(method,reduced_dim,transform,mu_1,cov_1,mu_2,cov_2)
    dim=length(mu_1);
    if nargin==5
        A=eye(dim);
        B=cov_1;
        X=mu_1;
    else if nargin==7
            if transform==false
                A=cov_1;
                B=cov_2;
                X=mu_2-mu_1;
            else
                transformer=cov_1^(-1/2);
                A=eye(dim);
                B=transformer*cov_2*transformer;
                X=transformer*(mu_2-mu_1);
            end
        end
    end
    manifold=grassmannfactory(dim,reduced_dim);
    problem.M=manifold;
    if method=="H"
        C=(A+B)/2;
        problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
        problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);
        [x, xcost, info, options] = trustregions(problem,[],options);
    end
end

