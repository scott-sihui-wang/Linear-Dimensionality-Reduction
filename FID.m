function proj=FID(means,covs,prior,method,dim,dim_r,m_iter,tol,solver)
    mu1=means(:,1,1);
    mu2=means(:,1,2);
    cov1=covs(:,:,1);
    cov2=covs(:,:,2);
    transformer=cov1^(-1/2);
    cov=transformer*cov2*transformer;
    X=transformer*(mu2-mu1);
    prev=zeros(dim,dim);
    proj=zeros(dim,dim_r);    
    options.maxiter=m_iter;
    options.tolgradnorm=tol;
    
    for i=1:dim_r
        manifold=spherefactory(dim,1);
        problem.M=manifold;
        if method=="KL"
            problem.costgrad = @(M) IKL(X,cov,prev,i-1,M);
        else if method=="rKL"
                problem.costgrad = @(M) IrKL(X,cov,prev,i-1,M);
            else if method=="SKL"
                    problem.costgrad = @(M) ISKL(X,cov,prev,i-1,M);
                else if method=="H"
                        problem.costgrad = @(M) IH(X,cov,prev,i-1,M);
                    else if method=="F"
                            problem.costgrad = @(M) IG(X,cov,prev,i-1,M);
                        else if method=="TV"
                                problem.costgrad = @(M) ITV(X,cov,prior,prev,i-1,M);
                            end
                        end
                    end
                end
            end
        end
        y=rand_start_point(prev,i-1,dim);
        if solver==1
            [x,xcost,info,options]=steepestdescent(problem,y,options);
        else if solver==2
                [x,xcost,info,options]=trustregions(problem,y,options);
            end
        end
        prev(:,i)=x;
        proj(:,i)=x;
    end
    proj=transformer*proj;
end

function [f,g]=IG(mu,cov,sol,n,M)
    f=-sqrt(M'*cov*M)-(1+0.5*mu'*M*M'*mu)/sqrt(M'*cov*M);
    egrad=-(M'*mu)/sqrt(M'*cov*M)*mu+(1+0.5*mu'*M*M'*mu)/(sqrt(M'*cov*M))^3*cov*M-cov*M/sqrt(M'*cov*M);
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function [f,g]=IKL(mu,cov,sol,n,M)
    f=-0.5*log(M'*cov*M)-0.5*(1+mu'*M*M'*mu)/(M'*cov*M);
    egrad=-(cov*M)/(M'*cov*M)-(M'*mu)*mu/(M'*cov*M)+(1+mu'*M*M'*mu)*cov*M/(M'*cov*M)^2;
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function [f,g]=IrKL(mu,cov,sol,n,M)
    f=0.5*log(M'*cov*M)-0.5*M'*cov*M-0.5*mu'*M*M'*mu;
    egrad=(cov*M)/(M'*cov*M)-cov*M-(M'*mu)*mu;
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function [f,g]=ISKL(mu,cov,sol,n,M)
    f=-0.5*M'*cov*M-0.5*mu'*M*M'*mu-0.5*(1+mu'*M*M'*mu)/(M'*cov*M);
    egrad=-cov*M-(mu'*M)*mu-(mu'*M)*mu/(M'*cov*M)+(1+mu'*M*M'*mu)*cov*M/(M'*cov*M)^2;
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function [f,g]=IH(mu,cov,sol,n,M)
    f=-log(M'*cov*M+1)+0.5*log(M'*cov*M)-0.5*(mu'*M*M'*mu)/(1+M'*cov*M);
    egrad=-(mu'*M*mu+2*cov*M)/(1+M'*cov*M)+cov*M/(M'*cov*M)+(mu'*M*M'*mu)*cov*M/(1+M'*cov*M)^2;
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function [f,g]=ITV(mu,cov,prior,sol,n,M)
    m=M'*mu;
    s2=M'*cov*M;
    s=sqrt(s2);
    k=prior(2)/prior(1);
    if s==1
        if m==0
            f=-abs(k-1);
            c1=0;
            c2=0;
        else
            x=m/2-s2/m*log(k/s);
            xx=(x-m)/s;
            if m>0
                f=-k*LF(xx,Inf)+LF(x,Inf)-LF(-Inf,x)+k*LF(-Inf,xx);
                c1=-LM(xx,Inf)/s;
                c2=-LS(xx,Inf)/s;
            else if m<0
                    f=-k*LF(-Inf,xx)+LF(-Inf,x)-LF(x,Inf)+k*LF(xx,Inf);
                    c1=-LM(-Inf,xx)/s;
                    c2=-LS(-Inf,xx)/s;
                end
            end
        end
    else
        d=m*m+2*(1-s2)*log(k/s);
        if d<=0
            f=-abs(k-1);
            c1=0;
            c2=0;
        else
            sq_d=sqrt(d);
            x1=(m-s*sq_d)/(1-s2);
            x2=(m+s*sq_d)/(1-s2);
            xx1=(x1-m)/s;
            xx2=(x2-m)/s;
            if s>1
                f=-LF(x2,x1)+k*LF(xx2,xx1)-k*LF(-Inf,xx2)+LF(-Inf,x2)-k*LF(xx1,Inf)+LF(x1,Inf);
                c1=LM(xx2,xx1)/s;
                c2=LS(xx2,xx1)/s;
            else
                f=-k*LF(xx1,xx2)+LF(x1,x2)-LF(-Inf,x1)+k*LF(-Inf,xx1)-LF(x2,Inf)+k*LF(xx2,Inf);
                c1=-LM(xx1,xx2)/s;
                c2=-LS(xx1,xx2)/s;
            end
        end
    end
    egrad=2*k*c1*mu+2*k/s*c2*cov*M;
    g=egrad-(egrad'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))*sol(:,i);
    end
end

function y=rand_start_point(sol,n,dim)
    y=randn(dim,1);
    for i=1:n
        y=y-(y'*sol(:,i))*sol(:,i);
    end
    y=y/sqrt(y'*y);
end
   
function y=LF(l,u)
    fun=@(x) 1/sqrt(2*pi)*exp(-1/2*x.*x);
    y=integral(fun,l,u);
end

function y=LM(l,u)
    fun=@(x) 1/sqrt(2*pi)*exp(-1/2*x.*x).*x;
    y=integral(fun,l,u);
end

function y=LS(l,u)
    fun=@(x) 1/sqrt(2*pi)*exp(-1/2*x.*x).*(x.*x-1);
    y=integral(fun,l,u);
end