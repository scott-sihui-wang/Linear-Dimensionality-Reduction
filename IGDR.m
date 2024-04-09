cov1=xlsread('D:\Sigma_1.xlsx');
mu1=xlsread('D:\Mean_1.xlsx');
cov2=xlsread('D:\Sigma_4.xlsx');
mu2=xlsread('D:\Mean_4.xlsx');

X=mu2-mu1;

n=100;
d=3;

prev=zeros(100,100);

r_sigma_1=zeros(100,1);
r_sigma_2=zeros(100,1);
r_mean_1=zeros(100,1);
r_mean_2=zeros(100,1);
dist=zeros(100,1);

options.maxiter=1000;
options.tolgradnorm=1e-6;

t=0;

for dim=1:d
    
    manifold=spherefactory(n,1);
    problem.M=manifold;
    
    problem.costgrad = @(M) IG_DR(cov1,cov2,X,prev,dim-1,M);
    
    if dim==1
        checkgradient(problem);
    end
    
    y=rand_start_point(prev,dim-1);
    
    [x, xcost, info, options] = steepestdescent(problem,y,options);
    
    prev(:,dim)=x;
    
    r_sigma_1(dim)=sqrt(prev(:,dim)'*cov1*prev(:,dim));
    r_sigma_2(dim)=sqrt(prev(:,dim)'*cov2*prev(:,dim));
    r_mean_1(dim)=prev(:,dim)'*mu1;
    r_mean_2(dim)=prev(:,dim)'*mu2;
    d1=sqrt((r_mean_1(dim)-r_mean_2(dim))^2/2+(r_sigma_1(dim)+r_sigma_2(dim))^2);
    d2=sqrt((r_mean_1(dim)-r_mean_2(dim))^2/2+(r_sigma_1(dim)-r_sigma_2(dim))^2);
    dist(dim)=(d1+d2)/(d1-d2);
    
    t=t+info(length(info)).time;
    
end

function [f,g]=IG_DR(A,B,X,sol,n,M)
    XXM=X*X'*M;
    MXXM=M'*XXM;
    AM=A*M;
    BM=B*M;
    MAM=M'*AM;
    MBM=M'*BM;
    f=-0.5*MXXM/(MAM*MBM)-MAM/MBM-MBM/MAM;
    egrad=-2*AM/MBM+2*MAM/(MBM)^2*BM-2*BM/MAM+2*MBM/(MAM)^2*AM-XXM/MAM/MBM+MXXM/(MAM)^2/MBM*AM+MXXM/MAM/(MBM)^2*BM;
    g=egrad-(egrad'*M)/(M'*M)*M;
    for i=1:n
        g=g-(g'*sol(:,i))/(sol(:,i)'*sol(:,i))*sol(:,i);
    end
end
function y=rand_start_point(sol,n)
    y=randn(100,1);
    for i=1:n
        y=y-(y'*sol(:,i))*sol(:,i);
    end
    y=y/sqrt(y'*y);
end