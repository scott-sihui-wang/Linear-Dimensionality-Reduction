cov=xlsread('D:\Sigma_1.xlsx');
mu=xlsread('D:\Mean_1.xlsx');
[eig_vec,lambda]=eig(cov);
alpha=eig_vec\mu;

n=100;
d=20;

prev=zeros(100,100);

proj=zeros(100,100);
r_sigma=zeros(100,1);
r_mean=zeros(100,1);
dist=zeros(100,1);

options.maxiter=1000;
options.tolgradnorm=1e-6;

t=0;

for dim=1:d
    
    manifold=spherefactory(n,1);
    problem.M=manifold;
    
    problem.costgrad = @(M) IG(alpha,lambda,prev,dim-1,M);
    
    y=rand_start_point(prev,dim-1);
    
    [x, xcost, info, options] = steepestdescent(problem,y,options);
    
    prev(:,dim)=x;
    
    proj(:,dim)=eig_vec*x;
    r_sigma(dim)=sqrt(proj(:,dim)'*cov*proj(:,dim));
    r_mean(dim)=proj(:,dim)'*mu;
    d1=sqrt((r_mean(dim))^2/2+(1+r_sigma(dim))^2);
    d2=sqrt((r_mean(dim))^2/2+(1-r_sigma(dim))^2);
    dist(dim)=(d1+d2)/(d1-d2);
    
    t=t+info(length(info)).time;
    
end

function [f,g]=IG(alpha,lambda,sol,n,M)
    f=-(1+0.5*alpha'*M*M'*alpha)/sqrt(M'*lambda*M)-sqrt(M'*lambda*M);
    egrad=-(M'*alpha)/sqrt(M'*lambda*M)*alpha+(1+1/2*alpha'*M*M'*alpha)/(sqrt(M'*lambda*M))^3*lambda*M-lambda*M/sqrt(M'*lambda*M);
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