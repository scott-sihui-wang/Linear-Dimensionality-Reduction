cov=xlsread('D:\Sigma_1.xlsx');
mu=xlsread('D:\Mean_1.xlsx');
[eig_vec,lambda]=eig(cov);
alpha=eig_vec\mu;

n=100;
d=1;

manifold=spherefactory(n,d);
problem.M=manifold;

problem.costgrad = @(M) IG_D1(alpha,lambda,M);

%checkgradient(problem);

options.maxiter=6000;
options.tolgradnorm=1e-12;

[x1, xcost, info, options] = steepestdescent(problem,[],options);

%figure;
%semilogy([info.iter], [info.gradnorm], '.-');
%xlabel('Iteration number');
%ylabel('Norm of the gradient of f');

proj1=eig_vec*x1;

r_sigma1=sqrt(proj1'*cov*proj1);
r_mean1=proj1'*mu;

d11=sqrt((r_mean1)^2/2+(1+r_sigma1)^2);
d12=sqrt((r_mean1)^2/2+(1-r_sigma1)^2);
dist1=(d11+d12)/(d11-d12);

manifold=spherefactory(n,d);
problem.M=manifold;

problem.costgrad = @(M) IG_D2(alpha,lambda,x1,M);

%checkgradient(problem);

options.maxiter=6000;
options.tolgradnorm=1e-12;

y=randn(100,1);
y=y-(y'*x1)*x1;
y=y/sqrt(y'*y);

[x2, xcost, info, options] = steepestdescent(problem,y,options);

%figure;
%semilogy([info.iter], [info.gradnorm], '.-');
%xlabel('Iteration number');
%ylabel('Norm of the gradient of f');

proj2=eig_vec*x2;

r_sigma2=sqrt(proj2'*cov*proj2);
r_mean2=proj2'*mu;

d21=sqrt((r_mean2)^2/2+(1+r_sigma2)^2);
d22=sqrt((r_mean2)^2/2+(1-r_sigma2)^2);
dist2=(d21+d22)/(d21-d22);

function [f,g]=IG_D1(alpha,lambda,M)
    f=-(1+0.5*alpha'*M*M'*alpha)/sqrt(M'*lambda*M)-sqrt(M'*lambda*M);
    egrad=-(M'*alpha)/sqrt(M'*lambda*M)*alpha+(1+1/2*alpha'*M*M'*alpha)/(sqrt(M'*lambda*M))^3*lambda*M-lambda*M/sqrt(M'*lambda*M);
    g=egrad-(egrad'*M)/(M'*M)*M;
end
function [f,g]=IG_D2(alpha,lambda,x,M)
    f=-(1+0.5*alpha'*M*M'*alpha)/sqrt(M'*lambda*M)-sqrt(M'*lambda*M);
    egrad=-(M'*alpha)/sqrt(M'*lambda*M)*alpha+(1+1/2*alpha'*M*M'*alpha)/(sqrt(M'*lambda*M))^3*lambda*M-lambda*M/sqrt(M'*lambda*M);
    g1=egrad-(egrad'*M)/(M'*M)*M;
    g=g1-(g1'*x)/(x'*x)*x;
end