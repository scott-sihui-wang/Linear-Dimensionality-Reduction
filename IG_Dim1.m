cov=xlsread('D:\Sigma_2.xlsx');
mu=xlsread('D:\Mean_2.xlsx');
[eig_vec,lambda]=eig(cov);
alpha=eig_vec\mu;

n=100;
d=1;

manifold=spherefactory(n,d);
problem.M=manifold;

problem.cost = @(M) -(1+0.5*alpha'*M*M'*alpha)/sqrt(M'*lambda*M)-sqrt(M'*lambda*M);
problem.egrad = @(M) -(M'*alpha)/sqrt(M'*lambda*M)*alpha+(1+1/2*alpha'*M*M'*alpha)/(sqrt(M'*lambda*M))^3*lambda*M-lambda*M/sqrt(M'*lambda*M);

checkgradient(problem);

options.maxiter=6000;
options.tolgradnorm=1e-12;

%[x, xcost, info, options] = trustregions(problem,[],options);
[x, xcost, info, options] = steepestdescent(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

proj=eig_vec*x;

r_sigma=sqrt(proj'*cov*proj);
r_mean=proj'*mu;

d1=sqrt((r_mean)^2/2+(1+r_sigma)^2);
d2=sqrt((r_mean)^2/2+(1-r_sigma)^2);
dist=(d1+d2)/(d1-d2);