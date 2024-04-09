n = 100;
d = 10;

A=xlsread('D:\Sigma1.xlsx');
%X=zeros(n,1);
%X=xlsread('mean_diff.xlsx');
X=xlsread('Mean_diff_positive.xlsx');
%X=xlsread('mean_diff.xlsx');
X=100*X;

manifold=grassmannfactory(n,d);
%manifold=stiefelfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*logdet(M'*A*M)+0.5*logdet(M'*M)-0.5*trace((M'*A*M)\(M'*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
problem.egrad = @(M) -(A*M)/(M'*A*M)+M/(M'*M)-M/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*M)+0.5*(A*M)*(M'*M)/(M'*A*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M);

checkgradient(problem);

%options.verbosity=0;

[x, xcost, info, options] = trustregions(problem);
%[x, xcost, info, options] = steepestdescent(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
