n = 100;
d = 10;
%C1=xlsread('D:\Sigma2.xlsx');
%C2=xlsread('D:\Sigma1.xlsx');
%A=C1^(-1/2)*C2*C1^(-1/2);
%A=C1\C2;
%A=C1/C2;

A=xlsread('D:\Sigma1.xlsx');

manifold=grassmannfactory(n,d);
%manifold=stiefelfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*logdet(M'*A*M)+0.5*logdet(M'*M)-0.5*trace((M'*A*M)\(M'*M));
problem.egrad = @(M) -(A*M)/(M'*A*M)+M/(M'*M)-M/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*M)+0.5*(A*M)*(M'*M)/(M'*A*M)/(M'*A*M);

%problem.cost = @(M) -0.5*logdet(M'*A*M)-0.5*trace(inv(M'*A*M));
%problem.egrad = @(M) -(A*M)/(M'*A*M)+(A*M)/(M'*A*M)/(M'*A*M);

checkgradient(problem);

%options.verbosity=0;

[x, xcost, info, options] = trustregions(problem);
%[x, xcost, info, options] = steepestdescent(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
