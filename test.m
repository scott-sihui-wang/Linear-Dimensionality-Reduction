n = 100;
d = 10;
A=xlsread('D:\Sigma1.xlsx');

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -trace(M'*A*M);
problem.egrad = @(M) -2*A'*M;

checkgradient(problem);

%options.verbosity=0;
%options.tolCost=1e-12;
[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
