n = 100;
d = 10;
A=xlsread('D:\Sigma1.xlsx');

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*M))-0.5*trace((M'*M)\(M'*A*M));
problem.egrad = @(M) -M/(M'*A*M)+0.5*A*M/(M'*A*M)/(M'*A*M)*M'*M+0.5*A*M*M'*M/(M'*A*M)/(M'*A*M)-A*M/(M'*M);

checkgradient(problem);
checkdiff(problem);

[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
