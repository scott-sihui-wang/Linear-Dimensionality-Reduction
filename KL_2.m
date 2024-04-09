n=100;
d=1;

A=eye(n);
B=xlsread('D:\Sigma_2.xlsx');

%X=(x1-x2)';

X=xlsread('D:\Mean_2.xlsx');

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.5*trace((M'*A*M)\(M'*B*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
problem.egrad = @(M) -(A*M)/(M'*A*M)+(B*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M);

checkgradient(problem);

options.maxiter=5000;
options.tolgradnorm=1e-9;

[x, xcost, info, options] = trustregions(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');