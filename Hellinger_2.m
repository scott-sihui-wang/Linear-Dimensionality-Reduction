n=100;
d=1;

A=eye(n);
B=xlsread('D:\Sigma_100.xlsx');

%X=(x1-x2)';

X=10000*randn(100,1);

manifold=grassmannfactory(n,d);
problem.M=manifold;

C=(A+B)/2;

problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);

checkgradient(problem);

options.maxiter=5000;
options.tolgradnorm=1e-9;

[x, xcost, info, options] = trustregions(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');