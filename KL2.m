n = 100;
d = 2;
A=xlsread('D:\Sigma2.xlsx');
B=xlsread('D:\Sigma1.xlsx');
%B=eye(n);

%P=B^(-1/2);
%A=P*A*P;
%B=P*B*P;

%manifold=stiefelfactory(n,d);
manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*logdet(M'*B*M)+0.5*logdet(M'*A*M)-0.5*trace((M'*B*M)\(M'*A*M));
%problem.egrad = @(M) -(B*M)/(M'*B*M)+(A*M)/(M'*A*M)-(A*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*A*M)+0.5*(B*M)*(M'*A*M)/(M'*B*M)/(M'*B*M);
problem.egrad = @(M) -(B*M)/(M'*B*M)+(A*M)/(M'*A*M)-(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M);

%problem.cost = @(M) -0.5*trace((M'*B*M)\(M'*A*M));
%problem.egrad = @(M) -(A*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*A*M)+0.5*(B*M)*(M'*A*M)/(M'*B*M)/(M'*B*M);

%problem.cost = @(M) -0.5*logdet(M'*B*M)+0.5*logdet(M'*A*M);
%problem.egrad = @(M) -B*M/(M'*B*M)+A*M/(M'*A*M);

checkgradient(problem);

%options.verbosity=0;

[x, xcost, info, options] = trustregions(problem,P*x_o);
%[x, xcost, info, options] = steepestdescent(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
