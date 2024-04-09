n = 100;
d = 10;
A=xlsread('D:\Sigma2.xlsx');
X=zeros([n,1]);
X(50)=1.0;

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*M))-0.5*trace((M'*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*M)*(M'*X);
%-0.5*trace(MAM\MBM)-0.5*trace(MBM\MAM)-0.5*XM/MAM*XM'-0.5*XM/MBM*XM'
%problem.egrad = @(M) -M/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*M)+0.5*(A*M)*(M'*M)/(M'*A*M)/(M'*A*M)-(A*M)/(M'*M)-(X*X')*M/(M'*A*M)+0.5*(A*M)*(M'*(X*X')*M)/(M'*A*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*(X*X')*M)-(X*X')*M/(M'*M)+0.5*(M*M')*(X*X')*M/(M'*M)/(M'*M)+0.5*M/(M'*M)/(M'*M)*(M'*(X*X')*M);
%-AM/MBM+BM/MBM*MAM/MBM-BM/MAM+AM/MAM*MBM/MAM+AM/MAM*MXXM/MAM-XXM/MAM+BM/MBM*MXXM/MBM-XXM/MBM
problem.egrad = @(M) -(A*M)/(M'*M)+M/(M'*M)*(M'*A*M)/(M'*M)-M/(M'*A*M)+(A*M)/(M'*A*M)*(M'*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+M/(M'*M)*(M'*(X*X')*M)/(M'*M)-(X*X')*M/(M'*M);

checkgradient(problem);

%[x, xcost, info, options] = steepestdescent(problem);
[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
