%n = 784;
%d = 3;
%m=xlsread('means.xlsx');
%c=xlsread('cov.xlsx');
%e=c';
%A=reshape(e(2,:),[784,784]);
%B=reshape(e(3,:),[784,784]);

%X1=randn(n,1);
%X2=randn(n,1);
%X=X1-X2;
%x=m(1,:)-m(4,:);
%X=x';
%X(10)=1.0;

%n=100;
%d=10;

%A=xlsread('D:\Sigma1.xlsx');
%B=xlsread('D:\Sigma2.xlsx');
%X=xlsread('Mean_diff_positive.xlsx');

n=7;
d=1;

%A=xlsread('icov2.xlsx');
%B=xlsread('icov3.xlsx');

%x1=xlsread('imean2.xlsx');
%x2=xlsread('imean3.xlsx');

X=(x1-x2)';

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
%problem.egrad = @(M) -(A*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*A*M)+0.5*(B*M)*(M'*A*M)/(M'*B*M)/(M'*B*M)-(B*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*B*M)+0.5*(A*M)*(M'*B*M)/(M'*A*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+0.5*(A*M)*(M'*(X*X')*M)/(M'*A*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*(X*X')*M)-(X*X')*M/(M'*B*M)+0.5*(B*M)*(M'*(X*X')*M)/(M'*B*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*(X*X')*M);
problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);

checkgradient(problem);
%checkdiff(problem);

options.tolgradnorm=1e-12;

%[x, xcost, info, options] = steepestdescent(problem);
[x, xcost, info, options] = trustregions(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
