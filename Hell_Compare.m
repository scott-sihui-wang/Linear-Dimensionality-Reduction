A=xlsread('D:\Sigma_1.xlsx');
B=xlsread('D:\Sigma_4.xlsx');
mu1=xlsread('D:\Mean_1.xlsx');
mu2=xlsread('D:\Mean_4.xlsx');

n=100;
d=3;

X=mu2-mu1;

manifold=grassmannfactory(n,d);
problem.M=manifold;

C=(A+B)/2;

problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);

checkgradient(problem);

options.maxiter=3000;
options.tolgradnorm=1e-6;

[x, xcost, info, options] = trustregions(problem,[],options);
%[x, xcost, info, options] = steepestdescent(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

r_sigma_1=zeros(1,d);
r_sigma_2=zeros(1,d);
r_mean_1=zeros(1,d);
r_mean_2=zeros(1,d);
dist_IG=zeros(1,d);

for i=1:d
    r_sigma_1(i)=sqrt(x(:,i)'*A*x(:,i));
    r_sigma_2(i)=sqrt(x(:,i)'*B*x(:,i));
    r_mean_1(i)=x(:,i)'*mu1;
    r_mean_2(i)=x(:,i)'*mu2;
    d1=sqrt((r_mean_1(i)-r_mean_2(i))^2/2+(r_sigma_1(i)+r_sigma_2(i))^2);
    d2=sqrt((r_mean_1(i)-r_mean_2(i))^2/2+(r_sigma_1(i)-r_sigma_2(i))^2);
    dist_IG(i)=(d1+d2)/(d1-d2);
end

%r_sigma=sqrt(x'*cov*x);
%r_mean=x'*mu;

%d1=sqrt((r_mean)^2/2+(1+r_sigma)^2);
%d2=sqrt((r_mean)^2/2+(1-r_sigma)^2);
%dist=(d1+d2)/(d1-d2);