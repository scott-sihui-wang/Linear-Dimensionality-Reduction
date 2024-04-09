%n = 784;
%d = 3;

%m=xlsread('means.xlsx');
%c=xlsread('cov.xlsx');
%e=c';
%A=reshape(e(1,:),[784,784]);
%B=reshape(e(2,:),[784,784]);
%C=reshape(e(3,:),[784,784]);

%X1=randn(n,1);
%X2=randn(n,1);
%X=X1-X2;
%x1=m(1,:)-m(2,:);
%X1=x1';
%x2=m(2,:)-m(3,:);
%X2=x2';
%x3=m(1,:)-m(3,:);
%X3=x3';

%X(10)=1.0;

n=7;
d=3;

%A=xlsread('icov1.xlsx');
%B=xlsread('icov2.xlsx');
%C=xlsread('icov3.xlsx');
%x1=xlsread('imean1.xlsx');
%x2=xlsread('imean2.xlsx');
%x3=xlsread('imean3.xlsx');
X1=(x1-x2)';
X2=(x2-x3)';
X3=(x1-x3)';

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X1'*M)/(M'*A*M)*(M'*X1)-0.5*(X1'*M)/(M'*B*M)*(M'*X1)-0.5*trace((M'*B*M)\(M'*C*M))-0.5*trace((M'*C*M)\(M'*B*M))-0.5*(X2'*M)/(M'*B*M)*(M'*X2)-0.5*(X2'*M)/(M'*C*M)*(M'*X2)-0.5*trace((M'*A*M)\(M'*C*M))-0.5*trace((M'*C*M)\(M'*A*M))-0.5*(X3'*M)/(M'*A*M)*(M'*X3)-0.5*(X3'*M)/(M'*C*M)*(M'*X3);
%problem.egrad = @(M) -(A*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*A*M)+0.5*(B*M)*(M'*A*M)/(M'*B*M)/(M'*B*M)-(B*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*B*M)+0.5*(A*M)*(M'*B*M)/(M'*A*M)/(M'*A*M)-(X1*X1')*M/(M'*A*M)+0.5*(A*M)*(M'*(X1*X1')*M)/(M'*A*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*(X1*X1')*M)-(X1*X1')*M/(M'*B*M)+0.5*(B*M)*(M'*(X1*X1')*M)/(M'*B*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*(X1*X1')*M)-(B*M)/(M'*C*M)+0.5*(C*M)/(M'*C*M)/(M'*C*M)*(M'*B*M)+0.5*(C*M)*(M'*B*M)/(M'*C*M)/(M'*C*M)-(C*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*C*M)+0.5*(B*M)*(M'*C*M)/(M'*B*M)/(M'*B*M)-(X2*X2')*M/(M'*B*M)+0.5*(B*M)*(M'*(X2*X2')*M)/(M'*B*M)/(M'*B*M)+0.5*(B*M)/(M'*B*M)/(M'*B*M)*(M'*(X2*X2')*M)-(X2*X2')*M/(M'*C*M)+0.5*(C*M)*(M'*(X2*X2')*M)/(M'*C*M)/(M'*C*M)+0.5*(C*M)/(M'*C*M)/(M'*C*M)*(M'*(X2*X2')*M)-(A*M)/(M'*C*M)+0.5*(C*M)/(M'*C*M)/(M'*C*M)*(M'*A*M)+0.5*(C*M)*(M'*A*M)/(M'*C*M)/(M'*C*M)-(C*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*C*M)+0.5*(A*M)*(M'*C*M)/(M'*A*M)/(M'*A*M)-(X3*X3')*M/(M'*A*M)+0.5*(A*M)*(M'*(X3*X3')*M)/(M'*A*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*(X3*X3')*M)-(X3*X3')*M/(M'*C*M)+0.5*(C*M)*(M'*(X3*X3')*M)/(M'*C*M)/(M'*C*M)+0.5*(C*M)/(M'*C*M)/(M'*C*M)*(M'*(X3*X3')*M);

problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X1*X1')*M)/(M'*A*M)-(X1*X1')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X1*X1')*M)/(M'*B*M)-(X1*X1')*M/(M'*B*M)-(B*M)/(M'*C*M)+(C*M)/(M'*C*M)*(M'*B*M)/(M'*C*M)-(C*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*C*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*(X2*X2')*M)/(M'*B*M)-(X2*X2')*M/(M'*B*M)+(C*M)/(M'*C*M)*(M'*(X2*X2')*M)/(M'*C*M)-(X2*X2')*M/(M'*C*M)-(A*M)/(M'*C*M)+(C*M)/(M'*C*M)*(M'*A*M)/(M'*C*M)-(C*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*C*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X3*X3')*M)/(M'*A*M)-(X3*X3')*M/(M'*A*M)+(C*M)/(M'*C*M)*(M'*(X3*X3')*M)/(M'*C*M)-(X3*X3')*M/(M'*C*M);



checkgradient(problem);
%checkdiff(problem);

%[x, xcost, info, options] = steepestdescent(problem);
[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

%xlswrite('1_2_3_dim3.xlsx',x);
%xlswrite('i_proj_3.xlsx',x);
