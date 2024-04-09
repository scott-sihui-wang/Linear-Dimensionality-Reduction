n = 784;
%dim=[25 50 75];
dim=3;
m=xlsread('means.xlsx');
c=xlsread('cov.xlsx');
e=c';
A=reshape(e(1,:),[784,784]);
B=reshape(e(2,:),[784,784]);

%X1=randn(n,1);
%X2=randn(n,1);
%X=X1-X2;
x=m(1,:)-m(2,:);
X=x';
%X(10)=1.0;
I1=[];
I2=[];
RES=[];

for i=1:length(dim)
    
    d=dim(i);
    
    manifold=grassmannfactory(n,d);
    problem.M=manifold;
    
    problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
    problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);
    
    
    
    
    checkgradient(problem);
    %checkdiff(problem);
    
    %[x, xcost, info, options] = steepestdescent(problem);
    [x, xcost, info, options] = trustregions(problem);
    
    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of f');
    
    img1=reshape(m(1,:)*(x*x'),[28,28]);
    img2=reshape(m(2,:)*(x*x'),[28,28]);
    
    I1=[I1 img1];
    I2=[I2, img2];
    
    name1=sprintf('1_dim%d.bmp',d);
    name2=sprintf('2_dim%d.bmp',d);
    
    %imwrite(img1,name1);
    %imwrite(img2,name2);
    
    name3=sprintf('1_2_dim%d.xlsx',d);
    %xlswrite(name3,x);
    
    RES=[RES,x];
    
end