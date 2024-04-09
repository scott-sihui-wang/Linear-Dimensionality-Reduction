n = 784;
dim=[3];
m=xlsread('means.xlsx');
c=xlsread('cov.xlsx');
e=c';

A=reshape(e(1,:),[784,784]);
B=reshape(e(2,:),[784,784]);
C=reshape(e(3,:),[784,784]);
D=reshape(e(4,:),[784,784]);
E=reshape(e(5,:),[784,784]);
F=reshape(e(6,:),[784,784]);
G=reshape(e(7,:),[784,784]);
H=reshape(e(8,:),[784,784]);
I=reshape(e(9,:),[784,784]);
J=reshape(e(10,:),[784,784]);
COV=[A B C D E F G H I J];

x1=(m(1,:))';
x2=(m(2,:))';
x3=(m(3,:))';
x4=(m(4,:))';
x5=(m(5,:))';
x6=(m(6,:))';
x7=(m(7,:))';
x8=(m(8,:))';
x9=(m(9,:))';
x10=(m(10,:))';
X=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10];

for i=1:10
    fn=sprintf('MNIST_Mean_%d.bmp',i);
    %imwrite(reshape(m(i,:),[28,28]),fn);
end

for i=1:length(dim)
    
    d=dim(i);
    
    manifold=grassmannfactory(n,d);
    %manifold=stiefelfactory(n,d);
    problem.M=manifold;
    
    problem.cost = @(M) totalCost(COV,X,M,10);
    problem.egrad = @(M) totalGrad(COV,X,M,10);
    
    checkgradient(problem);
    %checkdiff(problem);
    
    options.maxiter=200000;
    options.verbosity=0;
    options.stopfun= @ mystopfun;
    %[x, xcost, info, options] = steepestdescent(problem);
    [x, xcost, info, options] = trustregions(problem);
    
    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of f');
    
    img1=reshape(m(1,:)*(x*x'),[28,28]);
    img2=reshape(m(2,:)*(x*x'),[28,28]);
    img3=reshape(m(3,:)*(x*x'),[28,28]);
    img4=reshape(m(4,:)*(x*x'),[28,28]);
    img5=reshape(m(5,:)*(x*x'),[28,28]);
    img6=reshape(m(6,:)*(x*x'),[28,28]);
    img7=reshape(m(7,:)*(x*x'),[28,28]);
    img8=reshape(m(8,:)*(x*x'),[28,28]);
    img9=reshape(m(9,:)*(x*x'),[28,28]);
    img10=reshape(m(10,:)*(x*x'),[28,28]);
    
    name1=sprintf('MNIST1_dim%d.bmp',d);
    name2=sprintf('MNIST2_dim%d.bmp',d);
    name3=sprintf('MNIST3_dim%d.bmp',d);
    name4=sprintf('MNIST4_dim%d.bmp',d);
    name5=sprintf('MNIST5_dim%d.bmp',d);
    name6=sprintf('MNIST6_dim%d.bmp',d);
    name7=sprintf('MNIST7_dim%d.bmp',d);
    name8=sprintf('MNIST8_dim%d.bmp',d);
    name9=sprintf('MNIST9_dim%d.bmp',d);
    name10=sprintf('MNIST0_dim%d.bmp',d);
    
    %imwrite(img1,name1);
    %imwrite(img2,name2);
    %imwrite(img3,name3);
    %imwrite(img4,name4);
    %imwrite(img5,name5);
    %imwrite(img6,name6);
    %imwrite(img7,name7);
    %imwrite(img8,name8);
    %imwrite(img9,name9);
    %imwrite(img10,name10);
    
    namep=sprintf('MNIST_dim%d.xlsx',d);
    %xlswrite(namep,x);
    
end

function[f]= myCost(A,B,X,M)
    MAM=M'*A*M;
    MBM=M'*B*M;
    XM=X'*M;
    f=-0.5*trace(MAM\MBM)-0.5*trace(MBM\MAM)-0.5*XM/MAM*XM'-0.5*XM/MBM*XM';
end

function[g]= egrad(A,B,X,M)
    AM=A*M;
    BM=B*M;
    MAM=M'*A*M;
    MBM=M'*B*M;
    XXM=(X*X')*M;
    MXXM=M'*XXM;
    g=-AM/MBM+BM/MBM*MAM/MBM-BM/MAM+AM/MAM*MBM/MAM+AM/MAM*MXXM/MAM-XXM/MAM+BM/MBM*MXXM/MBM-XXM/MBM;
end

function[F]= totalCost(COV,X,M,n)
    F=0;
    for i=1:n
        for j=(i+1):n
            F=F+myCost(COV(:,(i-1)*784+1:i*784),COV(:,(j-1)*784+1:j*784),X(:,i)-X(:,j),M);
        end
    end     
end

function[G]= totalGrad(COV,X,M,n)
    G=0;
    for i=1:n
        for j=(i+1):n
            G=G+egrad(COV(:,(i-1)*784+1:i*784),COV(:,(j-1)*784+1:j*784),X(:,i)-X(:,j),M);
        end
    end
end

function stopnow=mystopfun(problem, x,info,last)
    stopnow=(last>=3 && info(last-2).cost-info(last).cost<1e-3);
end