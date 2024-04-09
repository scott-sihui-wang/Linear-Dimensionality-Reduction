n = 784;
d = 3;

m_tol=1e-6;

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

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) totalCost(COV,X,M,10);
problem.egrad = @(M) totalGrad(COV,X,M,10);

checkgradient(problem);
%checkdiff(problem);

%[x, xcost, info, options] = steepestdescent(problem);
[x, xcost, info, options] = trustregions(problem);

[c_info, r_info]=size(info);

while (info(r_info).gradnorm>m_tol && r_info>=options.maxiter)
    [x, xcost, info, options] = trustregions(problem,x);
    [c_info, r_info]=size(info);
end

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

%xlswrite('NEW_MNIST_dim3.xlsx',x);

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
    %g=-AM/MBM+0.5*BM/MBM/MBM*MAM+0.5*BM*MAM/MBM/MBM-BM/MAM+0.5*AM/MAM/MAM*MBM+0.5*AM*MBM/MAM/MAM-XXM/MAM+0.5*AM*MXXM/MAM/MAM+0.5*AM/MAM/MAM*MXXM-XXM/MBM+0.5*BM*MXXM/MBM/MBM+0.5*BM/MBM/MBM*MXXM;
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