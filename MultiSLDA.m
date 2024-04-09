n = 784;
d = 3;

m=xlsread('means.xlsx');
c=xlsread('cov.xlsx');
e=c';
A=reshape(e(1,:),[784,784]);
B=reshape(e(2,:),[784,784]);
C=reshape(e(3,:),[784,784]);

x1=m(1,:)-m(2,:);
X1=x1';
x2=m(2,:)-m(3,:);
X2=x2';
x3=m(1,:)-m(3,:);
X3=x3';

manifold=grassmannfactory(n,d);
problem.M=manifold;

problem.cost = @(M) myCost(A,B,X1,M)+myCost(B,C,X2,M)+myCost(A,C,X3,M);
problem.egrad = @(M) egrad(A,B,X1,M)+egrad(B,C,X2,M)+egrad(A,C,X3,M);

checkgradient(problem);
%checkdiff(problem);

%[x, xcost, info, options] = steepestdescent(problem);
[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

%xlswrite('NEW_1_2_3_dim3.xlsx',x);
%xlswrite('i_proj_3.xlsx',x);

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
