n=100;
r=10;
A=xlsread('D:\Sigma1.xlsx');
X=xlsread('Mean_diff_positive.xlsx');
%M=randn(100,10);
%M=x1(:,1:10);
M=x;
c=Cost(A,X,M);
df=Diff_Formula(A,X,M);
%df=inverse(A,M,1,1,n,r);
%dx=NInv(A,M,1,1,1e-13,n,r);
dx1=zeros(100,10);
dx2=zeros(100,10);
dx3=zeros(100,10);
dx4=zeros(100,10);
dx5=zeros(100,10);
for i=1:100
    for j=1:10
        dx1(i,j)=Diff_Sim(A,X,M,i,j,0.01);
        dx2(i,j)=Diff_Sim(A,X,M,i,j,0.0001);
        dx3(i,j)=Diff_Sim(A,X,M,i,j,0.000001);
        dx4(i,j)=Diff_Sim(A,X,M,i,j,0.00000001);
        dx5(i,j)=Diff_Sim(A,X,M,i,j,0.0000000001);
    end
end

function x = inverse(A,M,i,j,n,r)
    E=zeros(r,n);
    E(i,j)=1;
    MAM=M'*A*M;
    x=-MAM\(MAM\(M'*A*E'+E*A*M));
end

function x=NInv(A,M,i,j,t,n,r)
    E=zeros(r,n);
    E(i,j)=t;
    x=1/t*(inv((M'+E)*A*(M+E'))-inv(M'*A*M));
end

function x = Cost(A,X,M)
    x=-0.5*logdet(M'*A*M)+0.5*logdet(M'*M)-0.5*trace((M'*A*M)\(M'*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
end
function dx = Diff_Formula(A,X,M)
    %dx=-(A*M)/(M'*A*M)+M/(M'*M)-M/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*M)+0.5*(A*M)*(M'*M)/(M'*A*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+0.5*(A*M)*(M'*(X*X')*M)/(M'*A*M)/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*(X*X')*M);
    dx=-(A*M)/(M'*A*M)+M/(M'*M)-M/(M'*A*M)+0.5*(A*M)/(M'*A*M)/(M'*A*M)*(M'*M)+0.5*(A*M)*(M'*M)/(M'*A*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M);
end
function dx = Diff_Sim(A,X,M,i,j,d)
    [r c]=size(M);
    Del=zeros(r,c);
    Del(i,j)=d;
    dx=(Cost(A,X,M+Del)-Cost(A,X,M))/d;
end