clear;clc;

n=100;
r=10;

fac_tab=zeros(1,1000);
fac_tab(1,1)=1;
for i=2:length(fac_tab)
    fac_tab(1,i)=i*fac_tab(1,i-1);
end

A=xlsread('D:\Sigma1.xlsx');
%X=xlsread('mean_diff.xlsx');
X=xlsread('Mean_diff_positive.xlsx');
%X=zeros(n,1);
X=10*X;

[eig_vec,eig_val]=eig(A);
coeff=eig_vec\X;

choice=zeros(n,5);

for i=1:n
    choice(i,1)=i;
    choice(i,2)=coeff(i,1);
    choice(i,3)=eig_val(i,i);
    choice(i,4)=log(choice(i,3))+1/choice(i,3)+1/choice(i,3)*choice(i,2)*choice(i,2);
    choice(i,5)=choice(i,2)/choice(i,3);
end

choice=sortrows(choice,4,'descend');

x1=zeros(n,n);
orig=zeros(n,n);
for i=1:n
    x1(:,i)=eig_vec(:,choice(i,1));
    orig(:,i)=eig_vec(:,choice(i,1));
end

%swap=x1(:,1);
%x1(:,1)=x1(:,n);
%x1(:,n)=swap;
%swap=x1(:,2);
%x1(:,2)=x1(:,n-1);
%x1(:,n-1)=swap;

f=cost(A,X,x1,r);
fprintf('original: cost=%6.5f\n',f);

stepconst=0.1;
steplimit=0.000000000000001; %1e-15
cnt=0;

for i=1:10000000000 %1e+10
    
    M=x1(:,1:r);
    AM=A*M;
    MAM=M'*AM;
    MM=M'*M;
    XXM=(X*X')*M;
    MXXM=M'*XXM;
    G1=-AM/MAM+M/MM-M/MAM+AM/MAM*MM/MAM+AM/MAM*MXXM/MAM-XXM/MAM;
    %G1=AM/MAM*MXXM/MAM-XXM/MAM;
    %G1=-AM/MAM+M/MM-M/MAM+AM/MAM*MM/MAM;
    G=zeros(n,n);
    G(:,1:r)=G1;
    la=G'*x1-x1'*G; % -grad 
    
    %lie_alg=zeros(n,n);
    %t=min(0.000001,stepconst);
    %for j=1:r
    %    for k=r+1:n
    %        rot=eye(n);
    %        rot(j,j)=cos(t);
    %        rot(j,k)=sin(t);
    %        rot(k,j)=-sin(t);
    %        rot(k,k)=cos(t);
    %        diff=cost(A,X,x1*rot,r)-cost(A,X,x1,r);
    %        %fprintf('direction: %2.0f %2.0f: grad: %3.5f\n',j,k,diff);
    %        lie_alg(j,k)=diff/t;
    %        lie_alg(k,j)=-diff/t;
    %    end
    %end
    
    %la=x1*G'-G*x1';
    %la=G*x1'-x1*G';
    %la=zeros(n,n);
    %la(1:r,r+1:n)=rot_grad(1:r,r+1:n);
    %la(r+1:n,1:r)=rot_grad(r+1:n,1:r);
    
    %stepconst
    
    %step=1;
    
    step=stepconst/max(max(la));
    R=pwr(la,step,n,fac_tab);
    %R2=pwr(la,-step,n,fac_tab);
    
    t1=x1*R;
    %t2=x1*R2;
    f1=cost(A,X,t1,r);
    %f2=cost(A,X,t2,r);
    %if(f1>f && f2<f)
    %    fprintf("-\n");
    %    break;
    %    f=f2;
    %    x1=t2;
    %    cnt=cnt+1;
    %    if(cnt>2)
    %        stepconst=stepconst*2;
    %        cnt=0;
    %    end
    if(f1<f)
        %fprintf("+\n");
        f=f1;
        x1=t1;
        cnt=cnt+1;
        if(cnt>2)
            stepconst=stepconst*2;
            cnt=0;
            %fprintf('s+ ');
        end
    elseif(f1>=f)
        cnt=0;
        if(stepconst>steplimit)
            stepconst=stepconst*0.5;
            f1=f;
            %fprintf('s- ');
        else
            fprintf('Minimum Reached.');
            break;
        end
    %elseif (f1<f && f2<f)
    %    cnt=0;
    %    fprintf('Maximum Reached.');
    %    if(f1<f2)
    %        f=f1;
    %        x1=t1;
    %    else
    %        f=f2;
    %        x1=t2;
    %    end
    %else
    %    cnt=0;
    %    fprintf('Unidentified Cases');
    end
    fprintf('i=%6.0f; cost=%6.8f\n',i,f);
    
end

function R=pwr(X,d,n,tab)
    R=eye(n);
    X=d*X;
    for i=1:20
        R=R+X^i/tab(1,i);
        %R=R+X^i/factorial(i);
    end
end

function y=cost(A,X,N,r)
    M=N(:,1:r);
    y=-0.5*logdet(M'*A*M)+0.5*logdet(M'*M)-0.5*trace((M'*A*M)\(M'*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
    %y=-0.5*(X'*M)/(M'*A*M)*(M'*X);
    %y=-0.5*logdet(M'*A*M)+0.5*logdet(M'*M)-0.5*trace((M'*A*M)\(M'*M));
end