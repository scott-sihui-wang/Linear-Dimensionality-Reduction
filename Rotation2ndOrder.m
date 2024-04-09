clear;clc;

n=100;
r=10;

fac_tab=zeros(1,1000);
fac_tab(1,1)=1;
for i=2:length(fac_tab)
    fac_tab(1,i)=i*fac_tab(1,i-1);
end

A=xlsread('D:\Sigma1.xlsx');
B=xlsread('D:\Sigma2.xlsx');
X=xlsread('Mean_diff_positive.xlsx');

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

f=cost(A,B,X,x1,r);
fprintf('original: cost=%6.5f\n',f);

stepconst=0.1;
steplimit=0.000000000000001; %1e-15
cnt=0;

for i=1:10000 %1e+10
    
    M=x1(:,1:r);
    AM=A*M;
    MAM=M'*AM;
    BM=B*M;
    MBM=M'*BM;
    XXM=(X*X')*M;
    MXXM=M'*XXM;
    G1=-AM/MBM+BM/MBM*MAM/MBM-BM/MAM+AM/MAM*MBM/MAM+AM/MAM*MXXM/MAM-XXM/MAM+BM/MBM*MXXM/MBM-XXM/MBM;
    G=zeros(n,n);
    G(:,1:r)=G1;
    la=G'*x1-x1'*G; % -grad
    
    dG=zeros(n,r,n,n);
    
%     st=0.000001;
%     delta_G=zeros(n,n,n,n);
%     
%     for j=1:n
%         for k=1:n
%             rot=eye(n);
%             if j~=k
%                 rot(j,j)=cos(st);
%                 rot(j,k)=sin(st);
%                 rot(k,j)=-sin(st);
%                 rot(k,k)=cos(st);
%             end
%             XX=x1*rot;
%             MM=XX(:,1:r);
%             AMM=A*MM;
%             MMAMM=MM'*AMM;
%             BMM=B*MM;
%             MMBMM=MM'*BMM;
%             XXMM=(X*X')*MM;
%             MMXXMM=MM'*XXMM;
%             GG=-AMM/MMBMM+BMM/MMBMM*MMAMM/MMBMM-BMM/MMAMM+AMM/MMAMM*MMBMM/MMAMM+AMM/MMAMM*MMXXMM/MMAMM-XXMM/MMAMM+BMM/MMBMM*MMXXMM/MMBMM-XXMM/MMBMM;
%             GGG=zeros(n,n);
%             GGG(:,1:r)=GG;
%             delta_G(:,:,j,k)=((GGG'*XX-XX'*GGG)-(G'*x1-x1'*G))*1/st;
%         end
%     end
%     
%     comp=zeros(n,n);
    
    for j=1:n
        for k=1:r
            E=zeros(n,r);
            E(j,k)=1;
            AE=A*E;
            BE=B*E;
            MAE=M'*AE;
            MBE=M'*BE;
            MAE=MAE+MAE';
            MBE=MBE+MBE';
            XXE=X*X'*E;
            MXXE=M'*XXE;
            MXXE=MXXE+MXXE';
            P1=AM/MBM*MBE/MBM-AE/MBM;
            P2=BM/MAM*MAE/MAM-BE/MAM;
            P3=BE/MBM*MAM/MBM-BM/MBM*MBE/MBM*MAM/MBM+BM/MBM*MAE/MBM-BM/MBM*MAM/MBM*MBE/MBM;
            P4=AE/MAM*MBM/MAM-AM/MAM*MAE/MAM*MBM/MAM+AM/MAM*MBE/MAM-AM/MAM*MBM/MAM*MAE/MAM;
            P5=AE/MAM*MXXM/MAM-AM/MAM*MAE/MAM*MXXM/MAM+AM/MAM*MXXE/MAM-AM/MAM*MXXM/MAM*MAE/MAM;
            P6=BE/MBM*MXXM/MBM-BM/MBM*MBE/MBM*MXXM/MBM+BM/MBM*MXXE/MBM-BM/MBM*MXXM/MBM*MBE/MBM;
            P7=XXM/MAM*MAE/MAM-XXE/MAM;
            P8=XXM/MBM*MBE/MBM-XXE/MBM;
            dG(:,:,j,k)=P1+P2+P3+P4+P5+P6+P7+P8;
        end
    end
    
    DG=zeros(n,r,n,n);
    DM=zeros(n,n,n,n);
    H=zeros(n,n,n,n);
    
    for j=1:n
        for k=1:n
            for l=1:n
                DG(:,:,j,k)=DG(:,:,j,k)+x1(l,j)*dG(:,:,l,k)-x1(l,k)*dG(:,:,l,j);
            end
            DM(:,j,j,k)=-x1(:,k);
            DM(:,k,j,k)=x1(:,j);
            H1=zeros(n,n);
            dH=G1'*DM(:,:,j,k)+(DG(:,:,j,k))'*x1;
            H1(1:r,:)=dH;
            H1=H1-H1';
            H(:,:,j,k)=H1;
%             comp(j,k)=max(max(abs(H(:,:,j,k)-delta_G(:,:,j,k))));
        end
    end
    
    HMat=zeros(r*(n-r),r*(n-r));
    row=1;
    for j=1:r
        for k=r+1:n
            h1=H(:,:,j,k);
            h2=h1(1:r,r+1:n);
            h3=reshape(h2',1,r*(n-r));
            HMat(row,:)=h3;
            row=row+1;
        end
    end
    
    g1=la(1:r,r+1:n);
    g2=reshape(g1',r*(n-r),1);
    
    dr=HMat\g2;
    dr2=-reshape(dr,n-r,r);
    dr3=-dr2';
    
    rot=zeros(n,n);
    rot(1:r,r+1:n)=dr3;
    rot(r+1:n,1:r)=dr2;
    
    step=-stepconst/max(max(rot));
    R=pwr(rot,step,n,fac_tab);
    
    %step=stepconst/max(max(la));
    %R=pwr(la,step,n,fac_tab);
    
    t=x1*R;

    f1=cost(A,B,X,t,r);
    fprintf('f1 cost=%6.8f\n',f1);
    if(f1<f)
        f=f1;
        x1=t;
        cnt=cnt+1;
        if(cnt>2)
            stepconst=stepconst*2;
            cnt=0;
        end
    elseif(f1>=f)
        cnt=0;
        if(stepconst>steplimit)
            stepconst=stepconst*0.5;
            f1=f;
        else
            fprintf('Minimum Reached.');
            break;
        end
    end
    fprintf('i=%6.0f; cost=%6.8f\n',i,f);
    
end

function R=pwr(X,d,n,tab)
    R=eye(n);
    X=d*X;
    for i=1:20
        R=R+X^i/tab(1,i);
    end
end

function y=cost(A,B,X,N,r)
    M=N(:,1:r);
    y=-0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
end