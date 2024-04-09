clear all
clc

n=100;
r=10;

A=xlsread('D:\Sigma1.xlsx');
X=xlsread('Mean_diff_positive.xlsx');
%X=zeros(n,1);

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
selected=sortrows(choice(1:r,:),5,'descend');
candidate=sortrows(choice(r+1:n,:),2,'descend');
%selected_positive=sortrows(selected(find(selected(:,5)>=0),:),5,'descend');
%selected_negative=sortrows(selected(find(selected(:,5)<0),:),5,'ascend');
%candidate_positive=sortrows(candidate(find(candidate(:,2)>=0),:),2,'descend');
%candidate_negative=sortrows(candidate(find(candidate(:,2)<0),:),2,'ascend');

sum(selected(:,4))*0.5

x1=zeros(n,r);
for i=1:r
    x1(:,i)=eig_vec(:,selected(i,1));
end

sel_coeff=zeros(n,r);
for i=1:r
    sel_coeff(selected(i,1),i)=1;
end

rnk=zeros(r,1);
for i=1:r
    rnk(i,1)=i;
end
selected=[selected rnk];

figure
hold on

for i=1:length(candidate(:,1))
    l1=selected(1,3); l2=candidate(i,3); a1=selected(1,2); a2=candidate(i,2);
    fun=@(x)MyDivergence(x,l1,l2,a1,a2);
    px=0:0.01:1;
    py=fun(px);
    plot(px,py)
    %options=optimset('Algorithm','active-set','MaxIter',1e+01,'TolFun',1e-1,'TolCon',1e-1);
    [angle,fval]=fmincon(fun,0,[],[],[],[],0,1);
    selected(1,3)=l1*cos(angle)*cos(angle)+l2*sin(angle)*sin(angle);
    selected(1,2)=a1*cos(angle)+a2*sin(angle);
    selected(1,4)=-fval;
    selected(1,5)=selected(1,2)/selected(1,3);
    sel_coeff(:,selected(1,6))=sel_coeff(:,selected(1,6))*cos(angle);
    sel_coeff(candidate(i,1),selected(1,6))=sin(angle);
    j=2;
    comp=selected(1,:);
    while comp(5)<selected(j,5) && j<=r
        j=j+1;
    end
    if j>2
        selected(1:j-2,:)=selected(2:j-1,:);
        selected(j-1,:)=comp;
    end
end

%for i=1:length(candidate_negative(:,1))
%    l1=selected_negative(1,3); l2=candidate_negative(i,3); a1=selected_negative(1,2); a2=candidate_negative(i,2);
%    fun=@(x)MyDivergence(x,l1,l2,a1,a2);
%    [angle,fval]=fmincon(fun,0,[],[],[],[],0,1);
%    selected_negative(1,3)=l1*cos(angle)*cos(angle)+l2*sin(angle)*sin(angle);
%    selected_negative(1,2)=a1*cos(angle)+a2*sin(angle);
%    selected_negative(1,4)=-fval;
%    selected_negative(1,5)=selected_negative(1,2)*(1+selected_negative(1,3))/selected_negative(1,3);
%    selected_negative=sortrows(selected_negative,5,'descend');
%end

sum(selected(:,4))*0.5

x2=eig_vec*sel_coeff;

function cost = MyDivergence(theta,lambda_1,lambda_2,alpha_1,alpha_2)
    u1=cos(theta).*cos(theta).*lambda_1+sin(theta).*sin(theta).*lambda_2;
    u2=cos(theta).*alpha_1+sin(theta).*alpha_2;
    cost=-log(abs(u1))-1./u1-1./u1.*u2.*u2;
end