clear all;clc;

method.name="KL";
method.text="KL divergence";
method(2).name="rKL";
method(2).text="Backward KL divergence";
method(3).name="SKL";
method(3).text="Symmetric KL divergence";
method(4).name="H";
method(4).text="Hellinger Distance";
method(5).name="F";
method(5).text="Fisher Information Distance";
method(6).name="TV";
method(6).text="K-Generalized Total Variation";

optimizer.step=3000;
optimizer.tolerance=1e-9;
optimizer.method=1;

r_dim=3;

total_case=32400;

Err=zeros(2,total_case);
Ttl=zeros(2,total_case);
Mis=zeros(3,total_case);
PErr=zeros(1,total_case);
txt_file_name=string(missing);
txt_balance=string(missing);
txt_prior=string(missing);
txt_method=string(missing);
txt_dim=string(missing);

q=1;

for k=1:100
    filename="D:\\SynData_"+k+".xlsx";
    for l=1:9
        if l<=5
            row_bound=[(l-1)*200+1,1000;1001,2000];
        else
            row_bound=[1,1000;1001,l*200];
        end
        for m=0:1
            for i=1:6
                if i~=6
                    [means,covs,prior,data,row_bound,column_bound]=preprocessing(filename,"xlsx",101,"ascend",2,row_bound,[1,100],100,m,false);
                else
                    [means,covs,prior,data,row_bound,column_bound]=preprocessing(filename,"xlsx",101,"ascend",2,row_bound,[1,100],100,true,false);
                end
                proj=FID(means,covs,prior,method(i).name,100,r_dim,optimizer.step,optimizer.tolerance,optimizer.method);
                for j=1:r_dim
                    x=proj(:,1:j);
                    [r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,x]=criteria(means(:,1,1),means(:,1,2),covs(:,:,1),covs(:,:,2),x,0);
                    if i==6 && m==0
                        [err,ttl,mis_rate,p_err]=hypTest(data,row_bound,column_bound,x,r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,[1 1]);
                    else
                        [err,ttl,mis_rate,p_err]=hypTest(data,row_bound,column_bound,x,r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,prior);
                    end
                    Err(1,q)=err(1);Err(2,q)=err(2);
                    Ttl(1,q)=ttl(1);Ttl(2,q)=ttl(2);
                    Mis(1,q)=mis_rate(1);Mis(2,q)=mis_rate(2);Mis(3,q)=mis_rate(3);
                    PErr(q)=p_err;
                    txt_file_name(q)="SynData_"+k+".xlsx";
                    if l<=5
                        K=5/(6-l);
                    else
                        K=(l-5)/5;
                    end
                    txt_balance(q)="k="+K;
                    if m==0
                        txt_prior(q)="No Prior";
                    else
                        txt_prior(q)="With Prior";
                    end
                    txt_method(q)=method(i).text;
                    txt_dim(q)="Dimension="+j;
                    q
                    q=q+1;
                end
            end
        end
    end
end

Err=Err';
Ttl=Ttl';
Mis=Mis';
PErr=PErr';
txt_file_name=txt_file_name';
txt_balance=txt_balance';
txt_prior=txt_prior';
txt_method=txt_method';
txt_dim=txt_dim';
result=[Err Ttl Mis PErr txt_file_name txt_balance txt_prior txt_method txt_dim];
xlswrite("Test_Result_Syn.xlsx",result);