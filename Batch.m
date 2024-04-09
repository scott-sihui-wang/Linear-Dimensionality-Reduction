clear all;clc;

file.filename='.\abalone.data';
file.ext='data';
file.sorting_column=1;
file.sorting_method='ascend';
file.num_of_class=2;
file.row_bound=[1,1528;1529,2835];
file.column_bound=[2,9];
file.dim=8;
file.text="Dataset: Abalone Class 1 Vs Class 2";
file(2).filename='.\abalone.data';
file(2).ext='data';
file(2).sorting_column=1;
file(2).sorting_method='ascend';
file(2).num_of_class=2;
file(2).row_bound=[1,1528;2836,4177];
file(2).column_bound=[2,9];
file(2).dim=8;
file(2).text="Dataset: Abalone Class 1 Vs Class 3";
file(3).filename='.\abalone.data';
file(3).ext='data';
file(3).sorting_column=1;
file(3).sorting_method='ascend';
file(3).num_of_class=2;
file(3).row_bound=[1529,2835;2836,4177];
file(3).column_bound=[2,9];
file(3).dim=8;
file(3).text="Dataset: Abalone Class 2 Vs Class 3";
file(4).filename='.\banknote.txt';
file(4).ext='txt';
file(4).sorting_column=5;
file(4).sorting_method='ascend';
file(4).num_of_class=2;
file(4).row_bound=[1,762;763,1372];
file(4).column_bound=[1,4];
file(4).dim=4;
file(4).text="Dataset: Banknote";
file(5).filename='.\climate.dat';
file(5).ext='dat';
file(5).sorting_column=21;
file(5).sorting_method='ascend';
file(5).num_of_class=2;
file(5).row_bound=[1,46;47,540];
file(5).column_bound=[3,20];
file(5).dim=18;
file(5).text="Dataset: Climate";
file(6).filename='.\HTRU_2.csv';
file(6).ext='csv';
file(6).sorting_column=9;
file(6).sorting_method='ascend';
file(6).num_of_class=2;
file(6).row_bound=[1,16259;16260,17898];
file(6).column_bound=[1,8];
file(6).dim=8;
file(6).text="Dataset: HTRU_2";
file(7).filename='D:\ldr\OldResults\iris.xlsx';
file(7).ext='xlsx';
file(7).sorting_column=0;
file(7).sorting_method='ascend';
file(7).num_of_class=2;
file(7).row_bound=[1,50;51,100];
file(7).column_bound=[2,5];
file(7).dim=4;
file(7).text="Dataset: Iris Class 1 VS Class 2";
file(8).filename='D:\ldr\OldResults\iris.xlsx';
file(8).ext='xlsx';
file(8).sorting_column=0;
file(8).sorting_method='ascend';
file(8).num_of_class=2;
file(8).row_bound=[1,50;101,150];
file(8).column_bound=[2,5];
file(8).dim=4;
file(8).text="Dataset: Iris Class 1 VS Class 3";
file(9).filename='D:\ldr\OldResults\iris.xlsx';
file(9).ext='xlsx';
file(9).sorting_column=0;
file(9).sorting_method='ascend';
file(9).num_of_class=2;
file(9).row_bound=[51,100;101,150];
file(9).column_bound=[2,5];
file(9).dim=4;
file(9).text="Dataset: Iris Class 2 VS Class 3";
file(10).filename='.\magic.data';
file(10).ext='data';
file(10).sorting_column=0;
file(10).sorting_method='ascend';
file(10).num_of_class=2;
file(10).row_bound=[1,12332;12333,19020];
file(10).column_bound=[1,10];
file(10).dim=10;
file(10).text="Dataset: Magic";
file(11).filename='.\Seed.txt';
file(11).ext='txt';
file(11).sorting_column=0;
file(11).sorting_method='ascend';
file(11).num_of_class=2;
file(11).row_bound=[1,70;71,140];
file(11).column_bound=[1,7];
file(11).dim=7;
file(11).text="Dataset: Seed Class 1 Vs Class 2";
file(12).filename='.\Seed.txt';
file(12).ext='txt';
file(12).sorting_column=0;
file(12).sorting_method='ascend';
file(12).num_of_class=2;
file(12).row_bound=[1,70;141,210];
file(12).column_bound=[1,7];
file(12).dim=7;
file(12).text="Dataset: Seed Class 1 Vs Class 3";
file(13).filename='.\Seed.txt';
file(13).ext='txt';
file(13).sorting_column=0;
file(13).sorting_method='ascend';
file(13).num_of_class=2;
file(13).row_bound=[71,140;141,210];
file(13).column_bound=[1,7];
file(13).dim=7;
file(13).text="Dataset: Seed Class 2 Vs Class 3";

balance_set.prior=false;
balance_set.truncate=false;
balance_set.text="No Prior, No Truncate";
balance_set(2).prior=true;
balance_set(2).truncate=false;
balance_set(2).text="With Prior, No Truncate";
balance_set(3).prior=false;
balance_set(3).truncate=true;
balance_set(3).text="No Prior, With Truncate";

method.name="KL";
method.text="KL divergence";
method(2).name="rKL";
method(2).text="Backward KL divergence";
method(3).name="SKL";
method(3).text="Symmetric KL divergence";
method(4).name="H";
method(4).text="Hellinger Distance";

dim.val=1;
dim.text="Dimension=1";
dim(2).val=2;
dim(2).text="Dimension=2";
dim(3).val=3;
dim(3).text="Dimension=3";

transform.val=false;
transform.text="No Transformation";
transform(2).val=true;
transform(2).text="With Transformation";

optimizer.step=3000;
optimizer.tolerance=1e-9;
optimizer.method=2;

Err=zeros(2,936);
Ttl=zeros(2,936);
Mis=zeros(2,936);
PErr=zeros(1,936);
Txt=string(missing);

q=1;

for i=1:13 %dataset No
    for j=1:3 % prior and truncate
        for k=1:4 % divergences
            for l=1:3 % dimension
                for p=1:2 % transformation
                    [err,ttl,mis_rate,p_err]=process(file(i).filename,file(i).ext,file(i).sorting_column,file(i).sorting_method,file(i).num_of_class,file(i).row_bound,file(i).column_bound,balance_set(j).prior,balance_set(j).truncate,file(i).dim,dim(l).val,method(k).name,transform(p).val,optimizer.step,optimizer.tolerance,optimizer.method)
                    Err(1,q)=err(1);Err(2,q)=err(2);
                    Ttl(1,q)=ttl(1);Ttl(2,q)=ttl(2);
                    Mis(1,q)=mis_rate(1);Mis(2,q)=mis_rate(2);
                    PErr(q)=p_err;
                    Txt(q)=file(i).text+" "+balance_set(j).text+" "+method(k).text+" "+dim(l).text+" "+transform(p).text;
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
Txt=Txt';
result=[Err Ttl Mis PErr Txt];
xlswrite("Test_Result.xlsx",result);