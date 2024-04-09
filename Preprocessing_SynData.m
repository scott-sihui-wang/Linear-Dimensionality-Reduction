data=xlsread('D:\SynData.xlsx');
data=sortrows(data,101,'ascend');
d_1=data(1:1000,1:100);
d_2=data(1001:2000,1:100);
x1=mean(d_1);
x2=mean(d_2);
A=zeros(100);
B=zeros(100);
for i=1:1000
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./1000;
for i=1:1000
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./1000;