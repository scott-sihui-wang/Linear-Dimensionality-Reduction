data=load('.\HTRU_2.csv');
data=sortrows(data,9,'ascend');
d_1=data(1:16259,1:8);
d_1=d_1(randperm(16259,1639),:);
d_2=data(16260:17898,1:8);
x1=mean(d_1);
x2=mean(d_2);
A=zeros(8);
B=zeros(8);
for i=1:1639
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./1639;
for i=1:1639
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./1639;