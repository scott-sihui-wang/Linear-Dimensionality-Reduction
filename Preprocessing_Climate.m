clear all;clc;
data=load('.\climate.dat');
data=sortrows(data,21,'ascend');
data=data(:,3:20);
d_1=data(1:46,:);
d_2=data(47:540,:);
%d_2=d_2(randperm(494,46),:);
x1=mean(d_1);
x2=mean(d_2);
A=zeros(18);
B=zeros(18);
for i=1:46
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./46;
for i=1:494
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./494;
