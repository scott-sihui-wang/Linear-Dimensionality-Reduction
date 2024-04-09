data=textread('.\Seed.txt');
data=data(:,1:7);
d_1=data(1:70,:);
d_2=data(71:140,:);
d_3=data(141:210,:);
x1=mean(d_1);
x2=mean(d_2);
x3=mean(d_3);
A=zeros(7);
B=zeros(7);
C=zeros(7);
for i=1:70
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./70;
for i=1:70
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./70;
for i=1:70
    t=d_3(i,:)-x3;
    C=C+t'*t;
end
C=C./70;