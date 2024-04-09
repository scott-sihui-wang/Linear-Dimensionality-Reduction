data=load('.\abalone.data');
data=sortrows(data,1,'ascend');
d_1=data(1:1528,2:9);
d_2=data(1529:2835,2:9);
d_3=data(2836:4177,2:9);
x1=mean(d_1);
x2=mean(d_2);
x3=mean(d_3);
A=zeros(8);
B=zeros(8);
C=zeros(8);
for i=1:1528
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./1528;
for i=1:1307
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./1307;
for i=1:1342
    t=d_3(i,:)-x3;
    C=C+t'*t;
end
C=C./1342;