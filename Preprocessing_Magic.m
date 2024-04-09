data=load('.\magic.data');
d_1=data(1:12332,1:10);
d_2=data(12333:19020,1:10);
x1=mean(d_1);
x2=mean(d_2);
A=zeros(10);
B=zeros(10);
for i=1:12332
    t=d_1(i,:)-x1;
    A=A+t'*t;
end
A=A./12332;
for i=1:6688
    t=d_2(i,:)-x2;
    B=B+t'*t;
end
B=B./6688;