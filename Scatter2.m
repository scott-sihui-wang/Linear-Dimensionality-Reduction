images=loadMNISTImages('train-images.idx3-ubyte');
labels=loadMNISTLabels('train-labels.idx1-ubyte');
col=size(images,2);
%proj=xlsread('MNIST_dim3.xlsx');
proj=x;

e1=[1 0 0];%1-red
e2=[0 1 0];%2-green
e3=[0 0 1];%3-blue
e4=[1 1 0];%4-yellow
e5=[1 0 1];%5-purple red
e6=[0 1 1];%6-green blue
e7=[1 0.5 0];%7-orange
e8=[0.5 0.25 0];%8-brown
e9=[0.5 0 1];%9-purple blue
e10=[0.25 0.25 0.25];%0-dark gray
E=[e1; e2; e3; e4; e5; e6; e7; e8; e9; e10];
pts=zeros(col,4);

for i=1:col
    res=proj'*images(:,i);
    pts(i,1)=res(1);
    pts(i,2)=res(2);
    pts(i,3)=res(3);
    if(labels(i)==0)
        pts(i,4)=10;
    else
        pts(i,4)=labels(i);
    end
end

xlswrite('MNIST_scatter_pts_dim3.xlsx',pts);

figure
hold on

%for i=1:10
%    clr=E(i,:)
%    scatter(i-5,0,10,'MarkerEdgeColor',clr,'MarkerFaceColor',0.8*clr,'LineWidth',1.0);
%end

for i=1:col
    index=pts(i,4);
    clr=E(index,:);
    scatter3(pts(i,1),pts(i,2),pts(i,3),5,'MarkerEdgeColor',clr,'MarkerFaceColor',0.8*clr,'LineWidth',1.0);
    if(mod(i,300)==0)
        fprintf('%6.2f percent completed...\n',i*1.0/col*100);
    end
end