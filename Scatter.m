images=loadMNISTImages('train-images.idx3-ubyte');
labels=loadMNISTLabels('train-labels.idx1-ubyte');
col=size(images,2);
%proj=xlsread('1_2_3_dim3.xlsx');
%proj=xlsread('2V3_refined_proj.xlsx');
proj=xlsread('2V3_new1.xlsx');
%A=[];
%B=[];

%proj=xlsread('i_proj_3.xlsx');
%in=xlsread('iris.xlsx');
%in=in(:,2:5);
%sample1=in(1:50,:);
%sample2=in(51:100,:);
%sample3=in(101:150,:);

%figure
%hold on

%for i=1:50
%    res1=sample1(i,:)*x;
%    res2=sample2(i,:)*x;
%    res3=sample3(i,:)*x;
%    scatter3(res1(1),res1(2),res1(3),10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.8 0 0],'LineWidth',1.5);
%    scatter3(res2(1),res2(2),res2(3),10,'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 0.8 0],'LineWidth',1.5);
%    scatter3(res3(1),res3(2),res3(3),10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 0.8],'LineWidth',1.5);
%end

%A=[];
%B=[];
%C=[];

%for i=1:50
%    res1=sample1(i,:)*x;
%    res2=sample2(i,:)*x;
%    res3=sample3(i,:)*x;
%    A=[A res1];
%    B=[B res2];
%    C=[C res3];
%end

%h1=histogram(A)
%hold on
%h2=histogram(B)
%h3=histogram(C)

figure
hold on

for i=1:col
    if (labels(i)==2)
        res1=proj'*images(:,i);
        scatter3(res1(1),res1(2),res1(3),3,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.8 0 0],'LineWidth',0.2);
    end
    if (labels(i)==3)
        res2=proj'*images(:,i);
        scatter3(res2(1),res2(2),res2(3),3,'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 0.8 0],'LineWidth',0.2);
    end
    %if (labels(i)==3)
    %    res3=proj'*images(:,i);
    %    scatter3(res3(1),res3(2),res3(3),3,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 0.8],'LineWidth',0.2);
    %end
    if(mod(i,300)==0)
        fprintf("%6.2f percent completed...\n",i*1.0/col*100);
    end
end
%h1=histogram(A)