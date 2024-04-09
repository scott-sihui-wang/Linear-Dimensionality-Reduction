in=xlsread('iris.xlsx');
in=in(:,2:5);
sample1=in(1:50,:);
sample2=in(51:100,:);
sample3=in(101:150,:);
mean1=[mean(sample1(:,1)) mean(sample1(:,2)) mean(sample1(:,3)) mean(sample1(:,4))];
mean2=[mean(sample2(:,1)) mean(sample2(:,2)) mean(sample2(:,3)) mean(sample2(:,4))];
mean3=[mean(sample3(:,1)) mean(sample3(:,2)) mean(sample3(:,3)) mean(sample3(:,4))];
cov1=zeros(4);
cov2=zeros(4);
cov3=zeros(4);
for i=1:50
    cov1=cov1+(sample1(i,:)-mean1)'*(sample1(i,:)-mean1);
    cov2=cov2+(sample2(i,:)-mean2)'*(sample2(i,:)-mean2);
    cov3=cov3+(sample3(i,:)-mean3)'*(sample3(i,:)-mean3);
end
cov1=cov1/50;
cov2=cov2/50;
cov3=cov3/50;
xlswrite('imean1.xlsx',mean1);
xlswrite('imean2.xlsx',mean2);
xlswrite('imean3.xlsx',mean3);
xlswrite('icov1.xlsx',cov1);
xlswrite('icov2.xlsx',cov2);
xlswrite('icov3.xlsx',cov3);
