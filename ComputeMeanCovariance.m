images=loadMNISTImages('train-images.idx3-ubyte');
labels=loadMNISTLabels('train-labels.idx1-ubyte');
row=size(images,1);
col=size(images,2);
cnt=zeros(10,1);
means=zeros(10,row);
cov=zeros(10,row*row);
for i=1:col
    if (labels(i)==0)
        labels(i)=labels(i)+10;
    end
    cnt(labels(i))=cnt(labels(i))+1;
    means(labels(i),:)=means(labels(i),:)+(images(:,i))';
end
for i=1:10
    means(i,:)=means(i,:)/cnt(i);
end
for i=1:col
    MatCov=(images(:,i)-(means(labels(i),:))')*((images(:,i))'-means(labels(i),:));
    MatCov=reshape(MatCov,[1,row*row]);
    cov(labels(i),:)=cov(labels(i),:)+MatCov;
    if (mod(i,300)==0)
        fprintf("Computing covariance matrix, %6.2f percent completed...\n",i*1.0/col*100)
    end
end
for i=1:10
    cov(i,:)=cov(i,:)/cnt(i);
end
xlswrite('means.xlsx',means');
xlswrite('cov.xlsx',cov');