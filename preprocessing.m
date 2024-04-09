function [means,covs,prior,data,row_bound,column_bound] = preprocessing(filename,ext,sorting_column,sorting_method,num_of_class,row_bound,column_bound,dim,prior_enable,truncate_class)
    if(ext=="csv"||ext=="data"||ext=="txt"||ext=="dat")
        data=load(filename);
    else
        if (ext=="xls"||ext=="xlsx")
            data=xlsread(filename);
        end
    end
    if sorting_column > 0
        data=sortrows(data,sorting_column,sorting_method);
    end
    means=zeros(dim,1,num_of_class);
    covs=zeros(dim,dim,num_of_class);
    prior=zeros(num_of_class,1);
    if truncate_class ~=0
        if row_bound(1,2)-row_bound(1,1)>row_bound(2,2)-row_bound(2,1)
            truncate_class=2;
        else
            truncate_class=1;
        end
    end
    for i=1:num_of_class
        if truncate_class~=0 && truncate_class~=i
            shuffle=row_bound(i,1)+randperm(row_bound(i,2)-row_bound(i,1)+1)-1;
            data(row_bound(i,1):row_bound(i,2),column_bound(1):column_bound(2))=data(shuffle,column_bound(1):column_bound(2));
            row_bound(i,2)=row_bound(i,1)+row_bound(truncate_class,2)-row_bound(truncate_class,1);
        end
        means(:,:,i)=(mean(data(row_bound(i,1):row_bound(i,2),column_bound(1):column_bound(2))))';
        for j=row_bound(i,1):row_bound(i,2)
            t=data(j,column_bound(1):column_bound(2))-(means(:,:,i))';
            covs(:,:,i)=covs(:,:,i)+t'*t;
        end
        prior(i)=row_bound(i,2)-row_bound(i,1)+1;
        covs(:,:,i)=covs(:,:,i)/prior(i);
    end
    if prior_enable
        ttl=sum(prior);
        prior=prior/ttl;
    else
        prior=ones(num_of_class,1);
    end
end

