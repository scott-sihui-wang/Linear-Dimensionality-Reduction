data=xlsread('D:\ldr\OldResults\iris.xlsx');
label=ones(50,1);
data(1:50,1)=label;
data(51:100,1)=2*label;
data(101:150,1)=3*label;
avg=0;
for test=1:100
    r=randperm(size(data,1));
    data=data(r,:);
    train_size=round(length(data)*0.7);
    train_set=data(1:train_size,:);
    test_set=data(train_size+1:length(data),:);
    train_set=sortrows(train_set,1,'ascend');
    x1=zeros(1,4);
    x2=zeros(1,4);
    x3=zeros(1,4);
    A=zeros(4);
    B=zeros(4);
    C=zeros(4);
    i=1;
    cnt1=0;cnt2=0;cnt3=0;
    while train_set(i,1)==1
        x1=x1+train_set(i,2:5);
        cnt1=cnt1+1;
        i=i+1;
    end
    x1=x1./cnt1;
    while train_set(i,1)==2
        x2=x2+train_set(i,2:5);
        cnt2=cnt2+1;
        i=i+1;
    end
    x2=x2./cnt2;
    while i<=length(train_set) && train_set(i,1)==3
        x3=x3+train_set(i,2:5);
        cnt3=cnt3+1;
        i=i+1;
    end
    x3=x3./cnt3;
    for i=1:length(train_set)
        if train_set(i,1)==1
            t=train_set(i,2:5)-x1;
            A=A+t'*t;
        elseif train_set(i,1)==2
            t=train_set(i,2:5)-x2;
            B=B+t'*t;
        elseif train_set(i,1)==3
            t=train_set(i,2:5)-x3;
            C=C+t'*t;
        end
    end
    A=A./cnt1;B=B./cnt2;C=C./cnt3;

    n=4;
    d=1;

    X1=(x1-x2)';
    X2=(x2-x3)';
    X3=(x1-x3)';

    manifold=grassmannfactory(n,d);
    problem.M=manifold;

    problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X1'*M)/(M'*A*M)*(M'*X1)-0.5*(X1'*M)/(M'*B*M)*(M'*X1)-0.5*trace((M'*B*M)\(M'*C*M))-0.5*trace((M'*C*M)\(M'*B*M))-0.5*(X2'*M)/(M'*B*M)*(M'*X2)-0.5*(X2'*M)/(M'*C*M)*(M'*X2)-0.5*trace((M'*A*M)\(M'*C*M))-0.5*trace((M'*C*M)\(M'*A*M))-0.5*(X3'*M)/(M'*A*M)*(M'*X3)-0.5*(X3'*M)/(M'*C*M)*(M'*X3);
    problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X1*X1')*M)/(M'*A*M)-(X1*X1')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X1*X1')*M)/(M'*B*M)-(X1*X1')*M/(M'*B*M)-(B*M)/(M'*C*M)+(C*M)/(M'*C*M)*(M'*B*M)/(M'*C*M)-(C*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*C*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*(X2*X2')*M)/(M'*B*M)-(X2*X2')*M/(M'*B*M)+(C*M)/(M'*C*M)*(M'*(X2*X2')*M)/(M'*C*M)-(X2*X2')*M/(M'*C*M)-(A*M)/(M'*C*M)+(C*M)/(M'*C*M)*(M'*A*M)/(M'*C*M)-(C*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*C*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X3*X3')*M)/(M'*A*M)-(X3*X3')*M/(M'*A*M)+(C*M)/(M'*C*M)*(M'*(X3*X3')*M)/(M'*C*M)-(X3*X3')*M/(M'*C*M);

    options.maxiter=5000;

    [x, xcost, info, options] = trustregions(problem,[],options);

    c_1=x'*A*x;
    c_2=x'*B*x;
    c_3=x'*C*x;
    m_1=x1*x;
    m_2=x2*x;
    m_3=x3*x;
    s1=logdet(c_1);
    s2=logdet(c_2);
    s3=logdet(c_3);
    err1=0;err2=0;err3=0;
    for i=1:length(test_set)
        d=test_set(i,2:5)*x;
        l1=-s1-(d-m_1)/c_1*(d-m_1)';
        l2=-s2-(d-m_2)/c_2*(d-m_2)';
        l3=-s3-(d-m_3)/c_3*(d-m_3)';
        if test_set(i,1)==1
            if(l1<l2 || l1<l3)
                err1=err1+1;
            end
        elseif test_set(i,1)==2
            if(l2<l1 || l2<l3)
                err2=err2+1;
            end
        elseif test_set(i,1)==3
            if(l3<l1 || l3<l2)
                err3=err3+1;
            end
        end
    end
    percentage=(err1+err2+err3)/length(test_set);
    avg=avg+percentage;
end
avg=avg/test;