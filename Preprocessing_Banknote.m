data=load('.\banknote.txt');
avg=0;
for test=1:100
    r=randperm(size(data,1));
    data=data(r,:);
    train_size=round(length(data)*0.7);
    train_set=data(1:train_size,:);
    test_set=data(train_size+1:length(data),:);
    train_set=sortrows(train_set,5,'ascend');
    x1=zeros(1,4);
    x2=zeros(1,4);
    A=zeros(4);
    B=zeros(4);
    i=1;
    cnt1=0;cnt2=0;
    while train_set(i,5)==0
        x1=x1+train_set(i,1:4);
        cnt1=cnt1+1;
        i=i+1;
    end
    x1=x1./cnt1;
    while i<=length(train_set) && train_set(i,5)==1
        x2=x2+train_set(i,1:4);
        cnt2=cnt2+1;
        i=i+1;
    end
    x2=x2./cnt2;
    for i=1:length(train_set)
        if train_set(i,5)==0
            t=train_set(i,1:4)-x1;
            A=A+t'*t;
        elseif train_set(i,5)==1
            t=train_set(i,1:4)-x2;
            B=B+t'*t;
        end
    end
    A=A./cnt1;B=B./cnt2;

    n=4;
    d=1;

    X=(x1-x2)';

    manifold=grassmannfactory(n,d);
    problem.M=manifold;

    problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
    problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);

    options.maxiter=5000;

    [x, xcost, info, options] = trustregions(problem,[],options);

    c_1=x'*A*x;
    c_2=x'*B*x;
    m_1=x1*x;
    m_2=x2*x;
    s1=logdet(c_1);
    s2=logdet(c_2);
    err1=0;err2=0;
    for i=1:length(test_set)
        d=test_set(i,1:4)*x;
        l1=-s1-(d-m_1)/c_1*(d-m_1)';
        l2=-s2-(d-m_2)/c_2*(d-m_2)';
        if test_set(i,5)==0
            if(l1<l2)
                err1=err1+1;
            end
        elseif test_set(i,5)==1
            if(l2<l1)
                err2=err2+1;
            end
        end
    end
    percentage=(err1+err2)/length(test_set);
    avg=avg+percentage;
end
avg=avg/test;