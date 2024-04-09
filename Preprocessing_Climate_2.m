data=load('.\climate.dat');
data=data(:,3:21);
avg=0;
for test=1:100
    r=randperm(size(data,1));
    data=data(r,:);
    train_size=round(length(data)*0.7);
    train_set=data(1:train_size,:);
    test_set=data(train_size+1:length(data),:);
    train_set=sortrows(train_set,19,'ascend');
    x1=zeros(1,18);
    x2=zeros(1,18);
    A=zeros(18);
    B=zeros(18);
    i=1;
    cnt1=0;cnt2=0;
    while train_set(i,19)==0
        x1=x1+train_set(i,1:18);
        cnt1=cnt1+1;
        i=i+1;
    end
    x1=x1./cnt1;
    while i<=length(train_set) && train_set(i,19)==1
        x2=x2+train_set(i,1:18);
        cnt2=cnt2+1;
        i=i+1;
    end
    x2=x2./cnt2;
    for i=1:length(train_set)
        if train_set(i,19)==0
            t=train_set(i,1:18)-x1;
            A=A+t'*t;
        elseif train_set(i,19)==1
            t=train_set(i,1:18)-x2;
            B=B+t'*t;
        end
    end
    A=A./cnt1;B=B./cnt2;

    n=18;
    d=1;

    X=(x1-x2)';

    manifold=grassmannfactory(n,d);
    problem.M=manifold;

    C=(A+B)/2;

    problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
    problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);
    
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
        d=test_set(i,1:18)*x;
        l1=-s1-(d-m_1)/c_1*(d-m_1)';
        l2=-s2-(d-m_2)/c_2*(d-m_2)';
        if test_set(i,19)==0
            if(l1<l2)
                err1=err1+1;
            end
        elseif test_set(i,19)==1
            if(l2<l1)
                err2=err2+1;
            end
        end
    end
    percentage=(err1+err2)/length(test_set);
    avg=avg+percentage;
end
avg=avg/test;