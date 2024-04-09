methods=["KL","KLr","SKL","H"];
dim=[1,2,3];
l_m=length(methods);
l_d=length(dim);
results=zeros(l_m,l_d);
data=load('.\HTRU_2.csv');
for i=1:l_m
    for j=1:l_d
        for k=1:100
            r=randperm(size(data,1));
            data=data(r,:);
            train_size=round(length(data)*0.7);
            train_set=data(1:train_size,:);
            test_set=data(train_size+1:length(data),:);
            train_set=sortrows(train_set,9,'ascend');
            x1=zeros(1,8);
            x2=zeros(1,8);
            A=zeros(8);
            B=zeros(8);
            l=1;
            cnt1=0;cnt2=0;
            while train_set(l,9)==0
                x1=x1+train_set(l,1:8);
                cnt1=cnt1+1;
                l=l+1;
            end
            x1=x1./cnt1;
            while l<=length(train_set) && train_set(l,9)==1
                x2=x2+train_set(l,1:8);
                cnt2=cnt2+1;
                l=l+1;
            end
            x2=x2./cnt2;
            for l=1:length(train_set)
                if train_set(l,9)==0
                    t=train_set(l,1:8)-x1;
                    A=A+t'*t;
                elseif train_set(l,9)==1
                    t=train_set(l,1:8)-x2;
                    B=B+t'*t;
                end
            end
            A=A./cnt1;B=B./cnt2;

            n=8;
            d=dim(j);

            X=(x1-x2)';

            manifold=grassmannfactory(n,d);
            problem.M=manifold;

            C=(A+B)/2;
            
            if methods(i)=="KL"
                
                problem.cost = @(M) -0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.5*trace((M'*A*M)\(M'*B*M))-0.5*(X'*M)/(M'*A*M)*(M'*X);
                problem.egrad = @(M) -(A*M)/(M'*A*M)+(B*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M);

            elseif methods(i)=="KLr"
                
                problem.cost = @(M) -0.5*logdet(M'*B*M)+0.5*logdet(M'*A*M)-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*B*M)*(M'*X);
                problem.egrad = @(M) -(B*M)/(M'*B*M)+(A*M)/(M'*A*M)-(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);

            elseif methods(i)=="SKL"
                
                problem.cost = @(M) -0.5*trace((M'*A*M)\(M'*B*M))-0.5*trace((M'*B*M)\(M'*A*M))-0.5*(X'*M)/(M'*A*M)*(M'*X)-0.5*(X'*M)/(M'*B*M)*(M'*X);
                problem.egrad = @(M) -(A*M)/(M'*B*M)+(B*M)/(M'*B*M)*(M'*A*M)/(M'*B*M)-(B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*B*M)/(M'*A*M)+(A*M)/(M'*A*M)*(M'*(X*X')*M)/(M'*A*M)-(X*X')*M/(M'*A*M)+(B*M)/(M'*B*M)*(M'*(X*X')*M)/(M'*B*M)-(X*X')*M/(M'*B*M);
                
            elseif methods(i)=="H"
                
                problem.cost = @(M) -logdet(M'*C*M)+0.5*logdet(M'*A*M)+0.5*logdet(M'*B*M)-0.25*(X'*M)/(M'*C*M)*(M'*X);
                problem.egrad = @(M) -2*(C*M)/(M'*C*M)+(A*M)/(M'*A*M)+(B*M)/(M'*B*M)+0.5*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)-0.5*(X*X')*M/(M'*C*M);
                
            end

            options.maxiter=5000;

            [x, xcost, info, options] = trustregions(problem,[],options);

            c_1=x'*A*x;
            c_2=x'*B*x;
            m_1=x1*x;
            m_2=x2*x;
            s1=logdet(c_1);
            s2=logdet(c_2);
            err1=0;err2=0;
            for l=1:length(test_set)
                d=test_set(l,1:8)*x;
                l1=-s1-(d-m_1)/c_1*(d-m_1)';
                l2=-s2-(d-m_2)/c_2*(d-m_2)';
                if test_set(l,9)==0
                    if(l1<l2)
                        err1=err1+1;
                    end
                elseif test_set(l,9)==1
                    if(l2<l1)
                        err2=err2+1;
                    end
                end
            end
            percentage=(err1+err2)/length(test_set);
            results(i,j)=results(i,j)+percentage;
        end
        results(i,j)=results(i,j)/k;
    end
end