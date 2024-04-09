function [err,ttl,mis_rate,p_err]=hypTest(data,row_bound,column_bound,proj,r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,prior)
    err=zeros(2,1);
    ttl=zeros(2,1);
    mis_rate=zeros(3,1);
    ttl(1)=row_bound(1,2)-row_bound(1,1)+1;
    ttl(2)=row_bound(2,2)-row_bound(2,1)+1;
    for i=1:2
        for j=row_bound(i,1):row_bound(i,2)
            sample=data(j,column_bound(1):column_bound(2));
            t=proj'*sample';
            prob1=2*log(prior(1))-const1-(t-r_mu1)'/r_cov1*(t-r_mu1);
            prob2=2*log(prior(2))-const2-(t-r_mu2)'/r_cov2*(t-r_mu2);
            if prob1<prob2 && i==1
                err(1)=err(1)+1;
            end
            if prob1>prob2 && i==2
                err(2)=err(2)+1;
            end
        end
    end
    mis_rate(1)=err(1)/ttl(1);
    mis_rate(2)=err(2)/ttl(2);
    mis_rate(3)=(err(1)+err(2))/(ttl(1)+ttl(2));
    p_err=mis_rate(1)+mis_rate(2);
end