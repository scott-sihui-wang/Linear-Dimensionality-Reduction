function [r_mu1,r_mu2,r_cov1,r_cov2,const1,const2,proj]=criteria(mu1,mu2,cov1,cov2,proj,transformer)
    if transformer~=0
        proj=transformer*proj;
    end
    r_mu1=proj'*mu1;
    r_mu2=proj'*mu2;
    r_cov1=proj'*cov1*proj;
    r_cov2=proj'*cov2*proj;
    const1=logdet(r_cov1);
    const2=logdet(r_cov2);
end