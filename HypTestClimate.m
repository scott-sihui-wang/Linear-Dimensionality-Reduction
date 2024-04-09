c_1=x'*A*x;
c_2=x'*B*x;
m_1=x1*x;
m_2=x2*x;
s1=logdet(c_1);
s2=logdet(c_2);
err1=0;err2=0;
for i=1:46
    d=d_1(i,:)*x;
    l1=-s1-(d-m_1)/c_1*(d-m_1)';
    l2=-s2-(d-m_2)/c_2*(d-m_2)';
    if(l1<l2)
        err1=err1+1;
    end
end
for i=1:494
    d=d_2(i,:)*x;
    l1=-s1-(d-m_1)/c_1*(d-m_1)';
    l2=-s2-(d-m_2)/c_2*(d-m_2)';
    if(l2<l1)
        err2=err2+1;
    end
end