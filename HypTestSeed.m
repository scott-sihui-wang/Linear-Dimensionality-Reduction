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
for i=1:70
    d=d_1(i,:)*x;
    l1=-s1-(d-m_1)/c_1*(d-m_1)';
    l2=-s2-(d-m_2)/c_2*(d-m_2)';
    l3=-s3-(d-m_3)/c_3*(d-m_3)';
    if(l1<l2 || l1<l3)
        err1=err1+1;
    end
end
for i=1:70
    d=d_2(i,:)*x;
    l1=-s1-(d-m_1)/c_1*(d-m_1)';
    l2=-s2-(d-m_2)/c_2*(d-m_2)';
    l3=-s3-(d-m_3)/c_3*(d-m_3)';
    if(l2<l1 || l2<l3)
        err2=err2+1;
    end
end
for i=1:70
    d=d_3(i,:)*x;
    l1=-s1-(d-m_1)/c_1*(d-m_1)';
    l2=-s2-(d-m_2)/c_2*(d-m_2)';
    l3=-s3-(d-m_3)/c_3*(d-m_3)';
    if(l3<l1 || l3<l2)
        err3=err3+1;
    end
end