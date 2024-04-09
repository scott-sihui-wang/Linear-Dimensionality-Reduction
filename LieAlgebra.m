n=100;
lieA=zeros(n,n*n*(n-1)/2);
for i=1:n-1
    for j=i+1:n
        lieA(i,((2*n-i)*(i-1)/2+j-i-1)*n+j)=1;
        lieA(j,((2*n-i)*(i-1)/2+j-i-1)*n+i)=-1;
    end
end

        