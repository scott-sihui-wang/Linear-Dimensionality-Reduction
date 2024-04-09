n=18;
d=1;

X=(x1-x2)';

manifold=grassmannfactory(n,d);
problem.M=manifold;

r=0.99;

C=r*B+(1-r)*A;

problem.cost = @(M) -logdet(M'*C*M)-(r-1)*logdet(M'*A*M)+r*logdet(M'*B*M)+r*(r-1)*(X'*M)/(M'*C*M)*(M'*X);
problem.egrad = @(M) -2*(C*M)/(M'*C*M)-2*(r-1)*(A*M)/(M'*A*M)+2*r*(B*M)/(M'*B*M)-2*r*(r-1)*(C*M)/(M'*C*M)*(M'*(X*X')*M)/(M'*C*M)+2*r*(r-1)*(X*X')*M/(M'*C*M);

checkgradient(problem);

options.maxiter=5000;
options.tolgradnorm=1e-9;

[x, xcost, info, options] = trustregions(problem,[],options);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');