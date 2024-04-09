% fun1=@(x) abs(1/sqrt(2*pi)*exp(-1/2*x.*x)-1/sqrt(2*pi)/2.0000001*exp(-1/2/(2.0000001*2.0000001)*x.*x));
% fun2=@(x) abs(1/sqrt(2*pi)*exp(-1/2*x.*x)-1/sqrt(2*pi)/2*exp(-1/2/4*x.*x));
% a=integral(fun1,-Inf,Inf)
% b=integral(fun2,-Inf,Inf)
% (a-b)/0.0000001
% t=-10:0.01:10;
% y=1/sqrt(2*pi)*exp(-1/2*t.*t)-1/sqrt(2*pi)/2.0000001*exp(-1/2/(2.0000001*2.0000001)*t.*t);
% plot(t,y)
% g=@(x) 1/sqrt(2*pi)*exp(-1/2*x.*x)-1/sqrt(2*pi)/2.0000001*exp(-1/2/(2.0000001*2.0000001)*x.*x);
% ub=fzero(g,0)
% fun3=@(x) 1/sqrt(2*pi)/2.0000001*exp(-1/2/(2.0000001*2.0000001)*x.*x)-1/sqrt(2*pi)/2*exp(-1/2/4*x.*x);
% c=2*integral(fun3,ub,-ub)
fun =@(x) exp(-x.*x).*x;
integral(fun,-1,1)