function [u,xaxis,taxis,err] = Ch4_fd1dheat(theta,a,T,X,f,u0,M,N,u_precise,plt)
%solve u_t - a u_xx = f, t=0: u=u0
%with finite difference \theta formula

if nargin < 1
theta=1/2;
end
if nargin < 2
a=1;
end
if nargin < 3
T=.5;
end
if nargin < 4
X=1;
end
if nargin < 5
f=@(t,x)t-t;
end
if nargin < 6
u0=@(x)sin(pi*x);
end
if nargin < 7
M=10;
end
if nargin < 8
N=10;
end
if nargin < 9
u_precise=@(t,x)sin(pi*x).*exp(-pi^2*t);
end
if nargin < 10
plt=1;
end

tau=T/M;
h=X/N;
u=zeros(N-1,M);
xaxis=h:h:X-h;
taxis=0:tau:T;
[tgrid xgrid]=meshgrid(taxis,xaxis);
u(:,1)=u0(xaxis);
F=f(tgrid,xgrid);
e=ones(N-1,1);
D=spdiags([e -2*e e],-1:1,N-1,N-1)/(h^2);
I=speye(N-1);
% f0(1)=f0(1)+0;%dirichlet bd of 0
% f0(1)=f0(end)+0;%dirichlet bd of 0

if theta==0
for i=1:M
u(:,i+1)=(I+a*tau*D)*u(:,i)+tau*F(:,i);
end
else
for i=1:M
A=(I-theta*a*tau*D);
b=(I+(1-theta)*a*tau*D)*u(:,i)+tau*(theta*F(:,i+1)+(1-theta)*F(:,i));
u(:,i+1)=A\b;
end
end
uprecise=u_precise(tgrid,xgrid);
err=norm((u(:,end)-uprecise(:,end))./uprecise(:,end),2)
if(plt)
surf(taxis,xaxis,u);
%surf(taxis,xaxis,u_precise(tgrid,xgrid));
zlabel('u');
xlabel('t');
ylabel('x');
setfigure;
end
end

