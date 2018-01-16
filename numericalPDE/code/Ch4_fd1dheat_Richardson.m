function [u,xaxis,taxis,err,uprecise] = Ch4_fd1dheat_Richadson(theta,a,T,X,f,u0,M,N,u_precise)
%solve u_t - a u_xx = f, t=0: u=u0
%with finite difference \theta formula

if nargin < 1
theta=0;
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
M=100;
end
if nargin < 8
N=100;
end
if nargin < 9
u_precise=@(t,x)sin(pi*x).*exp(-pi^2*t);
end

tau=T/M;
h=X/N;
u=zeros(N-1,M);
xaxis=h:h:X-h;
taxis=0:tau:T;
dt=ones(1,N-1)*tau;
tzero=zeros(1,N-1);
[tgrid xgrid]=meshgrid(taxis,xaxis);
u(:,1)=u0(xaxis);
F=f(tgrid,xgrid);
e=ones(N-1,1);
D=spdiags([e -2*e e],-1:1,N-1,N-1)/(h^2);
I=speye(N-1);
% f0(1)=f0(1)+0;%dirichlet bd of 0
% f0(1)=f0(end)+0;%dirichlet bd of 0

%use some 2nd order method(e.g. Runge-Kutta) to compute the first t step
%(in order not to pollute the algorithm):
F1=f(tzero,xaxis);
u(:,2)=(I+tau*a*D)*u(:,1)+.5*tau*(F1+f(dt,xaxis+tau*F1))';

for i=1:M-1
A=(I-theta*a*2*tau*D);
b=(I+(1-theta)*a*2*tau*D)*u(:,i)+2*tau*(theta*F(:,i+1)+(1-theta)*F(:,i));
u(:,i+2)=A\b;
end
uprecise=u_precise(tgrid,xgrid);
err=norm(u-uprecise,Inf)
surf(taxis,xaxis,u);
%surf(taxis,xaxis,u_precise(tgrid,xgrid));
zlabel('u');
xlabel('t');
ylabel('x');
end

