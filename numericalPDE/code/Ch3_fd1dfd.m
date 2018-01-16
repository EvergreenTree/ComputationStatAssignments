function [u x]=Ch3_fd1dfd(b,c,f,N,u0,u1)
%finite difference method, 1 dimension, forward, dirichlet
%-u''+bu'+cu=f(x), 
%left-dirichlet u(0)=u0
%right-robin u(1)=0
h=1/N;
e=ones(N-1,1);
A=spdiags([-e 2*e -e],-1:1,N-1,N-1)/(h^2);
B=spdiags([-e e],0:1,N-1,N-1)/h;%forward
x=(h:h:(1-h))';
f0=f(x);
f0(end)=f0(end)-b/h*u0;%forward

f0(1)=f0(1)+u0/h^2;%dirichlet
f0(end)=f0(end)+u0/h^2;%dirichlet
% f0(end)=f0(end)+u1/(h^2);
% A(end,end)=A(end,end)-1/(h^2*(h+1));%robin


u=(A+b*B+c*speye(N-1))^(-1)*f0;
