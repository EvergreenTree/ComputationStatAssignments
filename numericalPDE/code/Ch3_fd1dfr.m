function [u x]=Ch3_fd1dfr(b,c,f,N,u0,u1,diff)
%finite difference method, 1 dimension, forward, robin
%-u''+bu'+cu=f(x), 
%left-dirichlet u(0)=u0
%right-robin u'(1)+u(1)=u1
h=1/N;
e=ones(N,1);
A=spdiags([-e 2*e -e],-1:1,N,N)/(h^2);
B=spdiags([-e e],0:1,N,N)/h;%forward
x=(h:h:1)';
f0=f(x);
f0(1)=f0(1)+u0/h^2;%dirichlet bd
M=A+b*B+c*speye(N);

M(end,(end-1):end)=[-1 1]/h+[0 1];
f0(end)=0;

% f0(end)=f0(end)+u1/h/(h+1);
% A(end,end)=A(end,end)-1/(h^2*(h+1));%robin

u=(M)^(-1)*f0;
end