function [U X Y err]=Ch3_fd2d(n,W,f)
%solve -u''+u=f in \Omega, u=1 on boundary
if(nargin<1)
n=100;
end
h=1/n;
if(nargin<2)
I=speye(n-1);
e=ones(n-1,1);
D=spdiags([-e 2*e -e],-1:1,n-1,n-1);
A=(kron(I,D)+kron(D,I))/h^2+speye((n-1)^2);%2d span on 1d
end
if(nargin<3) 
f=@(x, y) (1 + 8*pi^2) .* sin(2*pi*x) .* sin(2*pi*y)+1;
%f=@(x,y) 2.*y.*(1-y)+2.*x.*(1-x)+x.*(1-x).*y.*(1-y)+1;
end

h=1/n;
x=(1:1:(n-1))/n;
[X Y]=meshgrid(x);%2d domain
F = f(X,Y);
F([1 end],:)=F([1 end],:)+1/h^2;
F(:,[1 end])=F(:,[1 end])+1/h^2;
F = F(:);%2d span on 1d

uf = A\F;%solution
U = reshape(uf,n-1,n-1);%2d solution

%% plot
u0=@(x, y) sin(2*pi*x) .* sin(2*pi * y)+1;
%u0=@(x, y) x.*(1-x).*y.*(1-y)+1;
U0=u0(X,Y);
surf(X,Y,U);
hold on;
surf(X,Y,U0);
hold off

err=max(abs(U0(:)-U(:)))
%err=norm(U0-U)
end


