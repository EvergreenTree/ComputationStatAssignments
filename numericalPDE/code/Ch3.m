%% week8
f=@(x,k)-cos(k*x)/k^2+(cos(k)-1)/(k^2*sin(k))*sin(k*x)+1/k^2;
hold on;
for i=1:10
k=i;
f1=@(x)f(x,k);
xx=0:0.01:1;uu=f1(xx);
plot(xx,uu);
end
hold off;legend('k=1','k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10')

%% solve -u''+bu'+cu=1
%f=@(x)sin(2*pi*x);
f=@(x)-x./x;
hold on;
for i=1:10
[u x]=Ch3_fd1dfr(b,c,f,2^i,0,0);
plot(x,u);
end
setfigure;


hold on;
for k=1:3
b=-20+10*k;
[u x]=Ch3_fd1dfr(b,c,f,2^i,0,0);
plot(x,u);
end
legend('b=-10;c=0','b=0;c=0','b=10;c=0');setfigure;

b=10;c=10;
e=zeros(10,1);
h=2.^(-3:-1:-12);
for i=1:10
N=1/h(i);
[u x]=Ch3_fd1dfr(b,c,f,N,0,0);
u0=Ch3_preciseval(x,b,c,N);
e(i)=norm(u-u0);
end
%hold on;plot(x,u);plot(x,u0);
loglog(h,e);title('error-to-stepsize');

% [u x]=Ch3_fd1dbr(b,c,f,2^i,0,0);
% plot(x,u);
% [u x]=Ch3_fd1dcr(b,c,f,2^i,0,0);
% plot(x,u);

%% matlab built-in
twoode=@(x,y)[ y(2); b*y(2)+c*y(1)-1 ];
twobc=@(ya,yb)[ ya(1); yb(1) + 2 ];
solinit = bvpinit(linspace(0,4,5),[1 0]);
sol = bvp4c(twoode,twobc,solinit);
x = linspace(0,4);
y = deval(sol,x);
plot(x,y(1,:))

%% chebyshev (play)
f=inline('cos(k*acos(x))')
xx=linspace(-1,1);
for i=1:7
plot(xx,f(xx))
end
%%%


%% -u''+u=f
n=3
I=speye(n)
e=ones(n,1)
D=spdiags([-e 2*e -e],-1:1,n,n)
A=kron(I,D)+kron(D,I)

B=A+speye(n*n)


I0=speye(n^2)

[V,E]=eig(B)

%% fd2d
p=15;
N=100:100:100*p;
err=zeros(p,1);
for i=1:p
[~,~,~,err(i)]=Ch3_fd2d(N(i));
end
loglog(1./N,err);
setfigure;