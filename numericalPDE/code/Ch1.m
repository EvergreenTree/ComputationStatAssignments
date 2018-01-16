%%%%1.1
u=Ch1_Recur(.1,.1,3,-2,100);
semilogy(u,'*');
plot((u(2:end)-0.1)./(u(1:end-1)-0.1),'*')

%%%Ex2
u=Ch1_Recur(.1,.1,3/2,-1/2,100);
plot(u,'*');
plot((u(2:end)-0.1)./(u(1:end-1)-0.1),'*')

x=single(0.1)
format hex
x
format short
x

0.3-0.2-0.1

x=1;
while 1+x/2>1
    x=x/1;
end

f=inline('(1-cos(x))/x^2')
i=1
x(1)=1
while 1+x/2>1
    x(i+1)=x(i)/2;
end

c1=zeros(98,1);
for N=3:100
A1=2*diag(ones(N-2,1),-2)-3*diag(ones(N-1,1),-1)+eye(N);
A1(2,1)=0;
A2=A1;A2(1,N)=-1;
c1(N-2)=cond(A1);
c2(N-2)=cond(A2);
end
figure(1);clf;semilogy(c1,'-*');
figure(2);clf;plot(c2,'-*')


%%%Ex1
f=@(x)(1-cos(x))./x./x
x=(-100:100)*1e-8;
figure(1);plot(f(x));
x=(-10000:10000)*1e-10;
figure(2);plot(f(x));

%%%%1.2
var1=1;


%%%1.2.20 Newton-Rapson
F=@(x) x+log(x)
dF=@(x) 1+1./x;
Gtilde=@(x) x-1./dF(x)*F(x)
x=Ch1_SingleIter(Gtilde,0.1,10);
err=abs(x-x(end));
figure(1);clf;semilogy(err)
pause(1);leg_list{1}=sprintf('err_0 = ',err(1))

%%%Ex P24 2
F=@(x) x-exp(-x);
dF=@(x) 1+exp(-x);
Ftilde=@(x) x-F(x)./dF(x);
x=Ch1_SingleIter(@(x)x-F(x),0.5,24)
plot(x-lambertw(0, 1),'*');hold on;
y=Ch1_SingleIter(Ftilde,0.5,24);
plot(y-lambertw(0, 1),'o');
legend('FixedPt','Newton-Rapson');
title(sprintf('%d LALALALALALA',x(1)));

%%%Logistics Eq
for alpha = 1.4:0.01:3
    for i = 1:100
        
    end
    for i = 1:300
        
    end
end

%%%1.4.4
f=@(x,k) exp(sin(k*x)); 
f3=@(x) f(x,3);

method='uniform';
[pcoef,x]=Ch1_Interp(f3,-1,1,5,method);

xx=linspace(-1,1,100);
yy=f3(xx);
pval = polyval(pcoef,xx);

plot(xx,pval,'r',xx,f3(xx));

%%%1.4.5
a=[1 0 -2 3];
x=rand(30,1);
errlev=0.1

errcoef=Ch1_InterpFit(a,x,errlev);

%%%hw Runge
figure(1);clf();Ch1_Runge();


%%%hw Newton-Cotes Int
S=zeros(5,4);
T=zeros(5,1);
f=@(x) x+1./x;
Ch1_QuadNewton(f,4,1,4)
integral(f,1,4)
for i = 1:5
p=ones(1,i+1);
f=@(x)polyval(p,x);
T(i,1)=integral(f,1,4);
for j=1:4
S(i,j)=Ch1_QuadNewton(f,j,1,4);
end
end

%%%Gauss Int
f=@(x,n)x.^n;
S=zeros(7,5);
for n=1:7
g=@(x)f(x,n);
for nquad=0:4
    S(n,nquad+1)=QuadSimpson(g,1,4,nquad);
end
end