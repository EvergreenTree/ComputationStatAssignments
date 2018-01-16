f=@(t,u)u+u.^.5;
uinit=[0.5 1.5]
u0=uinit(1);
dt=8/50;
T=8;
t0=0;
u=Ch2_EulerExplicit(f,t0,u0,dt,T);


%%%hw2
U=@(t,t0)exp(t-t0)-t-1;
f=@(t,u)u+t;
u0=1;
dt=1/5;
T=1;
t0=-1;
[u,t]=Ch2_EulerExplicit(f,t0,u0,dt,T);
hold on;
plot(t,u);
[u,t]=Ch2_EulerImplicit(f,t0,u0,dt,T);
plot(t,u);
[u,t]=Ch2_EulerImproved(f,t0,u0,dt,T);
plot(t,u);
[u,t]=Ch2_EulerRevised(f,t0,u0,dt,T);
plot(t,u);
plot(t,U(t,t0),'.-');
title('Euler Method du/dt=u+t');
legend('Explicit','Implicit','Improved','revised','actual');
hold off;


%%% Convergence Domain of Improved Euler Method
%%$revised
t=exp((0:0.01*pi:2*pi)*1i);
n=length(t);
y=zeros(n,1);
z=zeros(n,1);
f=@(x,t)1+x+x.^2./2-t;
for i=1:n
temp=solve(@(x)f(x,t(i)));
y(i)=temp(1);z(i)=temp(2);
end
plot([y;z],'o')

f=@(x,t)(1+.5.*x)-(1-.5.*x).*t;
for i=1:n
y(i)=solve(@(x)f(x,t(i)));
end
plot(y,'o');

s=((t-1)./(t+1).*2);
fill([0 0 20 20],[-20 20 20 -20],'b');xlim([-10 10]);ylim([-10 10])



%%$ Taylor (week5)
syms u t;
f(t,u)=u-u.^3;
t0=0;u0=1.3;dt=1/5;T=2;
ff=@(x)exp(x)./sqrt(1/u0^2-1+exp(2*x));
[u2,t2]=Ch2_Taylor(f,t0,u0,dt,T,2);
[u3,t3]=Ch2_Taylor(f,t0,u0,dt,T,3);
[u1,t1]=Ch2_EulerExplicit(f,t0,u0,dt,T);
t4=0:dt:T;
u4=ff(t4);
hold on;plot(t4,u4,'.-');plot(t1,u1);plot(t2,u2);plot(t3,u3);hold off;
legend('actual','EulerExplicit','Taylor2','Taylor3');

%%%Runge Kutta 2~4 (w5)
f=@(t,u)u-u.^3;
t0=0;u0=.5;dt=.5;T=2;
ff=@(x)exp(x)./sqrt(1/u0^2-1+exp(2*x));
[u2,t2]=Ch2_RungeKutta(f,t0,u0,dt,T,2);%Explicit Improved Euler
[u3,t3]=Ch2_RungeKutta(f,t0,u0,dt,T,3);
[u4,t4]=Ch2_RungeKutta(f,t0,u0,dt,T,4);
t1=0:dt:T;
u1=ff(t1);
hold on;plot(t1,u1,'.-');plot(t2,u2);plot(t3,u3);plot(t4,u4);hold off;
legend('actual','RK2','RK3','RK4');

f=@(t,u)u-u.^3;
t0=0;u0=.8;T=10;
ff=@(x)exp(x)./sqrt(1/u0^2-1+exp(2*x));
N=11;%%Takes a while
e2=zeros(N,1);
e3=zeros(N,1);
e4=zeros(N,1);
for i=1:N
dt=2^(-i);
[u2,t2]=Ch2_RungeKutta(f,t0,u0,dt,T,2);%Explicit Improved Euler
e2(i)=abs(u2(end)-ff(t2(end)));
[u3,t3]=Ch2_RungeKutta(f,t0,u0,dt,T,3);
e3(i)=abs(u3(end)-ff(t3(end)));
[u4,t4]=Ch2_RungeKutta(f,t0,u0,dt,T,4);
e4(i)=abs(u4(end)-ff(t4(end)));
end
semilogy(e2);
hold on
semilogy(e3);
semilogy(e4);
hold off;
legend('RK2','RK3','RK4');



%%%Adams Error with different starting value estimation
f=@(t,u)2.*u;t0=0;T=10;u0=1;
uu=@(t)exp(2*t);
N=12;
ee=zeros(N,1);
eE=zeros(N,1);
e2=zeros(N,1);
e3=zeros(N,1);
e4=zeros(N,1);
for i=1:N
dt=2^(-i);
t=t0:dt:T;
u0e=uu(t(1:4));%exact
u0E=Ch2_EulerExplicit(f,t0,u0,dt,t0+3*dt);%EulerExplicit
u02=Ch2_RungeKutta(f,t0,u0,dt,t0+3*dt,2);%RungeKutta2
u03=Ch2_RungeKutta(f,t0,u0,dt,t0+3*dt,3);%RungeKutta3
u04=Ch2_RungeKutta(f,t0,u0,dt,t0+3*dt,4);%RungeKutta4
[ue,t]=Ch2_AdamsExplicit(f,t0,u0e,dt,T,4);
[uE,t]=Ch2_AdamsExplicit(f,t0,u0E,dt,T,4);
[u2,t]=Ch2_AdamsExplicit(f,t0,u02,dt,T,4);
[u3,t]=Ch2_AdamsExplicit(f,t0,u03,dt,T,4);
[u4,t]=Ch2_AdamsExplicit(f,t0,u04,dt,T,4);
uexact=uu(t);
ee(i)=abs(ue(end)-uexact(end));
eE(i)=abs(uE(end)-uexact(end));
e2(i)=abs(u2(end)-uexact(end));
e3(i)=abs(u3(end)-uexact(end));
e4(i)=abs(u4(end)-uexact(end));
end
semilogy(ee,'.-');hold on;
semilogy(eE);
semilogy(e2);
semilogy(e3);
semilogy(e4);
legend('u0 = exact','u0 = Euler','u0 = RK2','u0 = RK3','u0 = RK4')

%%
n=2^20
u=zeros(n,1);
u(1:2)=[0;0.1];
for i=3:n
u(i)=2*u(i-1) -u(i-2);
end
plot(u)

for i=3:n
u(i)=2*u(i-1);
end
plot(u-.1)


%%week 7 P118
t0=0;u0=0;dt=.1;T=10;
tt=linspace(t0,T,1000);
%u0=1;
U=@(t,lambda,u0)(sin(t)+lambda*cos(t))*lambda/(1+lambda^2)+ ...
    (u0-lambda^2/(1+lambda^2))*exp(-lambda*t);
f=@(t,u,lambda)lambda*(-u+cos(t));

for i=0:4
lambda=10^i;
uu=U(tt,lambda,u0);
plot(tt,uu,'r.');
hold on;
[u t]=Ch2_EulerExplicit(@(t,u)f(t,u,lambda),t0,u0,dt,T);
plot(t,u,'g-');
[u t]=Ch2_EulerImplicit(@(t,u)f(t,u,lambda),t0,u0,dt,T);
plot(t,u,'b-');
hold off;
legend('exact','EulerExplicit','EulerImplicit');
setfigure;ylim([-1 1]);
pause;
end

lambda=1000;
u0=1;
uu=U(tt,lambda,u0);
plot(tt,uu,'r.');hold on;
u0=U(t0:dt:(t0+3*dt),lambda,u0);
[u1 t1]=Ch2_Gear(@(t,u)f(t,u,lambda),t0,u0,dt,T,4);
plot(t1,u1,'g-');
[u2 t2]=Ch2_AdamsExplicit(@(t,u)f(t,u,lambda),t0,u0,dt,T,4);
plot(t2,u2,'b-');
hold off;legend('exact','Gear(Implicit)','AdamsExplicit')
ylim([-1 1]);

lambda=10;u0=1;
uu=U(tt,lambda,u0);
e1=zeros(7,1);
e2=zeros(7,1);
for k=1:7
dt=2^(-k);
u00=U(t0:dt:(t0+3*dt),lambda,u0);
[u1 t1]=Ch2_Gear(@(t,u)f(t,u,lambda),t0,u00,dt,T,4);
e1(k)=abs(u1(end)-uu(end));
[u2 t2]=Ch2_AdamsExplicit(@(t,u)f(t,u,lambda),t0,u00,dt,T,4);
e2(k)=abs(u2(end)-uu(end));
end
plot(e2,'o');hold on;plot(e1,'*');setfigure;legend('Adams','Gear');
ylim([0 1]);hold off

