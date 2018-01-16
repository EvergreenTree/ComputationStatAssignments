function [u,t] = Ch2_RungeKutta(f,t0,u0,dt,T,method)

t=t0:dt:T;n=length(t);
u=zeros(n,1);u(1)=u0;

for i=1:(n-1)
switch method
case 2
k=zeros(2,1);
k(1)=f(t,u(i));
k(2)=f(t+dt,u(i)+dt*k(1));
p=[1 1]/2;
case 3
k=zeros(3,1);
k(1)=f(t,u(i));
k(2)=f(t+dt/2,u(i)+dt/2*k(1));
k(3)=f(t+dt,u(i)-dt*k(1)+2*dt*k(2));
p=[1 4 1]/6;
case 4
k=zeros(4,1);
k(1)=f(t,u(i));
k(2)=f(t+dt/3,u(i)+dt/3*k(1));
k(3)=f(t+dt*2/3,u(i)-dt/3*k(1)+dt*k(2));
k(4)=f(t+dt,u(i)+dt*k(1)-dt*k(2)+dt*k(3));
p=[1 3 3 1]/8;
end
u(i+1)=u(i)+p*k*dt;
end

end
