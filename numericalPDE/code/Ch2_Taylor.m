function [u,t] = Ch2_Taylor(f,t0,u0,dt,T,p)
%Taylor ODE Method, input order p = 2, 3
dft=diff(f,'t');dfu=diff(f,'u');
d2ft=diff(dft,'t');d2fu=diff(dfu,'u');dftu=diff(dft,'u');
df=dft+f.*dfu;
d2f=d2ft+2.*f.*dftu+f.^2.*d2fu+dfu.*df;
i=1;
N=floor((T-t0)/dt);
u=zeros(1,N);u(1)=u0;
t=zeros(1,N);t(1)=t0;
switch p
case 2
while(t<=T)
t(i+1)=t(i)+dt;
u(i+1)=u(i)+dt.*f(t(i),u(i))+.5*dt.^2.*df(t(i),u(i));
i=i+1;
end
case 3
while(t<=T)
t(i+1)=t(i)+dt;
u(i+1)=u(i)+dt.*f(t(i),u(i))+.5*dt.^2.*df(t(i),u(i))+dt.^3.*d2f(t(i),u(i))/6;
i=i+1;
end
end
end

