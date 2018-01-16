function [u,t] = Ch2_Gear(f,t0,u0,dt,T,method)
%solve u'=f by differentiating u, single step
t=t0:dt:T;n=length(t);
u=zeros(n,1);
opt=optimset('Display','off');
switch method
case 2
u(1:2)=u0;
for i=2:(n-1)
u(i+1)=fsolve(@(x)x-2/3*([2 -.5].*[u(i),u(i-1)])-f(t(i+1),x)*dt,...
        u(i),opt);
end
case 3
u(1:3)=u0;
for i=3:(n-1)
u(i+1)=fsolve(@(x)x-6/11*([3 -1.5 1/3].*[u(i),u(i-1),u(i-2)])-f(t(i+1),x)*dt,...
        u(i),opt);
end
case 4
u(1:4)=u0;
for i=4:(n-1)
u(i+1)=fsolve(@(x)[-25/12 4 -3 4/3 -1/4].*[x,u(i),u(i-1),u(i-2),u(i-3)]+f(t(i+1),x)*dt,...
        u(i),opt);
end
end

