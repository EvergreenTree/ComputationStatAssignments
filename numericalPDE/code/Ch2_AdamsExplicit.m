function [u,t] = Ch2_AdamsExplicit(f,t0,u0,dt,T,method)
%solve u'=f by integrating f, singlestep
t=t0:dt:T;n=length(t);
u=zeros(n,1);
switch method
case 2
u(1:2)=u0;
for i=2:(n-1)
u(i+1)=u(i)+[1.5 -.5]*[f(t(i),u(i));f(t(i-1),u(i-1))]*dt;
end
case 3
u(1:3)=u0;
for i=3:(n-1)
u(i+1)=u(i)+[23 -16 5]*[f(t(i),u(i));f(t(i-1),u(i-1));f(t(i-2),u(i-2))]/12*dt;
end
case 4
u(1:4)=u0;
for i=4:(n-1)
u(i+1)=u(i)+[55 -59 37 -9]*[f(t(i),u(i));f(t(i-1),u(i-1));f(t(i-2),u(i-2));f(t(i-3),u(i-3))]/24*dt;
end
end

end
