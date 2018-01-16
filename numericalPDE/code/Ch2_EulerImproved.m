function [u,t] = Ch1_EulerImproved(f,t0,u0,dt,T)
%trapezoid (f_i + f_(i+1))/2
s=t0;
i=1;
N=floor((T-t0)/dt);
u=zeros(1,N);u(1)=u0;
t=zeros(1,N);t(1)=t0;
while(t<T)
t(i+1)=t(i)+dt;
u(i+1)=fsolve(@(x) u(i)+dt*.5*(f(t(i),u(i))+f(t(i+1),x))-x,u(i));
i=i+1;
end

end
