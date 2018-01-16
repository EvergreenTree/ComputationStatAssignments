function [u,t] = EulerExplicit(f,t0,u0,dt,T)
i=1;
N=floor((T-t0)/dt);
u=zeros(1,N);u(1)=u0;
t=zeros(1,N);t(1)=t0;
while(t<T)
u(i+1)=u(i)+dt*f(t(i),u(i));
t(i+1)=t(i)+dt;
i=i+1;
end

end
