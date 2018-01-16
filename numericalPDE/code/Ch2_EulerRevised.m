function [u,t] = Ch2_EulerRevised(f,t0,u0,dt,T)
%central diff applied to one explicit single step (using half time diff)
s=t0;
i=1;
N=floor((T-t0)/dt);
u=zeros(1,N);u(1)=u0;
t=zeros(1,N);t(1)=t0;
while(t<T)
t(i+1)=t(i)+dt;
u(i+1)=u(i)+dt*f(t(i),u(i)+dt/2*f(t(i),u(i)));
i=i+1;
end
end

