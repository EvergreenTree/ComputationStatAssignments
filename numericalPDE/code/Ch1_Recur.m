function u=Ch1_Recur(u0,u1,a,b,N)
u=zeros(N,1);
u(1)=u0;u(2)=u1;
for i=1:N-2
    u(i+2)=a*u(i+1)+b*u(i);
end
end