%%convergence analysis
err=zeros(4,1);
a=1
T=.5
X=1
f=@(t,x)t-t;
u0=@(x)sin(pi*x);
u_precise=@(t,x)sin(pi*x).*exp(-pi^2*t);
plt=0;

theta=1
%%%t
M=2.^((1:4)+7);
for i=1:4
[~,~,~,err(i)]=Ch4_fd1dheat(theta,a,T,X,f,u0,M(i),16,u_precise,plt);
end

loglog(T./M,err);
hold on;
loglog(2.^(-((2:5)+7)),2.^(-(2:5))/10);
legend('err','line of slope 1');
xlabel('t stepsize');
ylabel('error');
setfigure;

%%%x
N=2.^((1:4)+1)
for i=1:4
[~,~,~,err(i)]=Ch4_fd1dheat(theta,a,T,X,f,u0,10000,N(i),u_precise,plt);
end

loglog(X./N,err);
hold on;
loglog(2.^(-((2:5))),2.^(-(2:2:8)));
legend('err','line of slope 2');
xlabel('x stepsize');
ylabel('error');
setfigure;

%%convect

method={'ftbs','leapfrog','lax','lax-wendroff'}

f=@(x)(x>.25&x<.5).*4.*(x-.25)+(x>.5&x<.75).*4.*(.75-x);
%f=@(x)(x>.25&x<.75)*1;

u0=f(xx);


r=.5
for i=1:4
figure(i)
U=Ch4_FD1dConvec(u0,r,50,method{i},zeros(50));
surf(U)
xlabel('t')
ylabel('x')
title(sprintf('%s, %s%#0.1f',method{i},'r = ',r))
print(i,'-depsc',sprintf('%s%d','2',i));
end
