function errcoef = Ch1_InterpFit(a,x,errlev)
x=sort(x);
y=polyval(a,x);

dataerr=rand(size(y))-.5;
ydata=y + errlev*dataerr.*y;

pcoef=polyfit(x,ydata,length(a)-1);
pval=polyval(pcoef,x);
errcoef=norm(a-pcoef)/norm(a);

plot(x,ydata,'o');hold on;
plot(x,y,'r');hold off;
end

