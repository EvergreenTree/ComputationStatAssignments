function [] = Ch1_Runge()
f=@(x) 1./(1+25*x.^2);
xx=-1:1/200:1;
yy=f(xx);

for n=3:20
x_u=-1:2/n:1;
pc_u=polyfit(x_u,f(x_u),n);
pv_u=polyval(pc_u,xx);

x_c=-cos((0:n)*pi/n);
pc_c=polyfit(x_c,f(x_c),n);
pv_c=polyval(pc_c,xx);

plot(xx,yy,'k-',xx,pv_u,'.-',xx,pv_c,'.-');
legend('f(x)','Uniform Intersection','Chebyshev Intersecton');
title(sprintf('Runge Phenomenon (n=%d)',n));
saveas(gcf,strcat(num2str(n),'.jpg'));
end
end
