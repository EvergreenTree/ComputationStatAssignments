function [S] = QuadAimpson(f,a,b,n)
[xp,wt]=legengre_gauss(n);
x=(b-a)/2*xp+(b+a)/2;
S=wt*f(x)*(b-a)/2;
end

