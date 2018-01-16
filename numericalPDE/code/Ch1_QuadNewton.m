function [S] = Ch1_QuadNewtion(f,n,a,b)
x=linspace(a,b,n+1);
switch n
case 1 
wt=[1 1]/2;
case 2 
wt=[1 4 1]/6;
case 3 
wt=[1 3 3 1]/8;
case 4
wt=[7 32 12 32 7]/90;
end
S=sum(wt.*f(x))*(b-a);
end
