function u = Ch3_preciseval(x,b,c,N)
%precise value of Ch3_fd1dfr.m
if size(x,1)==1
x=x';
end
h=1/N;
l=roots([-1 b c]);
l1=l(1);l2=l(2);

u=[exp(l1*x) exp(l2*x)]*[1,1;(l1+1)*exp(l1),(l2+1)*exp(l2)]^(-1)*[1/c; 1/c]-1/c;

end

