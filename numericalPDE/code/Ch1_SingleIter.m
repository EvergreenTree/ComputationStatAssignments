function x=Ch1_SingleIter(G,x0,N)
%if nargin<4, stop=1; end
x=zeros(N+1,1);
x(1)=x0;
for i=1:N
x(i+1)=G(x(i));
end
end