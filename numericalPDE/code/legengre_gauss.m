function [xp,wt] = legengre_gauss(n)
b=1./sqrt(4-(1:n).^(-2));
[V,D]=eig(diag(b,1)+diag(b,-1));
[xp,idx]=sort(diag(D));
wt=2*V(1,idx).^2;
end

