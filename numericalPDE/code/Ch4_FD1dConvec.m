function U = Ch4_FD1dConvec(u0,r,Nt,method,ub)
%--------------------------------------------------------------------------
% U = Ch4_FD1dConvec(u0,r,Nt,method,ub)
% 
% solve the numerical solution of the convection equation
%   u_t + c*u_x = 0  in (0,T)*(0,1)
%   u(0,x) = u0
%   u(t,0) = ub
% where supp(u0)\subset (0,1), i=2:Nx; c>0, ub on left
% 
% Example by Wenbin Chen & Xinming Wu
%               2012-10-01
%--------------------------------------------------------------------------

iub = 0; if nargin==5, iub=1; end

Nx = length(u0)-1;

U = zeros(Nx+1,Nt+1);
U(:,1) = u0;

i = 2:Nx;
for it=1:Nt
    u0 = U(:,it);
    u1 = zeros(Nx+1,1);

    switch method
        case 'ftbs'
            i = 2:Nx+1;
            u1(i) = r*u0(i-1) + (1-r)*u0(i);
        case 'ftfs'
            i = 1:Nx;
            u1(i) = (1+r)*u0(i) - r*u0(i+1);
        case 'ftcs'
            u1(i) = r/2*u0(i-1) + u0(i) - r/2*u0(i+1);
        case 'lax'
            u1(i) = (1+r)/2*u0(i-1) + (1-r)/2*u0(i+1);
        case 'lax-wendroff'
            u1(i) = u0(i) - r/2*(u0(i+1)-u0(i-1)) + r^2/2*(u0(i+1)-2*u0(i)+u0(i-1));
        case 'beam-warming-left'
            i = 3:Nx+1;
            u1(i) = u0(i) - r/2*(3*u0(i)-4*u0(i-1)+u0(i-2)) + r^2/2*(u0(i)-2*u0(i-1)+u0(i-2));
        case 'beam-warming-right'
            i = 1:Nx-1;
            u1(i) = u0(i) - r/2*(-3*u0(i)+4*u0(i+1)-u0(i+2)) + r^2/2*(u0(i)-2*u0(i+1)+u0(i+2));
        case 'leapfrog'
            if it==1 % lax-wendroff
                u1(i) = u0(i) - r/2*(u0(i+1)-u0(i-1)) + r^2/2*(u0(i+1)-2*u0(i)+u0(i-1));
            else
                um = U(:,it-1);
                u1(i) = um(i) - r*(u0(i+1)-u0(i-1));
            end
    end
    
    U(:,it+1) = u1;
    if iub, U(1,it+1) = ub(it+1); end
    
end
