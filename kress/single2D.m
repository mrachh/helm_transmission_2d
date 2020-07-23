function S=single2D(k,u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the single layer potential in 2D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S\phi(x)=int(\partial D)\Phi(x,y)\phi(y)ds(y) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S\phi(t)=int_(0)^(2pi)M(t,\tau)\phi(\tau)d\tau %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M(t,\tau)=M_1{t,\tau)ln(4sin^2((t-\tau)/2))+M_2(t,\tau) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M_1(t,\tau)=-1/(2pi)besselj(0,kr(t,\tau))sqrt(x_1'(\tau)^2+x_2'(\tau)^2)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M_2(t,\tau)=M(t,\tau)-M_1(t,\tau)ln(4sin^2((t-\tau)/2)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M_2(t,t)=(i/2-(-psi(1)/pi)-1/(2pi)ln((k^2)/4(x_1'(t)^2+x_2'(t)^2)) %
% sqrt(x_1'(t)^2+x_2'(t)^2)                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M1 is not special %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \int_(0)^(2pi)ln(4sin^2((t-\tau)/2))f(\tau)d\tau=\sum_(j=0)^(2n-1) %
% R_j^(n)(t)f(t_j)                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_j^(n)(t)=-2pi/n\sum_(m=1)^(n-1)cos(m(t-t_j))/m-pi/(n^2) %
% cos(n(t-t_j))                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(u,2);
n1=size(v,2);
S=zeros(n,n1);
for p=1:n
    x11=coord1(u(p));
    x12=coord2(u(p));
    for l=1:n1
        x21=coord1(v(l));
        x22=coord2(v(l));
        dx21=dcoord1(v(l));
        dx22=dcoord2(v(l));       
        r=sqrt((x11-x21)^2+(x12-x22)^2);
        dr=dx21^2+dx22^2;
        m1=-1/(2*pi)*besselj(0,k*r)*sqrt(dr);
        if abs(u(p)-v(l))<=10^(-15)%I needed to put this, because a couple of values didn't match
            m2=(1i/2+psi(1)/pi-1/(2*pi)*log(k^2*dr/4))*sqrt(dr);
        else
            m2=(1i/2)*besselh(0,1,k*r)*sqrt(dr)-m1*log(4*(sin((u(p)-v(l))/2))^2);
        end
        S(p,l)=quadpotential2D(n/2,u(p),v(l))*m1+(2*pi/n)*m2;
    end
end

