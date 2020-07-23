function K=double2D(k,u,v)%S,coord_1,coord_2,dcoord_1,dcoord_2,ddcoord_1,ddcoord_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the double layer potential in 2D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K\phi(x)=int(\partial D)(\partial\Phi(x,y))\(\partial n(y))\phi(y)ds(y) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K\phi(t)=int_(0)^(2pi)L(t,\tau)\phi(\tau)d\tau %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L(t,\tau)=L_1{t,\tau)ln(4sin^2((t-\tau)/2))+L_2(t,\tau) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_1(t,\tau)=k/(2pi)besselj(1,kr(t,\tau))/r((t,\tau)*         %  
% [(x_1(t)-x_1(\tau))x_2'(\tau)-(x_2(t)-x_2(\tau))x_1'(\tau)]* %
% sqrt(x_1'(\tau)^2+x_2'(\tau)^2)                              %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_2(t,\tau)=L(t,\tau)-L_1(t,\tau)ln(4sin^2((t-\tau)/2)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_2(t,t)=1/(2pi)(x_1'(t))x_2"(t)-x_2'(t))x_1"(t))/(x_1'(t)^2+x_2'(t)^2) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_1(t,t)=0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \int_(0)^(2pi)ln(4sin^2((t-\tau)/2))f(\tau)d\tau=\sum_(j=0)^(2n-1) %
% R_j^(n)(t)f(t_j)                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_j^(n)(t)=-2pi/n\sum_(m=1)^(n-1)cos(m(t-t_j))/m-pi/(n^2) %
% cos(n(t-t_j))                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(u,2);
n1=size(v,2);
K=zeros(n,n1);
for p=1:n
    x11=coord1(u(p));
    x12=coord2(u(p));
    for l=1:n1       
        x21=coord1(v(l));
        x22=coord2(v(l));
        dx21=dcoord1(v(l));
        dx22=dcoord2(v(l));        
        ddx21=ddcoord1(v(l));
        ddx22=ddcoord2(v(l));                
        r=sqrt((x11-x21)^2+(x12-x22)^2);
        dr=dx21^2+dx22^2;
        if abs(u(p)-v(l))<10^(-15)%I needed to put this, because a couple of values didn't match
            l1=0;                     
            l2=(1/(2*pi))*(dx21*ddx22-dx22*ddx21)/(dr);            
        else
            l1=(k/(2*pi))*((x11-x21)*dx22-(x12-x22)*dx21)*besselj(1,k*r)/r;
            l2=(1i*k/2)*besselh(1,1,k*r)*((-x11+x21)*dx22-(-x12+x22)*dx21)/r - l1*log(4*(sin((u(p)-v(l))/2)^2));            
        end
        K(p,l)=-quadpotential2D(n1/2,u(p),v(l))*l1-(2*pi/n1)*l2;
    end
end
