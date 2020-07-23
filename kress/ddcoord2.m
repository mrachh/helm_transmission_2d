function ddx2=ddcoord2(t)
%type means the object used
%0->reference ball of radius 1
%1->kite
%2->pear
%3->peanut
%4->iterative Newton
%5->ellipse
%6->circle radius 1 center (1,1)
%7->circle with radius 2
%8->ellipse(20,1)
%9->gear
%10->iteration high frequancy
%11->iteration with center out of (0,0)
%12->iteration with spline
%13->iteration with spline and center out of (0,0)
global type_obj
%global sp
%test for speed
global uu
global DDSP
global H
global DH
global DDH
global NT
global radius
% global TA

switch type_obj
    case 0
        ddx2=-1*sin(t);
    case 1
        ddx2=-1.5*sin(t);
    case 2
        ddx2=-1.8*sin(3*t).*cos(t)-(2+3*cos(3*t)).*sin(t);
    case 3
        ddx2=-9*sin(t).*sin(t).*sin(t).*cos(t).*cos(t)./((3*cos(t).*cos(t)+1).^(1.5))+(3*sin(t).*sin(t).*sin(t)-9*sin(t).*cos(t).*cos(t)).*(1./sqrt(3*cos(t).*cos(t)+1))-sqrt(3*cos(t).*cos(t)+1).*sin(t);                
    case 4
        %test for speed
        ind=find(uu==t);
        ddx2=DDSP(2,ind);
        %ddx2=[0 1]*Bsplineval(ddsp,t);
    case 5
        ddx2=-4*sin(t);
    case 6
        ddx2=-1*sin(t);
    case 7
        ddx2=-2*sin(t);
    case 8
        ddx2=-sin(t);
    case 9
        ddx2=2*sin(t).*(-(1+0.2*cos(9*t))-9*1.8*cos(9*t))+2*cos(t).*(-1.8*sin(9*t)+(-0.2*9*sin(9*t)));
    case 10
        TA=trig_int(t,NT);
        ddx2=transpose(TA*(DDH-H)).*sin(t)+2*transpose(TA*DH).*cos(t);
    case 11
        TA=trig_int(t,NT);
        ddx2=transpose(TA*(DDH-H)).*sin(t)+2*transpose(TA*DH).*cos(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        ddx2=(der_val_spline(t,H,NT,2)-val_spline(t,H,NT)).*sin(t)+2*der_val_spline(t,H,NT,1).*cos(t);
    case 13
        ddx2=(der_val_spline(t,H,NT,2)-val_spline(t,H,NT)).*sin(t)+2*der_val_spline(t,H,NT,1).*cos(t);
    case 15
        ddx2=2*sin(t).*(-(1+0.2*cos(15*t))-45*cos(15*t))+4*cos(t).*(-3*sin(15*t));
    case 17
        ddx2=2*sin(t).*(-(1+0.2*cos(7*t))-0.2*7*7*cos(7*t))+4*cos(t).*(-1.4*sin(7*t));
    case 18
        ddx2=1.5*sin(t);
    case 20
        ddx2 = ( 0 - 16 * radius * cos(4*t) ).*sin(t) +  ( 0 - 4 * radius * sin(4*t) ).*cos(t)+...
               ( 0 - 4 * radius * sin(4*t) ).*cos(t) + ( 1 + radius * cos(4*t) ).*-sin(t);
end
