function ddx1=ddcoord1(t)
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
        ddx1=-1*cos(t);
    case 1
        ddx1=-cos(t)-2.6*cos(2*t);
    case 2
        ddx1=1.8*sin(3*t).*sin(t)-(2+3*cos(3*t)).*cos(t);        
    case 3
        ddx1=-9*sin(t).*sin(t).*cos(t).*cos(t).*cos(t)./((3*cos(t).*cos(t)+1).^(1.5))+(-3*cos(t).*cos(t).*cos(t)+9*sin(t).*sin(t).*cos(t))./sqrt(3*cos(t).*cos(t)+1)-sqrt(3*cos(t).*cos(t)+1).*cos(t);        
    case 4
        %test for speed
        ind=find(uu==t);
        ddx1=DDSP(1,ind);
        %ddx1=[1 0]*Bsplineval(ddsp,t);
    case 5
        ddx1=-3*cos(t);
    case 6
        ddx1=-1*cos(t);
    case 7
        ddx1=-2*cos(t);
    case 8
        ddx1=-20*cos(t);
    case 9
        ddx1=2*(-1.8*9*cos(9*t)-(1+0.2*cos(9*t))).*cos(t)-2*(-1.8*sin(9*t)+(-0.2*9*sin(9*t))).*sin(t);
    case 10
        TA=trig_int(t,NT);
        ddx1=transpose(TA*(DDH-H)).*cos(t)-2*transpose(TA*DH).*sin(t);
    case 11
        TA=trig_int(t,NT);
        ddx1=transpose(TA*(DDH-H)).*cos(t)-2*transpose(TA*DH).*sin(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        ddx1=(der_val_spline(t,H,NT,2)-val_spline(t,H,NT)).*cos(t)-2*der_val_spline(t,H,NT,1).*sin(t);
    case 13
        ddx1=(der_val_spline(t,H,NT,2)-val_spline(t,H,NT)).*cos(t)-2*der_val_spline(t,H,NT,1).*sin(t);
    case 20
        ddx1 = ( 0 - 16 * radius * cos(4*t) ).*cos(t) + ( 0 - 4 * radius * sin(4*t) ).*-sin(t) +  ...
            ( 0 - 4 * radius * sin(4*t) ).*-sin(t) + ( 1 + radius * cos(4*t) ).*-cos(t);        
end
