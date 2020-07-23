function dx1=dcoord1(t)
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
%15->object with 15 petals
global type_obj
%global sp
%test for speed
global uu
global DSP
global H
global DH
global NT
global radius
% global TA

switch type_obj
    case 0
        dx1=-1*sin(t);
    case 1
        dx1=-sin(t)-1.3*sin(2*t);
    case 2
        dx1=-(0.9*sin(3*t)).*cos(t)-(2+0.3*cos(3*t)).*sin(t);
    case 3
        dx1=-3*sin(t).*cos(t).*(1./sqrt(3*cos(t).*cos(t)+1)).*cos(t)-sqrt(3*cos(t).*cos(t)+1).*sin(t);
    case 4
        %test for speed
        ind=find(uu==t);
        dx1=DSP(1,ind);
        %dx1=[1 0]*Bsplineval(dsp,t);
    case 5
        dx1=-3*sin(t);
    case 6
        dx1=-1*sin(t);
    case 7
        dx1=-2*sin(t);
    case 8
        dx1=-20*sin(t);
    case 9
        dx1=2*(-1.8*sin(9*t)).*cos(t)-2*(1+0.2*cos(9*t)).*sin(t);
    case 10
        TA=trig_int(t,NT);
        dx1=transpose(TA*DH).*cos(t)-transpose(TA*H).*sin(t);
    case 11
        TA=trig_int(t,NT);
        dx1=transpose(TA*DH).*cos(t)-transpose(TA*H).*sin(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        dx1=der_val_spline(t,H,NT,1).*cos(t)-val_spline(t,H,NT).*sin(t);
    case 13
        dx1=der_val_spline(t,H,NT,1).*cos(t)-val_spline(t,H,NT).*sin(t);
    case 15
        dx1=2*(-3*sin(15*t)).*cos(t)-2*(1+0.2*cos(15*t)).*sin(t);
    case 17
        dx1=2*(-1.4*sin(7*t)).*cos(t)-2*(1+0.2*cos(7*t)).*sin(t);
    case 18
        dx1=+sin(t)+1.3*sin(2*t);
    case 19
        dx1=3/4*cos(2*t).*cos(t-pi/4)-3/4*(1+sin(2*t)/2).*sin(t-pi/4);
    case 20
        dx1 = ( 0 - 4 * radius * sin(4*t) ).*cos(t) + ( 1 + radius * cos(4*t) ).*-sin(t);        
end
