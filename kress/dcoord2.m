function dx2=dcoord2(t)
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
global DSP
global H
global DH
global NT
global radius
% global TA

switch type_obj
    case 0
        dx2=1*cos(t);
    case 1
        dx2=1.5*cos(t);
    case 2
        dx2=-0.9*sin(3*t).*sin(t)+(2+0.3*cos(3*t)).*cos(t);
    case 3
        dx2=-3*sin(t).*cos(t).*(1./sqrt(3*cos(t).*cos(t)+1)).*sin(t)+sqrt(3*cos(t).*cos(t)+1).*cos(t);
    case 4
        %test for speed
        ind=find(uu==t);
        dx2=DSP(2,ind);
        %dx2=[0 1]*Bsplineval(dsp,t);
    case 5
        dx2=4*cos(t);
    case 6
        dx2=1*cos(t);
    case 7
        dx2=2*cos(t);
    case 8
        dx2=cos(t);
    case 9
        dx2=2*(-1.8*sin(9*t)).*sin(t)+2*(1+0.2*cos(9*t)).*cos(t);
    case 10
        TA=trig_int(t,NT);
        dx2=transpose(TA*DH).*sin(t)+transpose(TA*H).*cos(t);
    case 11
        TA=trig_int(t,NT);
        dx2=transpose(TA*DH).*sin(t)+transpose(TA*H).*cos(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        dx2=der_val_spline(t,H,NT,1).*sin(t)+val_spline(t,H,NT).*cos(t);
    case 13
        dx2=der_val_spline(t,H,NT,1).*sin(t)+val_spline(t,H,NT).*cos(t);
    case 15
        dx2=2*(-3*sin(15*t)).*sin(t)+2*(1+0.2*cos(15*t)).*cos(t);
    case 17
        dx2=2*(-1.4*sin(7*t)).*sin(t)+2*(1+0.2*cos(7*t)).*cos(t);
    case 18
        dx2=-1.5*cos(t);
    case 19
        dx2=cos(2*t).*sin(t-pi/4)+(1+sin(2*t)/2).*cos(t-pi/4);
    case 20
        dx2 = ( 0 - 4 * radius * sin(4*t) ).*sin(t) + ( 1 + radius * cos(4*t) ).*cos(t);
end
