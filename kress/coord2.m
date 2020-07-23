function x2=coord2(t)
%type means the object used
%0->reference ball of radius 1
%1->kite
%2->pear
%3->peanut
%4->iterative Newton
%5->ellipse
%6->circle of radius 1 and center (1,1)
%7->circle with radius 2
%8->ellipse(20,1)
%9->gear
%10->iteration high frequency
%11->iteration with center out of (0,0)
%12->iteration with spline
%13->iteration with spline and center out of (0,0)
%15->object with 15 petals
%17->gear 7 petals
%18->kite rotacionada em 180o
%19->New figure(Kloeckner paper)
global type_obj
%global sp
%test for speed
global uu
global SP
global H
global NT
global Point_Center
global radius
% global TA

switch type_obj
    case 0
        x2=1*sin(t);
    case 1
        x2=1.5*sin(t);
    case 2
        x2=(2+0.3*cos(3*t)).*sin(t);
    case 3
        x2=sqrt(3*cos(t).*cos(t)+1).*sin(t);
    case 4
        %test for speed
        ind=find(uu==t);
        x2=SP(2,ind);
        %x2=[0 1]*Bsplineval(sp,t);
    case 5
        x2=4*sin(t);
    case 6
        x2=1+1*sin(t);
    case 7
        x2=2*sin(t);
    case 8
        x2=sin(t);
    case 9
        x2=2*(1+0.2*cos(9*t)).*sin(t);
    case 10
        TA=trig_int(t,NT);
        x2=transpose(TA*H).*sin(t);
    case 11
        TA=trig_int(t,NT);
        x2=Point_Center(2)+transpose(TA*H).*sin(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        x2=val_spline(t,H,NT).*sin(t);
    case 13
        x2=Point_Center(2)+val_spline(t,H,NT).*sin(t);
    case 15
        x2=2*(1+0.2*cos(15*t)).*sin(t);
    case 17
        x2=2*(1+0.2*cos(7*t)).*sin(t);
    case 18
        x2=-1.5*sin(t);
    case 19
        x2=sin(t-pi/4).*(1+sin(2*t)/2);
    case 20        
        a  = radius;
        x2 = ( 1 + a * cos(4*t) ).*sin(t);
end
