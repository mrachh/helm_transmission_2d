function x1=coord1(t)
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
        x1=1*cos(t);
    case 1
        x1=cos(t)+0.65*cos(2*t)-0.65;
    case 2
        x1=(2+0.3*cos(3*t)).*cos(t);
    case 3
        x1=sqrt(3*cos(t).*cos(t)+1).*cos(t);
    case 4
        %test for speed           
        ind=find(uu==t);
        x1=SP(1,ind);
        %x1=[1 0]*Bsplineval(sp,t);
    case 5
        x1=3*cos(t);
    case 6
        x1=1+1*cos(t);
    case 7
        x1=2*cos(t);
    case 8
        x1=20*cos(t);
    case 9
        x1=2*(1+0.2*cos(9*t)).*cos(t);
    case 10
        %Here I use H fr the coefs of the trig pol
        %and NT for the number f coeficients.
        TA=trig_int(t,NT);
        x1=transpose(TA*H).*cos(t);
    case 11
        TA=trig_int(t,NT);
        x1=Point_Center(1)+transpose(TA*H).*cos(t);
    case 12
        %Here I use H for the coefs of the b-spline
        %and NT for the order of the spline
        x1=val_spline(t,H,NT).*cos(t);
    case 13
        x1=Point_Center(1)+val_spline(t,H,NT).*cos(t);
    case 15
        x1=2*(1+0.2*cos(15*t)).*cos(t);
    case 17
        x1=2*(1+0.2*cos(7*t)).*cos(t);
    case 18
        x1=0.65-0.65*cos(2*t)-cos(t);
    case 19
        x1=(3/4)*cos(t-pi/4).*(1+sin(2*t)/2);
    case 20
        a  = radius;
        x1 = ( 1 + a * cos(4*t) ).*cos(t);
        
end
