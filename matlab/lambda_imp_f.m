function lambda_imp = lambda_imp_f(t)

lambda_imp = zeros(1,length(t));
%function 1
% lambda_imp(t<=pi/2) = 0.5;
% lambda_imp((t>pi/2) & (t<3*pi/2)) = 0.6;
% lambda_imp(t>=3*pi/2) = 0.5;

%function 2
% lambda_imp(t<=pi) = t(t<=pi)/pi+0.75;
% lambda_imp(t>pi) = - t(t>pi)/pi + 2+0.75;

%function 3
% lambda_imp(t<=pi/2) = .3;
% 
% lambda_imp((t>pi/2) & (t<=pi)) = 0.5/(pi/2) * t((t>pi/2) & (t<=pi)) + (0.3-0.5/(pi/2)*pi/2);
% 
% lambda_imp((t>pi) & (t<=3*pi/2)) = (-0.5/(pi/2)) * t((t>pi) & (t<=3*pi/2)) + (0.8-(-0.5/(pi/2))*pi);
% 
% lambda_imp(t>3*pi/2) = .3;

%function 4
% lambda_imp = 1 + 0.1 * cos(t);
% lambda_imp = 0.7;

% function paper
lambda_imp(t<=pi) = -0.1*t(t<=pi)/pi + 0.6;
lambda_imp(t>pi) = 0.1*t(t>pi)/pi + 0.4;

    