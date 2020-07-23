function lambda_imp = lambda_imp_f(t)

lambda_imp = zeros(1,length(t));
% lambda_imp(t<=pi/2) = 0.5;
% lambda_imp((t>pi/2) & (t<3*pi/2)) = 0.6;
% lambda_imp(t>=3*pi/2) = 0.5;

lambda_imp(t<=pi) = t(t<=pi)/pi;
lambda_imp(t>pi) = - t(t>pi)/pi + 2;

    