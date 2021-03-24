function lambda_imp_gau = lambda_imp_gau_f(t,rd)
    tt = t(:);
    n = ceil(10/rd);
    if(mod(n,2) == 0) 
        n = n+1;
    end
    nind = -n:1:n;
    coefs(1:2:2*n+1) = 4.0./(nind(1:2:2*n+1).^2)/2/pi;
    coefs(n+1) = 3*pi/2;
    coefs = coefs.*exp(-nind.^2*rd^2)/pi;
    z = exp(1j*tt*nind);
    disp(size(z))
    disp(size(coefs))
    lambda_imp_gau = real(z*coefs'); 
    
    lambda_imp_gau = reshape(lambda_imp_gau,size(t));
end