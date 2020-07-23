load('impedance_0005.mat')

figure; hold on;
h1 = plot(0,0);
for ik=1:length(khv)
    n_bd  = ceil(32*kh(end));
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 64)
        n_bd = 64;
    end
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
    var_imp = lambda_vecs(ik).coefs;
    N_imp = (length(var_imp)-1)/2;
    lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*var_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_bd',1:N_imp))*var_imp(N_imp+2:end)')';
    lambda_imp_orig = lambda_imp_f(t_bd);    
    h0 = plot(t_bd, lambda_imp_orig,'b');
    delete(h1)        
    h1 = plot(t_bd, lambda_imp,'r');
    err(ik) = norm(lambda_imp-lambda_imp_orig)*sqrt(pi/length(lambda_imp));
%     pause(1);
    
end

figure; semilogy(khv,err)