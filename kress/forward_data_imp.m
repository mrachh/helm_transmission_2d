function [umeas,ubd,dubd,inv_Fw] = forward_data_imp(kh,bd_data,src,tgt,t_bd,lambda_imp)

    %boundary
    n_bd = length(t_bd);
    
    %generating operators
    S  = slmat(kh,src,t_bd);
    D  = dlmat(kh,src,t_bd);
    Sp = sprimelmat(kh,src,t_bd);
    T  = dprimelmat(kh,src,t_bd);
    
    %solving the system
    eta    = kh;    
    Fw_mat = (T + 1i* eta * (Sp - eye(n_bd)/2) + ... %du part
        1i * kh * bsxfun(@times,lambda_imp',D+eye(n_bd)/2+1i*eta*S));%i lambda u part  
    inv_Fw = inv(Fw_mat);
    pot    = inv_Fw * bd_data;

    %calculating scattered field at target
    S_tgt = slmat_out(kh,src,tgt);
    D_tgt = dlmat_out(kh,src,tgt);
    umeas = (D_tgt + 1i * eta * S_tgt)*pot;
    
    %calculating u
    ubd  = (D + eye(n_bd)/2 + 1i* eta *S ) * pot;
    dubd = (T + 1i * eta * ( Sp - eye(n_bd)/2 ) ) * pot;
    
    return
