    N_bd          = 8;
    coefs_bd      = zeros(1,2*N_bd+1);
    coefs_bd(1)   = 1.;
    coefs_bd(4)   = 0.2;
    coefs_bd(5)   = 0.02;
    coefs_bd(7)   = 0.1;
    coefs_bd(9)   = 0.1;
    n_bd = 1000;
    if mod(n_bd,2)
        n_bd = n_bd+1;
    end
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
    h_bd     = 2*pi/n_bd;
    c_t   = cos(t_bd);
    s_t   = sin(t_bd);
    p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
        sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
    dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
        bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
    xs    = p_t .* c_t;
    ys    = p_t .* s_t;
    dxs   = dp_t .* c_t + p_t .* (-s_t);
    dys   = dp_t .* s_t + p_t .* c_t;
    ds    = sqrt(dxs.^2 + dys.^2);
    Lplane    = length(ds)*h_bd;
    d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';	
    d2xs = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
    d2ys = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);
    ds   = sqrt(dxs.^2 + dys.^2);
	H    = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
    bd_ref = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
        'dys',dys,'ds',ds,'Lplane',Lplane,'d2xs',d2xs,'d2ys',d2ys,'H',H,...
        'n_bd',n_bd,'N_bd',N_bd,'coefs_bd',coefs_bd);
