function [D] = specdiffmat(N,srcinfo)
%
%  Note that the routine expects N to be even
%  Do not run with odd N
%
    tj = 2*pi/N*(0:N-1);
    D = (-1).^(1:N) .* cot(tj/2) / 2;   % note overall - sgn due to 1st row not col
    D(1) = 0;              % kill the Inf
    D = circulant(D);
    D = diag(1.0./(srcinfo(5,:)))*D;
end
