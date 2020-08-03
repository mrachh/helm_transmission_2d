function [der1,der2] = spectral_derivatives(input)

N = length(input);
input_hat = fft(input);
der_func1 = 1i*[0:N/2-1 0 -N/2+1:-1] .* input_hat;
der_func2 = 1i*[0:N/2-1 0 -N/2+1:-1] .* der_func1;
der1 = real(ifft(der_func1));
der2 = real(ifft(der_func2));
