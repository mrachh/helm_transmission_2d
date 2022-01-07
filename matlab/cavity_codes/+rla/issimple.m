function out = issimple(xs,ys)
    nn = numel(xs);
    x_ft = fft(xs);
    y_ft = fft(ys);
    x_big_ft = zeros([5*nn,1]);
    y_big_ft = zeros([5*nn,1]);
    x_big_ft(1:(nn/2)) = x_ft(1:(nn/2));
    y_big_ft(1:(nn/2)) = y_ft(1:(nn/2));
    x_big_ft(5*nn -nn/2 +1:end) = x_ft((nn/2+1):end);
    y_big_ft(5*nn -nn/2 +1:end) = y_ft((nn/2+1):end);
    xbig = real(ifft(x_big_ft))*(5);
    ybig = real(ifft(y_big_ft))*(5);
    %plot(xs,ys,xbig,ybig)
    
    P = polyshape(xbig,ybig);
%     fprintf('NumRegions=%d\n',P.NumRegions)
    if P.NumRegions>1
        out = false;
    else
        out = true;
    end
    
    return
    