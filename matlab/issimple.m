function out = issimple(xs,ys)
    P = polyshape(xs,ys);
%     fprintf('NumRegions=%d\n',P.NumRegions)
    if P.NumRegions>1
        out = false;
    else
        out = true;
    end
    
    return
    