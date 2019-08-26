function [srcinfo,h] = ellipse(a,b,n)
    srcinfo = zeros(5,n);
    ts = (1:n)/n*2*pi;
    srcinfo(1,:) = a*cos(ts);
    srcinfo(2,:) = b*sin(ts);
    srcinfo(5,:) = sqrt((b*cos(ts)).^2 + (a*sin(ts)).^2);
    srcinfo(3,:) = b*cos(ts)./srcinfo(5,:);
    srcinfo(4,:) = a*sin(ts)./srcinfo(5,:);
    h = 2*pi/n;
end
