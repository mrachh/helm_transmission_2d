xtrue = src_info.xs;
ytrue = src_info.ys;
htrue = src_info.h;
dxtrue= src_info.dxs;
dytrue= src_info.dys;
strue = src_info.ds*htrue;

xfit = src_info_out.xs;
yfit = src_info_out.ys;
dxfit= src_info_out.dxs;
dyfit= src_info_out.dys;
hfit = src_info_out.h;
sfit = src_info_out.ds*hfit;


if (1<0)
ks = -10:0.1:10;
[KX,KY,XX] = meshgrid(ks,ks,xtrue);
[~,~,YY] = meshgrid(ks,ks,ytrue);
[~,~,SS] = meshgrid(ks,ks,strue);

ttrue = sum(exp(1i*(KX.*XX+KY.*YY)).*SS,3);

ks = -10:0.1:10;
[KX,KY,XX] = meshgrid(ks,ks,xfit);
[~,~,YY] = meshgrid(ks,ks,yfit);
[~,~,SS] = meshgrid(ks,ks,sfit);

tfit = sum(exp(1i*(KX.*XX+KY.*YY)).*SS,3);

end

ks = -100:0.5:100;
[KX,KY,XX] = meshgrid(ks,ks,xtrue);
[~,~,YY] = meshgrid(ks,ks,ytrue);
[~,~,DX] = meshgrid(ks,ks,dxtrue);
[~,~,DY] = meshgrid(ks,ks,dytrue);
[~,~,SS] = meshgrid(ks,ks,strue);

DR = sqrt(DX.^2+DY.^2);
DX = DX./DR;
DY = DY./DR;
DKR = (KX.^2+KY.^2);
ttrue = sum((DY.*KX-DX.*KY)./DKR.*exp(1i*(KX.*XX+KY.*YY)).*SS,3);

ks = -100:0.5:100;
[KX,KY,XX] = meshgrid(ks,ks,xfit);
[~,~,YY] = meshgrid(ks,ks,yfit);
[~,~,SS] = meshgrid(ks,ks,sfit);
[~,~,DX] = meshgrid(ks,ks,dxfit);
[~,~,DY] = meshgrid(ks,ks,dyfit);

DR = sqrt(DX.^2+DY.^2);
DX = DX./DR;
DY = DY./DR;
DKR = (KX.^2+KY.^2);
tfit = sum((DY.*KX-DX.*KY)./DKR.*exp(1i*(KX.*XX+KY.*YY)).*SS,3);

[KKX,KKY] = meshgrid(ks,ks);

