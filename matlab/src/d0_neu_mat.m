function [xmat] = d0_neu_mat(norder,h,srcinfo)
%
%  Representation:
%    u = D_{0} [\sigma] 
%
%  Data returned:
%    Neumann data (du/dn)
%
%
%  Input: 
%    norder - Alpert quadrature rule order
%    srcinfo(5,n) - source info
%       srcinfo(1:2,:) - locations
%       srcinfo(3:4,:) - normals
%       srcinfo(5,:) - dsdt
%
%    h - step size in parameter space
%
%  Output:
%    xmat - complex(n,n)
%      D_{0}'
%       

  [m,n] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (6,n)');
  D = specdiffmat(n,srcinfo);
  xmat2 = lap_dlp_mat(norder,h,srcinfo);
  xmat2 = xmat2 - 0.5*eye(n);
  xmat = -xmat2*D + 0j;
end
%  
%
%
%
%
