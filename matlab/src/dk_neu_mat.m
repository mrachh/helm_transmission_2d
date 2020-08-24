function [xmat] = dk_neu_mat(zpars,norder,h,srcinfo,d0)
%
%  Representation:
%    u = D_{0} [\sigma] 
%
%  Data returned:
%    Neumann data (du/dn)
%
%
%  Input: 
%    zpars(1) - 
%      zpars(1) - Helmholtz paramter (k1)
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
  assert(m==5,'srcinfo must be of shape (5,n)');
  if (nargin == 4)
    d0mat = d0_neu_mat(norder,h,srcinfo);
  elseif (nargin == 5)
    d0mat = d0;
  else
     fprintf('invalid number of arguments');
     return;
  end        
  dkdiffmat = ddiff0_neu_mat(zpars,norder,h,srcinfo);
  xmat = dkdiffmat - d0mat;
end
%  
%
%
%
%
