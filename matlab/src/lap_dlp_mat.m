function [xmat] = lap_dlp_mat(norder,h,srcinfo)
%
%  Representation:
%    u = D_{0} [\sigma] 
%
%  Data returned:
%    dirichlet data (du/dn)
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
%      D_{0}
%       

  [m,n] = size(srcinfo);
  assert(m==5,'srcinfo must be of shape (5,n)');
  xmat = complex(zeros(n),0);
  mex_id_ = 'lap_dlp_mat(i size_t[x], i size_t[x], i double[x], i double[xx], io dcomplex[xx])';
[xmat] = kern_mats(mex_id_, n, norder, h, srcinfo, xmat, 1, 1, 1, 5, n, n, n);
end
%  
%
%
%
%
%
%
