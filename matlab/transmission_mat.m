function [xmat] = transmission_mat(zks,a,b,norder,h,srcinfo)
%
%
%  Input: 
%    zks(2) - complex
%       Helmholtz parameters
%    a(2) - complex
%       scaling for jump in u
%    b(2) - complex
%      scaling for jump in dudn
%    norder - Alpert quadrature rule order
%    srcinfo(5,n) - source info
%       srcinfo(1:2,:) - locations
%       srcinfo(3:4,:) - normals
%       srcinfo(5,:) - dsdt
%
%    h - step size in parameter space
%
%  Output:
%    xmat - complex(2*n,2*n)
%       transmission matrix

  [m,n] = size(srcinfo);
  assert(m==5,'srcinfo must be of shape (5,n)');
  xmat = complex(zeros(2*n),0);
  nsys = 2*n;
  mex_id_ = 'trans_mat(i int[x], i int[x], i double[x], i double[xx], i dcomplex[x], i dcomplex[x], i dcomplex[x], io dcomplex[xx])';
[xmat] = kern_mats(mex_id_, n, norder, h, srcinfo, zks, a, b, xmat, 1, 1, 1, 5, n, 2, 2, 2, nsys, nsys);
end
  
