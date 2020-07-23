function [xmat] = rpcomb_neu_ext_mat(zpars,norder,h,srcinfo)
%
%  This subroutine generates a discretization of the representation
%  
%   u = S_{k} + i*\alpha*D_{k}(S_{ik})(\sigma) (Note this would work well
%     only for k real)
%  
%   correpsonding to the impendence boundary condition
%
%    du/dn = f
%
%  Input: 
%    zpars(3) - double complex
%      zpars(1) - Helmholtz paramter (k)
%      zpars(2) - \alpha in equation above
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
%       combined field matrix

  [m,n] = size(srcinfo);
  assert(m==5,'srcinfo must be of shape (5,n)');
  zk = zpars(1);
  alpha = zpars(2);
  spmat = sprime_ext_mat(zk,norder,h,srcinfo);
  zk2 = 1j*zk;
  spikmat = sprime_ext_mat(zk2,norder,h,srcinfo);
  spikmat = spikmat + 0.5*eye(n);
  sikmat = slp_mat(zk2,norder,h,srcinfo);
  zpars2 = complex(zeros(2,1));
  zpars2(1) = zk;
  zpars2(2) = 1j*zk;
  ddiffmat = ddiff_neu_mat(zpars2,norder,h,srcinfo);
  xmat=  spmat + 1j*alpha*(ddiffmat*sikmat ...
     +spikmat*spikmat-eye(n)/4);
end
%  
%
%
