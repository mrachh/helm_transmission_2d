function [xmat] = rpcomb_imp_ext_mat(zpars,norder,h,srcinfo)
%
%  This subroutine generates a discretization of the representation
%  
%   u = S_{k} + i*\alpha*D_{k}(S_{ik})(\sigma) (Note this would work well
%     only for k real)
%  
%   correpsonding to the impendence boundary condition
%
%    u + ik \eta du/dn = f
%
%  Input: 
%    zpars(3) - double complex
%      zpars(1) - Helmholtz paramter (k)
%      zpars(2) - \alpha in equation above
%      zpars(3) - \eta the impendence strength
%    norder - Alpert quadrature rule order
%    srcinfo(5,n) - source info
%       srcinfo(1:2,:) - locations
%       srcinfo(3:4,:) - normals
%       srcinfo(5,:) - dsdt
%       srcinfo(6,:) - curvature
%
%    h - step size in parameter space
%
%  Output:
%    xmat - complex(n,n)
%       combined field matrix

  [m,n] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (6,n)');
  zk = zpars(1);
  eta = zpars(3);
  alpha = zpars(2);
  smat = slp_mat(zk,norder,h,srcinfo);
  zk2 = 1j*abs(zk);
  sikmat = slp_mat(zk2,norder,h,srcinfo);
  dmat = dlp_ext_mat(zk,norder,h,srcinfo);
  xmat=  smat+1j*alpha*dmat*sikmat;
  
  spmat = sprime_ext_mat(zk,norder,h,srcinfo);
  spikmat = sprime_ext_mat(zk2,norder,h,srcinfo);
  spikmat = spikmat + 0.5*eye(n);
  zpars2 = complex(zeros(2,1));
  zpars2(1) = zk;
  zpars2(2) = 1j*abs(zk);
  ddiffmat = ddiff_neu_mat(zpars2,norder,h,srcinfo);
  xmat = xmat+1j*zk*eta*(spmat + 1j*alpha*(ddiffmat*sikmat ...
     +spikmat*spikmat-eye(n)/4));
  
end
%  
%
%
