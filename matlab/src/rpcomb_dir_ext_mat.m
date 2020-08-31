function [xmat] = rpcomb_dir_ext_mat(zpars,norder,h,srcinfo)
%
%  This subroutine generates a discretization of the representation
%  
%   u = S_{k} + i*\alpha*D_{k}(S_{ik})(\sigma) (Note this would work well
%     only for k real)
%  
%   correpsonding to the impendence boundary condition
%
%    u = f
%
%  Input: 
%    zpars(2) - double complex
%      zpars(1) - Helmholtz paramter (k)
%      zpars(2) - \alpha in equation above
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

  [m,~] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (5,n)');
  zk = zpars(1);
  smat = slp_mat(zk,norder,h,srcinfo);
  zk2 = 1j*abs(zk);
  sikmat = slp_mat(zk2,norder,h,srcinfo);
  dmat = dlp_ext_mat(zk,norder,h,srcinfo);
  xmat=  smat+1j*zpars(2)*dmat*sikmat;
end
%  
%
%
