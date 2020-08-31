function [varargout] = resample_curve(xy,nb,eps,nmax)
%
%  This subroutine resamples a given curve whose xy 
%  and attempts to pass a bandlimited curve thorugh them
%
%  Input arguments:
%    xy(2,n) -
%      x,y coordinates of the input points
%    nb - 
%      parameter related to the bandlimit of the curve
%    eps - 
%      (optional) desired accuracy for interpolation of
%      output curve thorugh input points. Default value is
%      1e-12.
%    nmax - 
%      (optional) max number of iterations for attempting
%      to pass bandlimited curve with increasing 
%      storage and computation cost (both these parameters
%      double with with every iterate used). Default value
%      is 4.
%  
%  Output arguemnts:
%    srcinfo(6,nout) - 
%      resampled curve
%    h - 
%      spacing between points in parameter space
%    wsave -
%      (optional, if requested) an array which will
%      allow evaluation of x,y at arbitrary points
%    ts(n) -
%      (optional, if requested) coordinated in parameter
%      space closest to input points
%
%
  [~,n] = size(xy);

  if( nargout < 2 || nargout > 4)
    fprintf('invalid number of output arguments\n');
    fprintf('out arguments must be 2,3 or 4\n');
    vargout(1:nargout) = {0};
    return;
  end

  if (nargin == 2)
    epsuse = 1e-12;
    nmaxuse = 4;
  elseif (nargin == 3)
    epsuse = eps;
    nmaxuse = 4;
  elseif (nargin == 4)
    epsuse = eps;
    nmaxuse = nmax;
  else
    fprintf('invalid number of arguments\n');
  end
  nlarge = n;
  nout = n;
  lused = 0;
  lsave = 0;
  ierm = 0;
  mex_id_ = 'simple_curve_resampler_mem(i size_t[x], i double[xx], i size_t[x], i double[x], i size_t[x], io size_t[x], io size_t[x], io size_t[x], io size_t[x], io size_t[x])';
[nlarge, nout, lsave, lused, ierm] = curve_resampler(mex_id_, n, xy, nb, epsuse, nmaxuse, nlarge, nout, lsave, lused, ierm, 1, 2, n, 1, 1, 1, 1, 1, 1, 1, 1);
  if (ierm == 4)
    fprintf('band limit too small\n');
    fprintf('returning without computing anything\n');
    fprintf('run routine again with larger value of nb\n');
    varargout(1:nargout) = {0};
  else if (ierm == 2 && nargin >= 3)
    fprintf('interpolation error not to user specified accuracy\n');
    fprintf('try running the routine by setting nmax > 4\n');
  end

  srcinfo = zeros(6,nout);
  wsave = zeros(lsave,1);
  ts = zeros(n+1,1);
  nn = n+1;
  h = 0.0;
  curvelen = 0.0;
  ier = 0;
  mex_id_ = 'simple_curve_resampler_guru(i size_t[x], i double[xx], i size_t[x], i size_t[x], i size_t[x], i size_t[x], i size_t[x], io double[xx], io double[x], io double[x], io double[x], io double[x], io size_t[x])';
[srcinfo, h, curvelen, wsave, ts, ier] = curve_resampler(mex_id_, n, xy, nb, nlarge, lsave, lused, nout, srcinfo, h, curvelen, wsave, ts, ier, 1, 2, n, 1, 1, 1, 1, 1, 6, nout, 1, 1, lsave, nn, 1);
  varargout{1} = srcinfo;
  varargout{2} = h;

  if(nargout>=3)
    varargout{3} = wsave;
  end

  if(nargout == 4)
    varargout{4} = ts(1:n);
  end

end
%
%
%
%
%
