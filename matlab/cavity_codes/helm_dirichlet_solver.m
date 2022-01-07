function [fields] = helm_dirichlet_solver(kh,src_info,sensor_info,cparams)
%HELM_DIRICHLET_SOLVER solves the helmholtz Dirichlet value problem
%for an analytic solution and a collection of plane waves for a rounded
%polygon
%
% Syntax: [u,varargout] = helm_dirichlet_solver(verts,zk,targs,angs,xyin,cparams)
% Input:
%   verts = xy coordinates of vertices describing the polygon (2,nverts)
%   zk - Helmholtz wave number
%   targs - xy coordinates of targets in exterior for evaluating potential
%   (2,nt)
%   angs - angles of incidence for the incoming plane waves (nangs,1)
%   xyin - location of interior point 
% Optional input:
%   cparams - options structure
%       cparams.rounded = true if corner rounding is to be used.
%                         false if no rounding is used (true)
%       cparams.autowidths = automatically compute widths (true)
%       cparams.autowidthsfac = if using autowidths, set widths
%                             to autowidthsfac*minimum of adjoining
%                             edges (0.4)
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
% output:
%   u - potential at targets (nt,nangs)
%   optional output arguments:
%   chnkr - discretized chunker used in solver
%   bd_sol - solution of integral equation corresponding to incident plane
%       waves (n,nangs)
%   F - compressed linear system corresponding to integral equation using
%       rskelf
%   err - estimated error in solving dirichlet problem
%

fields = [];
verts = [src_info.xs;src_info.ys];

if(nargin<=3)
    cparams_use = [];
    cparams_use.rounded = true;
    cparams_use.autowidths = true;
    cparams_use.autowidthsfac = 0.4;
    cparams_use.eps = 1e-6;
else
    cparams_use = cparams;
end
pref.k = 16;
chnkr = chunkerpoly(verts,cparams_use,pref);

ref_opts = []; ref_opts.maxchunklen = (2*pi)/abs(kh);
chnkr = refine(chnkr,ref_opts);

fkern = @(s,t) chnk.helm2d.kern(kh,s,t,'C',1);
dval = 0.5;
opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.rank_or_tol = 1e-10;
F = chunkerflam(chnkr,fkern,dval,opts_flam);


%%%%%%%%%%%%%%%%
   [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
   [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
   nt_uni = length(tgt_uni(:,1));
   induse = itgt + (idir-1)*nt_uni;
   
   t_dir_uni = t_dir_uni(:);
   x_dir = cos(t_dir_uni)';
   y_dir = sin(t_dir_uni)';
   n_dir = length(x_dir);
   %xs = src_info.xs(:)';
   %ys = src_info.ys(:)';
%%%%   
   rval = chnkr.r;
   rval = reshape(rval,2,chnkr.k*chnkr.nch);
   xs = rval(1,:);
   ys = rval(2,:);
   
fields.uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
tic
bd_sol = rskelf_sv(F,-fields.uinc);
toc
nt = length(tgt_uni);
nangs = length(t_dir_uni);
u = complex(zeros(nt,nangs));
%size(tgt_uni)
%sensor_info

eval_srcinfo = []; eval_srcinfo.r = chnkr.r(:,:,:); 
    eval_srcinfo.d = chnkr.d(:,:,:); eval_srcinfo.d2 = chnkr.d2(:,:,:);
eval_targinfo = []; eval_targinfo.r = tgt_uni';
kernmat = fkern(eval_srcinfo,eval_targinfo);

tic;
u = (kernmat*bd_sol);
toc;

%for i=1:nangs
    %size(chunkerkerneval(chnkr,fkern,bd_sol(:,i),tgt_uni'))
    %tic
    %opt_cnk = []; opt_cnk.forcesmooth = true;    
    %u(:,i) = chunkerkerneval(chnkr,fkern,bd_sol(:,i),tgt_uni',opt_cnk);
    %toc
%end
fields.uscat_tgt = u;

   
end