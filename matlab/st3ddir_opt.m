function [U, threads] = st3ddir_opt(srcinfo, targ, ifppregtarg)
%ST3DDIR_OPT Optimized direct Stokeslet sums, with any number of densities.
% U = st3ddir_opt(srcinfo,targ,ifppregtarg) has the Stokeslet inputs and
% output shapes of st3ddir. Default flag 3 preserves that interface; pass 1
% for optimized velocity only. Flags 2/3 use the original Fortran routine.
% Stresslets, rotlets and doublets are unsupported; use st3ddir for these.
% 'make matlab' builds its MEX gateway from the same library as fmm3d.
if nargin < 3, ifppregtarg = 3; end
validateattributes(ifppregtarg, {'double'}, {'scalar','integer','>=',0,'<=',3});
unsupported = {'strslet','strsvec','rotlet','rotvec','doublet','doubvec'};
assert(~any(isfield(srcinfo,unsupported)), 'st3ddir_opt:kernel', ...
    'st3ddir_opt supports Stokeslets only; use st3ddir for other kernels.');
nd = 1;
if isfield(srcinfo,'nd'), nd = srcinfo.nd; end
validateattributes(nd, {'double'}, {'scalar','integer','positive','finite'});
sources = srcinfo.sources;
validateattributes(sources, {'double'}, {'real','2d','nrows',3,'finite'});
validateattributes(targ, {'double'}, {'real','2d','nrows',3,'finite'});
ns = size(sources,2); nt = size(targ,2);
stoklet = srcinfo.stoklet;
validateattributes(stoklet, {'double'}, {'real','finite'});
if nd == 1
    assert(isequal(size(stoklet),[3,ns]), 'Stoklet must have shape [3,ns].');
else
    assert(size(stoklet,1)==nd && size(stoklet,2)==3 && ...
        size(stoklet,3)==ns && numel(stoklet)==3*nd*ns, ...
        'Stoklet must have shape [nd,3,ns].');
end
[pot,pre,grad,threads] = st3ddir_opt_mex(sources, ...
    reshape(stoklet,[3*nd,ns]),targ,ifppregtarg,nd);
U = struct('pottarg',[],'pretarg',[],'gradtarg',[]);
if ifppregtarg >= 1, U.pottarg = squeeze(reshape(pot,[nd,3,nt])); end
if ifppregtarg >= 2, U.pretarg = pre; end
if ifppregtarg >= 3, U.gradtarg = squeeze(reshape(grad,[nd,3,3,nt])); end
end
