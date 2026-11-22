function [U, threads] = st3ddir_opt(srcinfo, targ, ifppregtarg)
%ST3DDIR_OPT Optimized direct Stokeslet/stresslet sums, any number of densities.
% U = st3ddir_opt(srcinfo,targ,ifppregtarg) has the Stokeslet and (type I)
% stresslet inputs and the output shapes of st3ddir. Give srcinfo.stoklet,
% srcinfo.strslet with srcinfo.strsvec, or both (contributions are added).
% Default flag 3 preserves the st3ddir interface; pass 1 for optimized
% velocity only. Flags 2/3 use the original Fortran routines.
% Rotlets and doublets are unsupported; use st3ddir for these.
% 'make matlab' builds its MEX gateway from the same library as fmm3d.
if nargin < 3, ifppregtarg = 3; end
validateattributes(ifppregtarg, {'double'}, {'scalar','integer','>=',0,'<=',3});
unsupported = {'rotlet','rotvec','doublet','doubvec'};
assert(~any(isfield(srcinfo,unsupported)), 'st3ddir_opt:kernel', ...
    'st3ddir_opt supports Stokeslets and stresslets only; use st3ddir for other kernels.');
ifstok = isfield(srcinfo,'stoklet');
ifstrs = isfield(srcinfo,'strslet');
assert(ifstrs == isfield(srcinfo,'strsvec'), 'st3ddir_opt:kernel', ...
    'Stresslets require both strslet and strsvec.');
assert(ifstok || ifstrs, 'st3ddir_opt:kernel', ...
    'Give stoklet and/or strslet with strsvec.');
nd = 1;
if isfield(srcinfo,'nd'), nd = srcinfo.nd; end
validateattributes(nd, {'double'}, {'scalar','integer','positive','finite'});
sources = srcinfo.sources;
validateattributes(sources, {'double'}, {'real','2d','nrows',3,'finite'});
validateattributes(targ, {'double'}, {'real','2d','nrows',3,'finite'});
ns = size(sources,2); nt = size(targ,2);
stoklet = []; strslet = []; strsvec = [];
if ifstok, stoklet = packed(srcinfo.stoklet, 'Stoklet', nd, ns); end
if ifstrs
    strslet = packed(srcinfo.strslet, 'Strslet', nd, ns);
    strsvec = packed(srcinfo.strsvec, 'Strsvec', nd, ns);
end
[pot,pre,grad,threads] = st3ddir_opt_mex(sources, stoklet, strslet, ...
    strsvec, targ, ifppregtarg, nd);
U = struct('pottarg',[],'pretarg',[],'gradtarg',[]);
if ifppregtarg >= 1, U.pottarg = squeeze(reshape(pot,[nd,3,nt])); end
if ifppregtarg >= 2, U.pretarg = pre; end
if ifppregtarg >= 3, U.gradtarg = squeeze(reshape(grad,[nd,3,3,nt])); end
end

function a = packed(a, name, nd, ns)
% Check a strength array of shape [3,ns] (nd=1) or [nd,3,ns] and return it
% as [3*nd,ns] for the MEX gateway.
validateattributes(a, {'double'}, {'real','finite'});
if nd == 1
    assert(isequal(size(a),[3,ns]), '%s must have shape [3,ns].', name);
else
    assert(size(a,1)==nd && size(a,2)==3 && size(a,3)==ns && ...
        numel(a)==3*nd*ns, '%s must have shape [nd,3,ns].', name);
end
a = reshape(a,[3*nd,ns]);
end
