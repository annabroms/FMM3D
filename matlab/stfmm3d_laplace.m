function [velocity, timing] = stfmm3d_laplace(sources, forces, targets, tolerance)
%STFMM3D_LAPLACE Target Stokeslet velocities via four Laplace charge densities.
% [velocity,timing] = stfmm3d_laplace(sources,forces,targets,tolerance)
% sources and forces: real finite double 3-by-Ns; targets: 3-by-Nt.
% velocity: 3-by-Nt, kernel (I/r + r*r'/r^3)/(8*pi), viscosity=1.
% tolerance defaults to 1e-9. No wall, source outputs or derivative outputs.
% This is an array-based interface, not the general stfmm3d interface.
% As with FMM3D, interactions below the backend domain-relative distance
% cutoff are omitted. Callers must supply nonsingular, resolved targets.
% Timing includes validation, packing, the lfmm3d call, and reconstruction.
% Call without arguments for a deterministic axial correctness check.
if nargin==0, self_test(); velocity=[]; timing=[]; return; end
total_timer=tic;
if nargin<4, tolerance=1e-9; end
validateattributes(sources,{'double'},{'real','finite','2d','nrows',3});
validateattributes(forces,{'double'},{'real','finite','2d','nrows',3,'ncols',size(sources,2)});
validateattributes(targets,{'double'},{'real','finite','2d','nrows',3});
validateattributes(tolerance,{'double'},{'real','finite','scalar','>',0,'<',1});
ns=size(sources,2); nt=size(targets,2);
timing=struct('backend_call_seconds',0,'packing_seconds',0,'assembly_seconds',0);
velocity=zeros(3,nt);
if ns==0 || nt==0
    timing.total_seconds=toc(total_timer); return;
end
if exist('lfmm3d','file')~=2 || exist('fmm3d','file')~=3
    error('stfmm3d_laplace:MissingDependency','Add the compiled FMM3D matlab folder to the path.');
end
timer=tic;
% Recenter all coordinates to avoid origin-dependent moment cancellation.
lower=min(min(sources,[],2),min(targets,[],2));
upper=max(max(sources,[],2),max(targets,[],2));
center=lower/2+upper/2;
Y=sources-center; X=targets-center;
src=struct('sources',Y,'nd',4, ...
    'charges',[forces/2; sum(Y.*forces,1)/2]);
timing.packing_seconds=toc(timer);
timer=tic;
[L,ier]=lfmm3d(tolerance,src,0,X,2);
timing.backend_call_seconds=toc(timer);
if ier~=0, error('stfmm3d_laplace:BackendFailure','Stokes Laplace backend returned ier=%d.',ier); end
timer=tic;
if numel(L.pottarg)~=4*nt || numel(L.gradtarg)~=12*nt
    error('stfmm3d_laplace:BackendFailure','Unexpected four-density Laplace output dimensions.');
end
potential=reshape(L.pottarg,4,nt);
gradient=reshape(L.gradtarg,4,3,nt);
if any(~isfinite(potential(:))) || any(~isfinite(gradient(:)))
    error('stfmm3d_laplace:BackendFailure','Nonfinite four-density Laplace output.');
end
for i=1:3
    velocity(i,:)=potential(i,:)+reshape(gradient(4,i,:),1,nt);
    for j=1:3
        velocity(i,:)=velocity(i,:)-X(j,:).*reshape(gradient(j,i,:),1,nt);
    end
end
if any(~isfinite(velocity(:)))
    error('stfmm3d_laplace:BackendFailure','Nonfinite reconstructed Stokes velocity.');
end
timing.assembly_seconds=toc(timer);
timing.total_seconds=toc(total_timer);
end

function self_test()
% Compare the four-potential representation against an explicit axial Stokeslet.
[u,t]=stfmm3d_laplace([0;0;0],[0;0;1],[0;0;2],1e-12);
assert(norm(u-[0;0;1/(8*pi)])<1e-13);
assert(t.total_seconds>=t.backend_call_seconds);
end
