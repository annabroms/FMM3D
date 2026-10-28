function [u, timing] = stfmm3d_laplace(sources, forces, targ, tol)
% STFMM3D_LAPLACE    Stokeslet target velocities using the Laplace FMM.
%
% [u, timing] = stfmm3d_laplace(sources, forces, targ, tol)
%
% Computes u(x) = sum_j [f_j/r + R*(R'*f_j)/r^3]/(8*pi),
% where R = x-y_j and r = norm(R). The viscosity is assumed to be 1.
%
% Inputs (finite real doubles):
%   sources  - source locations, 3-by-ns
%   forces   - Stokeslet strengths, 3-by-ns
%   targ     - target locations, 3-by-nt
%   tol      - requested FMM precision (optional, default 1e-9)
%
% Outputs:
%   u        - velocity at targets, 3-by-nt
%   timing   - packing, backend call, assembly and total times in seconds
%
% Uses four Laplace charge densities and target potentials/gradients only.
% For pressure, velocity gradients or other Stokes kernels, use STFMM3D.
% Interactions below the FMM3D distance cutoff are omitted, as in LFMM3D.
% Call without inputs for a small correctness check.
%
% See also: STFMM3D, LFMM3D, STFMM3D_LAPLACEPERFTEST

if nargin == 0
    self_test();
    u = []; timing = [];
    return
end
total_timer = tic;
if nargin < 4, tol = 1e-9; end
validateattributes(sources, {'double'}, {'real','finite','2d','nrows',3});
validateattributes(forces, {'double'}, {'real','finite','size',size(sources)});
validateattributes(targ, {'double'}, {'real','finite','2d','nrows',3});
validateattributes(tol, {'double'}, {'real','finite','scalar','>',0,'<',1});

nt = size(targ,2);
u = zeros(3,nt);
timing = struct('packing_seconds',0, 'backend_call_seconds',0, 'assembly_seconds',0);
if isempty(sources) || isempty(targ)
    timing.total_seconds = toc(total_timer);
    return
end

% Shift the common origin to reduce cancellation in the fourth potential.
t = tic;
lower = min(min(sources,[],2), min(targ,[],2));
upper = max(max(sources,[],2), max(targ,[],2));
center = lower/2 + upper/2;
y = sources - center;
x = targ - center;

% psi_1,...,psi_3 have charges f/2; psi_4 has charges (y dot f)/2.
% LFMM3D includes the Laplace normalization 1/(4*pi).
srcinfo.sources = y;
srcinfo.nd = 4;
srcinfo.charges = [forces/2; sum(y.*forces,1)/2];
timing.packing_seconds = toc(t);
t = tic;
[L, ier] = lfmm3d(tol,srcinfo,0,x,2);  % no source fields; target pot and grad
timing.backend_call_seconds = toc(t);
if ier ~= 0
    error('stfmm3d_laplace:BackendFailure','LFMM3D returned ier=%d.',ier);
end

% grad(j,i,m) is d_i psi_j at target m. Reshape also handles nt=1.
t = tic;
if numel(L.pottarg) ~= 4*nt || numel(L.gradtarg) ~= 12*nt
    error('stfmm3d_laplace:BackendFailure','Unexpected LFMM3D output dimensions.');
end
pot = reshape(L.pottarg,4,nt);
grad = reshape(L.gradtarg,4,3,nt);
if any(~isfinite(pot(:))) || any(~isfinite(grad(:)))
    error('stfmm3d_laplace:BackendFailure','Nonfinite LFMM3D output.');
end

% u_i = psi_i - sum_j x_j d_i psi_j + d_i psi_4.
for i = 1:3
    u(i,:) = pot(i,:) + reshape(grad(4,i,:),1,nt);
    for j = 1:3
        u(i,:) = u(i,:) - x(j,:).*reshape(grad(j,i,:),1,nt);
    end
end
if any(~isfinite(u(:)))
    error('stfmm3d_laplace:BackendFailure','Nonfinite Stokes velocity.');
end
timing.assembly_seconds = toc(t);
timing.total_seconds = toc(total_timer);
end

function self_test()
% Check the normalization against an axial Stokeslet at distance 2.
[u,timing] = stfmm3d_laplace([0;0;0],[0;0;1],[0;0;2],1e-12);
assert(norm(u-[0;0;1/(8*pi)]) < 1e-13);
assert(timing.total_seconds >= timing.backend_call_seconds);
end
