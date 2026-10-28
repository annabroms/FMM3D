function [U, timing] = stfmm3d_laplace(eps,srcinfo,ifppreg,targ,ifppregtarg)
% STFMM3D_LAPLACE    FMM in 3D for Stokeslet velocities via Laplace kernels.
%
% U = stfmm3d_laplace(eps,srcinfo,ifppreg,targ,ifppregtarg)
% [U,timing] = stfmm3d_laplace(eps,srcinfo,ifppreg,targ,ifppregtarg)
%
%  Stokeslet FMM in R^3: evaluate all pairwise particle interactions
%  (ignoring self-interactions) and, possibly, interactions with targets,
%  using the Laplace fast multipole method with requested precision eps.
%
%  This routine computes the sum for the velocity vector,
%
%       u_i(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
%
%  for each requested evaluation point x and i=1,2,3, where sigma^{(m)}
%  is the Stokeslet force (a 3 vector per source point y^{(m)}).
%  Repeated indices are summed over 1,2,3 (Einstein convention).
%  The dynamic viscosity is assumed to be 1, as in STFMM3D.
%
%  For a source y and target x, let r_i = x_i-y_i (note sign)
%  and let r = sqrt(r_1^2 + r_2^2 + r_3^2). The Stokeslet is
%
%       G_{ij}(x,y) = (delta_{ij}/r + r_i r_j/r^3)/(8*pi).
%
%  Relation to Laplace kernels:
%  Let K(x,y) = 1/(4*pi*r), the charge kernel used by LFMM3D.
%  Since d_i K = -r_i/(4*pi*r^3), where d_i differentiates with
%  respect to the evaluation coordinate x_i, the Stokeslet satisfies
%
%       G_{ij}(x,y) = [delta_{ij} K - (x_j-y_j) d_i K]/2.
%
%  For each force density, form four scalar Laplace potentials:
%
%       psi_j(x) = sum_m [sigma^{(m)}_j/2] K(x,y^{(m)}),  j=1,2,3,
%       psi_4(x) = sum_m [(y^{(m)} dot sigma^{(m)})/2] K(x,y^{(m)}).
%
%  Thus the four charge densities passed to LFMM3D are sigma_1/2,
%  sigma_2/2, sigma_3/2 and (y dot sigma)/2. Combining the returned
%  potentials and their evaluation-coordinate gradients gives
%
%       u_i(x) = psi_i(x) - sum_{j=1}^3 x_j d_i psi_j(x) + d_i psi_4(x).
%
%  In particular, the last two terms combine to
%
%       sum_m R_i^{(m)} (R^{(m)} dot sigma^{(m)})/(8*pi*r_m^3),
%
%  where R^{(m)} = x-y^{(m)} and r_m = norm(R^{(m)}). The potential
%  psi_i supplies the remaining sum_m sigma^{(m)}_i/(8*pi*r_m) term
%  of the Stokeslet velocity.
%  The factors of 1/2 in the charges convert the Laplace normalization
%  1/(4*pi) to the Stokeslet normalization 1/(8*pi).
%
%  The code uses a common shifted origin c: x is replaced by x-c and
%  y by y-c, including in the fourth charge density and the assembly
%  formula. The displacement is unchanged and the c-dependent terms
%  cancel exactly. This shift reduces cancellation in the assembly.
%  In particular, translating all points by a leaves psi_1,...,psi_3
%  unchanged but adds a dot psi to psi_4. Its gradient must cancel the
%  extra sum_j a_j d_i psi_j term in the velocity reconstruction.
%  Large absolute coordinates therefore produce large intermediate terms
%  whose errors can dominate the final velocity. Recentering is unnecessary
%  in exact arithmetic, but protects against this loss of accuracy.
%
%  Evaluation points can be the sources (ifppreg = 1), other targets
%  (ifppregtarg = 1), or both. When x=y^{(m)}, the corresponding term
%  is dropped. Interactions below the FMM3D distance cutoff are also
%  omitted, as in LFMM3D.
%
%  Only Stokeslet velocities are computed. Pressure, velocity gradients,
%  stresslets, rotlets and doublets are not supported; requests for these
%  raise an error. Use STFMM3D for those quantities and kernels.
%  Internally, four Laplace charge densities per force density and their
%  potentials/gradients are used. These gradients are intermediate
%  Laplace fields, not returned gradients of the Stokes velocity.
%
%  Args:
%
%  -  eps: double
%        precision requested from LFMM3D, 0 < eps < 1
%  -  srcinfo: structure containing the following source information:
%     *  srcinfo.sources: double(3,ns)
%           source locations, y^{(m)}
%     *  srcinfo.nd: integer
%           number of force densities for the same points (optional,
%           default nd = 1)
%     *  srcinfo.stoklet: double(nd,3,ns), or double(3,ns) for nd = 1
%           Stokeslet force strengths, sigma^{(m)} (optional;
%           if omitted, velocities are zero)
%        Coordinates and force strengths must be finite real doubles.
%  -  ifppreg: integer
%        source evaluation flag: 0 for no source output, 1 for velocity
%        (flags 2 and 3 for pressure/gradient are not supported)
%  -  targ: double(3,nt)
%        finite real target locations (optional)
%  -  ifppregtarg: integer
%        target evaluation flag: 0 for no target output, 1 for velocity
%        (optional, default 0; flags 2 and 3 are not supported)
%
%  Returns (same velocity shapes and field names as STFMM3D):
%
%  -  U.pot: velocity at source locations, if requested
%  -  U.pottarg: velocity at target locations, if requested
%        Velocities have shape (nd,3,n), with singleton dimensions
%        removed as in STFMM3D; n is ns or nt, respectively.
%  -  U.pre, U.grad, U.pretarg, U.gradtarg: always empty
%        Unrequested velocity fields are also empty.
%  -  timing: optional structure with elapsed times in seconds
%        Timing is collected only when a second output is requested.
%     *  timing.packing_seconds: coordinate shifting and charge packing
%     *  timing.backend_call_seconds: LFMM3D call
%     *  timing.assembly_seconds: velocity assembly from Laplace fields
%     *  timing.total_seconds: full call, including input validation
%
%  Example (same call and velocity outputs as STFMM3D):
%     srcinfo = struct('sources',sources,'stoklet',forces);
%     U = stfmm3d_laplace(1e-9,srcinfo,0,targ,1);
%     [U,timing] = stfmm3d_laplace(1e-9,srcinfo,0,targ,1);
%
%  Call without inputs for a small correctness check.
%
% See also: STFMM3D, LFMM3D, STFMM3D_LAPLACEPERFTEST

if nargin == 0
    self_test();
    U = []; timing = [];
    return
end
collect_timing = nargout > 1;
timing = [];
if collect_timing
    total_timer = tic;
    timing = struct('packing_seconds',0, 'backend_call_seconds',0, ...
        'assembly_seconds',0, 'total_seconds',0);
end
narginchk(3,5);
if nargin < 4, targ = zeros(3,0); end
if nargin < 5, ifppregtarg = 0; end
validateattributes(eps, {'double'}, {'real','finite','scalar','>',0,'<',1});
validateattributes(srcinfo, {'struct'}, {'scalar'});
validateattributes(ifppreg, {'numeric'}, {'real','finite','scalar','integer','>=',0,'<=',3});
validateattributes(ifppregtarg, {'numeric'}, {'real','finite','scalar','integer','>=',0,'<=',3});
if ifppreg > 1 || ifppregtarg > 1
    error('stfmm3d_laplace:UnsupportedOutput', ...
        'Only velocity flags 0 and 1 are supported. Use STFMM3D for pressure or gradients.');
end
unsupported = {'strslet','strsvec','rotlet','rotvec','doublet','doubvec'};
if any(isfield(srcinfo,unsupported))
    error('stfmm3d_laplace:UnsupportedKernel', ...
        'Only Stokeslet forces are supported. Use STFMM3D for other Stokes kernels.');
end
sources = srcinfo.sources;
validateattributes(sources, {'double'}, {'real','finite','2d','nrows',3});
validateattributes(targ, {'double'}, {'real','finite','2d','nrows',3});
ns = size(sources,2);
nt = size(targ,2);
nd = 1;
if isfield(srcinfo,'nd'), nd = srcinfo.nd; end
validateattributes(nd, {'numeric'}, {'real','finite','scalar','integer','positive'});
if isfield(srcinfo,'stoklet')
    forces = srcinfo.stoklet;
    validateattributes(forces, {'double'}, {'real','finite'});
    if nd == 1
        valid_shape = isequal(size(forces),[3,ns]) || ...
            (ndims(forces) <= 3 && size(forces,1) == 1 && ...
            size(forces,2) == 3 && size(forces,3) == ns);
    else
        valid_shape = ndims(forces) <= 3 && size(forces,1) == nd && ...
            size(forces,2) == 3 && size(forces,3) == ns;
    end
    if ~valid_shape
        error('stfmm3d_laplace:InvalidStrengths', ...
            'stoklet must have shape (nd,3,ns), or (3,ns) for nd = 1.');
    end
    forces = reshape(forces,nd,3,ns);
else
    forces = zeros(nd,3,ns);
end

U = struct('pot',[], 'pre',[], 'grad',[], ...
    'pottarg',[], 'pretarg',[], 'gradtarg',[]);
if ns == 0 || (ifppreg == 0 && (ifppregtarg == 0 || nt == 0)) || ...
        ~isfield(srcinfo,'stoklet')
    if ifppreg == 1, U.pot = squeeze(zeros(nd,3,ns)); end
    if ifppregtarg == 1, U.pottarg = squeeze(zeros(nd,3,nt)); end
    if collect_timing, timing.total_seconds = toc(total_timer); end
    return
end
if ifppregtarg == 1 && nt == 0, U.pottarg = squeeze(zeros(nd,3,0)); end

% Shift the origin to reduce cancellation between x_j*d_i psi_j and d_i psi_4.
if collect_timing, t = tic; end
lower = min(sources,[],2);
upper = max(sources,[],2);
if ifppregtarg == 1 && nt > 0
    lower = min(lower,min(targ,[],2));
    upper = max(upper,max(targ,[],2));
end
center = lower/2 + upper/2;
y = sources - center;
if ifppregtarg == 1
    x = targ - center;
else
    x = zeros(3,0);
end

% psi_1,...,psi_3 have charges f/2; psi_4 has charges (y dot f)/2.
% Density index varies fastest. LFMM3D includes the factor 1/(4*pi).
lapinfo.sources = y;
lapinfo.nd = 4*nd;
moment = sum(reshape(y,1,3,ns).*forces,2);
lapinfo.charges = reshape(cat(2,forces/2,moment/2),4*nd,ns);
if collect_timing
    timing.packing_seconds = toc(t);
    t = tic;
end
if ifppregtarg == 1 && nt > 0
    [L, ier] = lfmm3d(eps,lapinfo,2*ifppreg,x,2);
else
    [L, ier] = lfmm3d(eps,lapinfo,2*ifppreg);
end
if collect_timing, timing.backend_call_seconds = toc(t); end
if ier ~= 0
    error('stfmm3d_laplace:BackendFailure','LFMM3D returned ier=%d.',ier);
end

if collect_timing, t = tic; end
if ifppreg == 1
    U.pot = assemble_velocity(L.pot,L.grad,y,nd);
end
if ifppregtarg == 1 && nt > 0
    U.pottarg = assemble_velocity(L.pottarg,L.gradtarg,x,nd);
end
if collect_timing
    timing.assembly_seconds = toc(t);
    timing.total_seconds = toc(total_timer);
end
end

function u = assemble_velocity(pot,grad,x,nd)
% grad(d,j,i,m) is d_i psi_j for force density d at evaluation point m.
n = size(x,2);
if numel(pot) ~= 4*nd*n || numel(grad) ~= 12*nd*n || ...
        any(~isfinite(pot(:))) || any(~isfinite(grad(:)))
    error('stfmm3d_laplace:BackendFailure','Invalid LFMM3D output.');
end
% u_i = psi_i - sum_j x_j d_i psi_j + d_i psi_4.
% Assemble each component once, avoiding repeated output updates and slices.
if nd == 1
    pot = reshape(pot,4,n);
    grad = reshape(grad,12,n);
    u = zeros(3,n);
    for i = 1:3
        k = 4*(i-1);
        u(i,:) = pot(i,:) + grad(k+4,:) ...
            - x(1,:).*grad(k+1,:) - x(2,:).*grad(k+2,:) - x(3,:).*grad(k+3,:);
    end
else
    pot = reshape(pot,nd,4,n);
    grad = reshape(grad,nd,12,n);
    x1 = reshape(x(1,:),1,1,n);
    x2 = reshape(x(2,:),1,1,n);
    x3 = reshape(x(3,:),1,1,n);
    u = zeros(nd,3,n);
    for i = 1:3
        k = 4*(i-1);
        u(:,i,:) = pot(:,i,:) + grad(:,k+4,:) ...
            - x1.*grad(:,k+1,:) - x2.*grad(:,k+2,:) - x3.*grad(:,k+3,:);
    end
end
if any(~isfinite(u(:)))
    error('stfmm3d_laplace:BackendFailure','Nonfinite Stokes velocity.');
end
u = squeeze(reshape(u,nd,3,n));
end

function self_test()
% Check the normalization against an axial Stokeslet at distance 2.
srcinfo = struct('sources',[0;0;0],'stoklet',[0;0;1]);
[U,timing] = stfmm3d_laplace(1e-12,srcinfo,0,[0;0;2],1);
assert(norm(U.pottarg(:)-[0;0;1/(8*pi)]) < 1e-13);
assert(timing.total_seconds >= timing.backend_call_seconds);
end
