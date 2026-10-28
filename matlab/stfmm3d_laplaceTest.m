function stfmm3d_laplaceTest()
% STFMM3D_LAPLACETEST    Regression checks for Laplace-based Stokeslet velocities.
%
% stfmm3d_laplaceTest()
%
%  Run with the MATLAB interface directory on the path and a working FMM3D
%  MEX binary. Both STFMM3D and LFMM3D are needed. 
%
%  The reference for numerical comparisons is STFMM3D, requested at the
%  same tolerance and with only Stokeslet velocity outputs enabled. These
%  are comparisons between two FMM formulations.
%
%  Coverage:
%    * One and two force densities; source-only, target-only and combined
%      evaluation, including omitted optional arguments and coincident points.
%    * Output field names, unrequested fields and singleton dimensions,
%      matching the STFMM3D interface and its squeeze convention.
%    * Identical velocities with and without the optional timing output;
%      finite, nonnegative stage times contained within the total time.
%      These checks do not measure speed or assert a performance threshold.
%    * Empty point sets, absent forces, disabled evaluations and explicit
%      errors for unsupported kernels/outputs or malformed force arrays.
%    * Translation invariance for a unit-sized cloud far from the origin,
%      with direct Stokeslet sums as the target velocity reference.
%
%  check_velocity requires matching array sizes and uses the mixed
%  absolute/relative bound ||actual-expected||_F <=
%  100*tol*max(1,||expected||_F). The factor 100 allows differences between
%  the two FMM approximations; the unit floor handles near-zero fields.
%
% See also: STFMM3D_LAPLACE, STFMM3D, STFMM3D_LAPLACEPERFTEST

rng(1);
tol = 1e-10;
ns = 513;
targ = rand(3,79)-0.5;

% Exercise the usual (3,ns) input and the (nd,3,ns) density layout.
for nd = [1,2]
    srcinfo = struct('sources',rand(3,ns)-0.5,'nd',nd, ...
        'stoklet',randn(nd,3,ns));
    if nd == 1, srcinfo.stoklet = reshape(srcinfo.stoklet,3,ns); end
    reference = stfmm3d(tol,srcinfo,1,targ,1);
    % Each column is [source flag; target flag]: sources, targets, then both.
    for flags = [1,0; 0,1; 1,1]'
        U = stfmm3d_laplace(tol,srcinfo,flags(1),targ,flags(2));
        [timed,timing] = stfmm3d_laplace(tol,srcinfo,flags(1),targ,flags(2));
        % Asking for timing must not change any returned velocity or field.
        assert(isequal(U,timed));
        seconds = cell2mat(struct2cell(timing));
        assert(all(isfinite(seconds)) && all(seconds >= 0));
        % Packing, backend and assembly intervals must fit within total time.
        assert(timing.total_seconds >= sum(seconds(1:3)));
        % Pressure/gradient fields exist for compatibility but remain empty.
        assert(isequal(fieldnames(U),fieldnames(reference)));
        assert(isempty(U.pre) && isempty(U.grad));
        assert(isempty(U.pretarg) && isempty(U.gradtarg));
        if flags(1)
            check_velocity(U.pot,reference.pot,tol);
        else
            assert(isempty(U.pot));
        end
        if flags(2)
            check_velocity(U.pottarg,reference.pottarg,tol);
        else
            assert(isempty(U.pottarg));
        end
    end
    % Three arguments request source evaluation without supplying targets.
    source_only = stfmm3d_laplace(tol,srcinfo,1);
    check_velocity(source_only.pot,reference.pot,tol);
    % Supplying targets without a target flag must still default that flag to 0.
    default_flag = stfmm3d_laplace(tol,srcinfo,1,targ);
    assert(isequal(default_flag,source_only));
    % Evaluating at the source coordinates through the target path must omit
    % self-interactions and agree with the source output from STFMM3D.
    coincident = stfmm3d_laplace(tol,srcinfo,0,srcinfo.sources,1);
    check_velocity(coincident.pottarg,reference.pot,tol);

    % Singleton output dimensions match STFMM3D for each density count.
    singleton = stfmm3d_laplace(tol,srcinfo,0,targ(:,1),1);
    direct = stfmm3d(tol,srcinfo,0,targ(:,1),1);
    check_velocity(singleton.pottarg,direct.pottarg,tol);
    if nd == 1
        % A single density may also be supplied with an explicit leading 1.
        srcinfo.stoklet = reshape(srcinfo.stoklet,1,3,ns);
        explicit_density = stfmm3d_laplace(tol,srcinfo,1);
        check_velocity(explicit_density.pot,reference.pot,tol);
    end
    % Restrict both locations and strengths to one source, then check the
    % zero self-velocity and the velocity at one distinct target.
    srcinfo.sources = srcinfo.sources(:,1);
    srcinfo.stoklet = reshape(srcinfo.stoklet,nd,3,ns);
    srcinfo.stoklet = srcinfo.stoklet(:,:,1);
    if nd == 1, srcinfo.stoklet = reshape(srcinfo.stoklet,3,1); end
    single_source = stfmm3d_laplace(tol,srcinfo,1,targ(:,1),1);
    direct = stfmm3d(tol,srcinfo,1,targ(:,1),1);
    check_velocity(single_source.pot,direct.pot,tol);
    check_velocity(single_source.pottarg,direct.pottarg,tol);
end

% The fourth charge uses y dot f. Without a common coordinate shift, moving
% a small cloud far from the origin introduces large gradient terms that
% must cancel in the reconstructed velocity. Check both density layouts.
% Coordinates are multiples of 2^-12 and the translation is a power of two:
% the translated coordinates retain exactly the same source-target offsets.
% This isolates reconstruction errors from loss of precision in input data.
shift_sources = round(4096*(rand(3,ns)-0.5))/4096;
shift_targets = round(4096*(rand(3,17)-0.5))/4096;
offset = 2^30*[1;-1;1];
assert(isequal((shift_sources+offset)-offset,shift_sources));
assert(isequal((shift_targets+offset)-offset,shift_targets));
for nd = [1,2]
    forces = randn(nd,3,ns);
    srcinfo = struct('sources',shift_sources,'nd',nd,'stoklet',forces);
    direct = direct_velocity(shift_sources,forces,shift_targets);
    centered = stfmm3d_laplace(tol,srcinfo,1,shift_targets,1);
    srcinfo.sources = shift_sources+offset;
    translated = stfmm3d_laplace(tol,srcinfo,1,shift_targets+offset,1);
    check_velocity(translated.pot,centered.pot,tol);
    check_velocity(translated.pottarg,direct,tol);
end

% With no sources, requested target velocities are zero and the recorded
% backend duration is zero because no LFMM3D call is needed.
srcinfo = struct('sources',zeros(3,0),'stoklet',zeros(3,0));
[U,timing] = stfmm3d_laplace(tol,srcinfo,1,targ,1);
assert(isempty(U.pot) && isequal(U.pottarg,zeros(size(targ))));
assert(timing.backend_call_seconds == 0);
% An omitted stoklet field also represents a zero velocity field.
srcinfo = struct('sources',rand(3,5));
U = stfmm3d_laplace(tol,srcinfo,1,targ,1);
assert(isequal(U.pot,zeros(3,5)) && isequal(U.pottarg,zeros(size(targ))));
srcinfo.stoklet = randn(3,5);
% Empty targets produce an empty target velocity array; they must not
% suppress a source evaluation requested in the same call.
U = stfmm3d_laplace(tol,srcinfo,0,zeros(3,0),1);
assert(isempty(U.pot) && isequal(U.pottarg,zeros(3,0)));
U = stfmm3d_laplace(tol,srcinfo,1,zeros(3,0),1);
direct = stfmm3d(tol,srcinfo,1);
check_velocity(U.pot,direct.pot,tol);
% Disabling both evaluations leaves every output field empty.
U = stfmm3d_laplace(tol,srcinfo,0);
assert(all(structfun(@isempty,U)));

% Unsupported physics must not silently produce partial results.
% Flags 2 and 3 request pressure/velocity gradients in STFMM3D; neither is
% supported here, on either the source or the target evaluation path.
for flag = [2,3]
    check_error(@() stfmm3d_laplace(tol,srcinfo,flag), ...
        'stfmm3d_laplace:UnsupportedOutput');
    check_error(@() stfmm3d_laplace(tol,srcinfo,0,targ,flag), ...
        'stfmm3d_laplace:UnsupportedOutput');
end
% Reject each unsupported kernel field, including an orientation supplied
% without its matching strength, rather than ignoring it silently.
fields = {'strslet','strsvec','rotlet','rotvec','doublet','doubvec'};
for k = 1:numel(fields)
    unsupported = srcinfo;
    unsupported.(fields{k}) = zeros(3,5);
    check_error(@() stfmm3d_laplace(tol,unsupported,1), ...
        'stfmm3d_laplace:UnsupportedKernel');
end
% There are five sources, so a force array with four columns is invalid.
srcinfo.stoklet = zeros(3,4);
check_error(@() stfmm3d_laplace(tol,srcinfo,1), ...
    'stfmm3d_laplace:InvalidStrengths');
% The wrapper's no-input self-test verifies the 1/(8*pi) normalization
% against an axial unit Stokeslet evaluated at distance 2.
stfmm3d_laplace();
fprintf('STFMM3D_LAPLACE interface, velocity and timing checks passed.\n');
end

function check_velocity(actual,expected,tol)
% Require the STFMM3D output shape as well as numerical agreement. Flattening
% gives a Frobenius norm over all components, points and force densities.
assert(isequal(size(actual),size(expected)));
assert(norm(actual(:)-expected(:)) <= 100*tol*max(1,norm(expected(:))));
end

function u = direct_velocity(sources,forces,targ)
% Independent Stokeslet formula using displacements, with no coordinate-
% weighted Laplace potentials. Reference targets are distinct from sources.
nd = size(forces,1);
ns = size(sources,2);
nt = size(targ,2);
u = zeros(nd,3,nt);
for k = 1:nt
    R = reshape(targ(:,k)-sources,1,3,ns);
    r = sqrt(sum(R.^2,2));
    assert(all(r(:) > 0));
    u(:,:,k) = sum(forces./r + R.*(sum(R.*forces,2)./r.^3),3)/(8*pi);
end
u = squeeze(u);
end

function check_error(call,identifier)
% Require the specific public error identifier. A different exception must
% fail the test, and returning normally must fail it as well.
try
    call();
catch err
    assert(strcmp(err.identifier,identifier));
    return
end
error('Expected error %s was not raised.',identifier);
end
