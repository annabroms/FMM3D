% Compare complete STFMM3D and Laplace-based Stokeslet velocity calls.
% Edit the parameters below, then run this script. No files are written.
disp('Stokeslet velocity: STFMM3D versus four-potential LFMM3D')
disp('Total = complete single-output interface call, including all wrapper work.')
disp('Both methods use the same srcinfo, targets, tolerance and velocity outputs.')

source_counts = [4096 32768 131072];
target_counts = source_counts;   % collocation points; may be changed independently
tol = 1e-8;
nrepeat = 7;
ncheck = 16;               % direct checks sum all sources at these targets
show_laplace_breakdown = false; % separate diagnostic calls, outside comparison

rng(0)
assert(numel(source_counts) == numel(target_counts));
ncase = numel(source_counts);
summary = zeros(ncase,5);
raw_total_seconds = zeros(nrepeat,2,ncase);
raw_laplace_breakdown_seconds = zeros(nrepeat,4,ncase*show_laplace_breakdown);
for problem = 1:ncase
    ns = source_counts(problem);
    nt = target_counts(problem);
    sources = rand(3,ns) - 0.5;
    forces = randn(3,ns);
    targ = rand(3,nt) - 0.5;    % independent targets in the same unit cube

    % Warm up both interfaces before timing.
    srcinfo = struct('sources',sources,'stoklet',forces);
    stfmm3d(tol,srcinfo,0,targ,1);
    stfmm3d_laplace(tol,srcinfo,0,targ,1);

    % Each row is a paired repetition; columns are STFMM3D and Laplace.
    % Alternate execution order to reduce systematic timing bias.
    total_seconds = zeros(nrepeat,2);
    for repeat = 1:nrepeat
        order = [1 2];
        if mod(repeat,2) == 0, order = [2 1]; end
        for backend = order
            % Shared inputs are already built. Time only the complete call;
            % both return the same velocity structure without internal timers.
            total_timer = tic;
            if backend == 1
                U = stfmm3d(tol,srcinfo,0,targ,1);  % target velocity only
            else
                U = stfmm3d_laplace(tol,srcinfo,0,targ,1);
            end
            total_seconds(repeat,backend) = toc(total_timer);
            % Common output extraction and verification are outside timing.
            assert(all(isfinite(U.pottarg(:))));
            if backend == 1
                u_stokes = reshape(U.pottarg,3,nt);
            else
                u_laplace = reshape(U.pottarg,3,nt);
            end
        end
    end

    t_total = median(total_seconds,1);
    fprintf('\n%d sources to %d targets, tol=%.3g; median of %d repeats\n',ns,nt,tol,nrepeat)
    fprintf('                 total (s)\n')
    fprintf('STFMM3D          %9.4g\n',t_total(1))
    fprintf('Laplace          %9.4g\n',t_total(2))
    fprintf('STFMM3D/Laplace total speedup: %.3gx\n',t_total(1)/t_total(2))

    if show_laplace_breakdown
        % Internal timers are enabled only in these extra diagnostic calls.
        for repeat = 1:nrepeat %#ok<UNRCH> Enabled by the diagnostic flag above.
            [~,timing] = stfmm3d_laplace(tol,srcinfo,0,targ,1);
            other = timing.total_seconds - timing.packing_seconds ...
                - timing.backend_call_seconds - timing.assembly_seconds;
            raw_laplace_breakdown_seconds(repeat,:,problem) = ...
                [timing.packing_seconds, timing.backend_call_seconds, ...
                timing.assembly_seconds, other];
        end
        stages = median(raw_laplace_breakdown_seconds(:,:,problem),1);
        fprintf('Separate Laplace diagnostic (s): packing %.4g, LFMM3D %.4g, assembly %.4g, other %.4g\n',stages)
    end

    % Check both velocities against the direct Stokeslet formula.
    % This work is outside the timings above. Only ncheck targets are checked.
    selected = unique(round(linspace(1,nt,min(nt,ncheck))));
    u_direct = zeros(3,numel(selected));
    scale = zeros(1,numel(selected));
    for m = 1:numel(selected)
        for j = 1:ns
            R = targ(:,selected(m)) - sources(:,j);  % target minus source
            r = norm(R);
            f = forces(:,j);
            contribution = (f/r + R*dot(f,R)/r^3)/(8*pi);
            u_direct(:,m) = u_direct(:,m) + contribution;
            scale(m) = scale(m) + norm(contribution);
        end
    end

    % Absolute: max target error norm. Relative: Frobenius error / reference norm.
    % Contribution-normalized error: at each checked target, divide the
    % velocity error norm by the sum of the norms of individual source
    % contributions, then take the maximum over the checked targets:
    % max_x ||u_FMM(x)-u_direct(x)||_2 / sum_j ||G(x,y_j)*f_j||_2.
    % The denominator measures contributions before cancellation, so this
    % error does not grow just because they cancel to a small net velocity.
    % realmin protects the division when all contributions are zero.
    fprintf('\nDirect checks at %d targets (all %d sources):\n',numel(selected),ns)
    fprintf('                 max absolute    relative        contribution-normalized error\n')
    velocities = {u_stokes(:,selected), u_laplace(:,selected)};
    names = {'STFMM3D','Laplace'};
    for backend = 1:2
        difference = velocities{backend} - u_direct;
        absolute = max(vecnorm(difference));
        relative = NaN;           % undefined for a negligible reference
        if norm(u_direct,'fro') > 100*eps*norm(scale)
            relative = norm(difference,'fro')/norm(u_direct,'fro');
        end
        contribution_normalized_error = max(vecnorm(difference)./max(scale,realmin));
        fprintf('%-16s %12.3g    %12.3g    %12.3g\n', ...
            names{backend},absolute,relative,contribution_normalized_error)
        assert(contribution_normalized_error <= 100*tol+1e-12);
        assert(isnan(relative) || relative <= 100*tol+1e-11);
    end
    fprintf('Max STFMM3D/Laplace velocity difference (all targets): %.3g\n', ...
        max(vecnorm(u_stokes-u_laplace)))

    % Keep the paired samples and one summary row for each geometry.
    raw_total_seconds(:,:,problem) = total_seconds;
    summary(problem,:) = [ns nt t_total t_total(1)/t_total(2)];
end

% Ratios above one mean the Laplace evaluator was faster.
speedup_table = array2table(summary, 'VariableNames', ...
    {'sources','targets','stokes_total_s','laplace_total_s','total_speedup'});
fprintf('\nSpeedup across all geometries (median complete-call timings):\n')
disp(speedup_table)

%% Translation accuracy: with and without coordinate recentering
% This separate accuracy experiment is outside all performance timings.
% Move the same unit-sized cloud by (a,-a,a). Direct Stokeslet sums depend
% only on displacements, so their exact values should not change.
% Coordinates are multiples of 2^-12: for the offsets below, translations
% are exactly representable. Verify this so input rounding cannot explain
% an observed loss of accuracy in the coordinate-weighted reconstruction.
rng(73)
shift_ns = 4096;
shift_nt = 32;
shift_tol = 1e-9;
shift_offsets = [0 1e2 1e4 1e6 1e8 1e10 1e12];
shift_sources = round(4096*(rand(3,shift_ns)-0.5))/4096;
shift_targets = round(4096*(rand(3,shift_nt)-0.5))/4096;
shift_forces = randn(3,shift_ns);

% Independent reference: sum the Stokeslet kernel at every test target,
% without using Laplace potentials or absolute-coordinate factors.
shift_direct = zeros(3,shift_nt);
for target = 1:shift_nt
    R = shift_targets(:,target)-shift_sources;
    r = sqrt(sum(R.^2,1));
    assert(all(r > 0),'Translation test targets must differ from sources.');
    shift_direct(:,target) = sum(shift_forces./r ...
        + R.*(sum(R.*shift_forces,1)./r.^3),2)/(8*pi);
end
shift_reference_norm = norm(shift_direct,'fro');
assert(shift_reference_norm > 0);
shift_relative_errors = zeros(numel(shift_offsets),4);

for offset_index = 1:numel(shift_offsets)
    offset_vector = shift_offsets(offset_index)*[1;-1;1];
    test_srcinfo = struct('sources',shift_sources+offset_vector, ...
        'stoklet',shift_forces);
    translated_targets = shift_targets+offset_vector;
    assert(isequal(test_srcinfo.sources-offset_vector,shift_sources) && ...
        isequal(translated_targets-offset_vector,shift_targets), ...
        'Chosen translation loses precision in the input coordinates.');

    % Native STFMM3D uses the original coordinates. For its shifted version,
    % recenter both input point sets before calling the same native routine.
    U_native = stfmm3d(shift_tol,test_srcinfo,0,translated_targets,1);
    lower = min(min(test_srcinfo.sources,[],2),min(translated_targets,[],2));
    upper = max(max(test_srcinfo.sources,[],2),max(translated_targets,[],2));
    test_center = lower/2+upper/2;
    centered_srcinfo = test_srcinfo;
    centered_srcinfo.sources = test_srcinfo.sources-test_center;
    U_native_shifted = stfmm3d(shift_tol,centered_srcinfo,0, ...
        translated_targets-test_center,1);

    % The local helper omits recentering deliberately. The production
    % Laplace wrapper applies its normal shift internally.
    u_laplace_unshifted = unshifted_laplace_velocity(shift_tol, ...
        test_srcinfo,translated_targets);
    U_laplace_shifted = stfmm3d_laplace(shift_tol,test_srcinfo,0,translated_targets,1);
    shift_velocities = {U_native.pottarg,U_native_shifted.pottarg, ...
        u_laplace_unshifted,U_laplace_shifted.pottarg};
    for method = 1:4
        velocity = reshape(shift_velocities{method},3,shift_nt);
        assert(all(isfinite(velocity(:))));
        shift_relative_errors(offset_index,method) = ...
            norm(velocity-shift_direct,'fro')/shift_reference_norm;
    end
end

fprintf('\nTranslation accuracy: %d sources, %d targets, tol=%.3g\n', ...
    shift_ns,shift_nt,shift_tol)
disp('Offset a means translation by (a,-a,a); all source-target separations are unchanged.')
disp('Errors are ||u-u_direct||_F / ||u_direct||_F over all test targets.')
disp('Shifted columns use the same bounding-box center for sources and targets.')
translation_error_table = array2table([shift_offsets(:),shift_relative_errors], ...
    'VariableNames',{'offset_a','stokes_unshifted_relerr','stokes_shifted_relerr', ...
    'laplace_unshifted_relerr','laplace_shifted_relerr'});
disp(translation_error_table)
% The unshifted errors intentionally illustrate cancellation and may exceed
% the requested tolerance. Only the two shifted paths must meet this bound.
assert(all(all(shift_relative_errors(:,[2,4]) <= 100*shift_tol+1e-12)));

function u = unshifted_laplace_velocity(tol,srcinfo,targ)
% Demonstration helper for one density and target velocities only. Use the
% same four Laplace charges and reconstruction as STFMM3D_LAPLACE, but keep
% the original coordinates to expose cancellation at large origin offsets.
lapinfo.sources = srcinfo.sources;
lapinfo.nd = 4;
lapinfo.charges = [srcinfo.stoklet/2; ...
    sum(srcinfo.sources.*srcinfo.stoklet,1)/2];
[L,ier] = lfmm3d(tol,lapinfo,0,targ,2);
if ier ~= 0
    error('stfmm3d_laplacePerfTest:BackendFailure','LFMM3D returned ier=%d.',ier);
end
nt = size(targ,2);
pot = reshape(L.pottarg,4,nt);
grad = reshape(L.gradtarg,12,nt);
u = zeros(3,nt);
for i = 1:3
    k = 4*(i-1);
    u(i,:) = pot(i,:) + grad(k+4,:) ...
        - targ(1,:).*grad(k+1,:) - targ(2,:).*grad(k+2,:) - targ(3,:).*grad(k+3,:);
end
end
