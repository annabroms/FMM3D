% Performance and accuracy test: compare STFMM3D and LFMM3D velocities.
% Edit the parameters below, then run this script. No files are written.
disp('Stokeslet velocity: STFMM3D versus four-potential LFMM3D')

source_counts = [512 4096 32768 131072];
target_counts = source_counts;   % collocation points; may be changed independently
tol = 1e-9;
nrepeat = 7;
ncheck = 16;               % direct checks sum all sources at these targets

rng(0)
assert(numel(source_counts) == numel(target_counts));
ncase = numel(source_counts);
summary = zeros(ncase,8);
raw_call_seconds = zeros(nrepeat,2,ncase);
raw_total_seconds = zeros(nrepeat,2,ncase);
for problem = 1:ncase
    ns = source_counts(problem);
    nt = target_counts(problem);
    sources = rand(3,ns) - 0.5;
    forces = randn(3,ns);
    targ = rand(3,nt) - 0.5;    % independent targets in the same unit cube

    % Warm up both interfaces before timing.
    srcinfo = struct('sources',sources,'stoklet',forces);
    stfmm3d(tol,srcinfo,0,targ,1);
    stfmm3d_laplace(sources,forces,targ,tol);

    % Each row is a paired repetition; columns are STFMM3D and Laplace.
    % Alternate execution order to reduce systematic timing bias.
    call_seconds = zeros(nrepeat,2);
    total_seconds = zeros(nrepeat,2);
    for repeat = 1:nrepeat
        order = [1 2];
        if mod(repeat,2) == 0, order = [2 1]; end
        for backend = order
            if backend == 1
                total_timer = tic;
                srcinfo = struct('sources',sources,'stoklet',forces);
                t = tic;
                U = stfmm3d(tol,srcinfo,0,targ,1);  % target velocity only
                call_seconds(repeat,1) = toc(t);
                u_stokes = reshape(U.pottarg,3,nt);
                assert(all(isfinite(u_stokes(:))));
                total_seconds(repeat,1) = toc(total_timer);
            else
                [u_laplace,timing] = stfmm3d_laplace(sources,forces,targ,tol);
                call_seconds(repeat,2) = timing.backend_call_seconds;
                total_seconds(repeat,2) = timing.total_seconds;
            end
        end
    end

    % Overhead is measured for each invocation before taking medians.
    overhead_seconds = total_seconds - call_seconds;
    t_call = median(call_seconds,1);
    t_total = median(total_seconds,1);
    t_overhead = median(overhead_seconds,1);
    fprintf('\n%d sources to %d targets, tol=%.3g; median of %d repeats\n',ns,nt,tol,nrepeat)
    fprintf('                 FMM call (s)    overhead (s)    total (s)\n')
    fprintf('STFMM3D          %12.4g    %12.4g    %9.4g\n',t_call(1),t_overhead(1),t_total(1))
    fprintf('Laplace          %12.4g    %12.4g    %9.4g\n',t_call(2),t_overhead(2),t_total(2))
    fprintf('STFMM3D/Laplace: call %.3gx, total %.3gx\n',t_call(1)/t_call(2),t_total(1)/t_total(2))

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
    % Constituent: max target error / sum of individual contribution magnitudes.
    fprintf('\nDirect checks at %d targets (all %d sources):\n',numel(selected),ns)
    fprintf('                 max absolute    relative        constituent\n')
    velocities = {u_stokes(:,selected), u_laplace(:,selected)};
    names = {'STFMM3D','Laplace'};
    for backend = 1:2
        difference = velocities{backend} - u_direct;
        absolute = max(vecnorm(difference));
        relative = NaN;           % undefined for a negligible reference
        if norm(u_direct,'fro') > 100*eps*norm(scale)
            relative = norm(difference,'fro')/norm(u_direct,'fro');
        end
        constituent = max(vecnorm(difference)./max(scale,realmin));
        fprintf('%-16s %12.3g    %12.3g    %12.3g\n',names{backend},absolute,relative,constituent)
        assert(constituent <= 100*tol+1e-12);
        assert(isnan(relative) || relative <= 100*tol+1e-11);
    end
    fprintf('Max STFMM3D/Laplace velocity difference (all targets): %.3g\n', ...
        max(vecnorm(u_stokes-u_laplace)))

    % Keep the paired samples and one summary row for each geometry.
    raw_call_seconds(:,:,problem) = call_seconds;
    raw_total_seconds(:,:,problem) = total_seconds;
    summary(problem,:) = [ns nt t_call t_total t_call(1)/t_call(2) t_total(1)/t_total(2)];
end

% Ratios above one mean the Laplace evaluator was faster.
speedup_table = array2table(summary, 'VariableNames', ...
    {'sources','targets','stokes_call_s','laplace_call_s', ...
    'stokes_total_s','laplace_total_s','call_speedup','total_speedup'});
fprintf('\nSpeedup across all four geometries (median timings):\n')
disp(speedup_table)
