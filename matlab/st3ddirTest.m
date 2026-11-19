% Direct Stokeslet velocities: st3ddir_opt versus st3ddir.
% Compares three ways of computing the velocity induced by Stokeslets:
%
%   full       st3ddir(srcinfo,targ,3)       velocity, pressure and gradient
%   velocity   st3ddir(srcinfo,targ,1)       velocity only
%   opt        st3ddir_opt(srcinfo,targ,1)   optimized velocity only
%
% Sections 1 and 2 check correctness on small problems; section 3 times
% the three methods. Build the MEX files with 'make matlab' in the FMM3D
% root first. OpenMP reads OMP_NUM_THREADS when MATLAB starts, so set it
% before starting MATLAB, e.g. from the FMM3D root:
%
%   OMP_NUM_THREADS=8 matlab -batch "addpath('matlab'); st3ddirTest"
%
% Results are saved as st3ddirTest-threads-<n>.mat, .json and .csv in
% matlab/results/st3ddirTest (ignored by git).

clear; close all;
rng(927)

tol = 1e-12;
source_counts = [32768 131072];
source_counts = [4096 32768];
target_counts = source_counts;
nrepeat = 3;
output_folder = fullfile(fileparts(mfilename('fullpath')),'results','st3ddirTest');

%% 1. st3ddir_opt agrees with st3ddir
% Several numbers of densities nd and all three output flags:
%   1 = velocity,  2 = velocity and pressure,  3 = all three outputs.
% Only flag 1 is optimized; for flags 2 and 3, st3ddir_opt calls the
% original routine. Outputs must also have the same shapes.
clear srcinfo
ns = 13;
nt = 9;
for nd = [1 2 3 4 8 16]
    srcinfo.sources = rand(3,ns);
    srcinfo.nd = nd;
    srcinfo.stoklet = squeeze(randn(nd,3,ns));

    % Two source-target pairs must be skipped: target 1 sits on source 1,
    % and target 2 lies 1e-16 from source 2, inside the exclusion threshold.
    targ = rand(3,nt);
    targ(:,1) = srcinfo.sources(:,1);
    targ(:,2) = srcinfo.sources(:,2) + [1e-16;0;0];

    for flag = 1:3
        U1 = st3ddir(srcinfo,targ,flag);
        % st3ddir_opt also returns the number of OpenMP threads it used.
        [U2,threads] = st3ddir_opt(srcinfo,targ,flag);
        for output = {'pottarg','pretarg','gradtarg'}
            name = output{1};
            assert(isequal(size(U2.(name)),size(U1.(name))) && ...
                relative_difference(U2.(name),U1.(name)) < tol, ...
                'st3ddir_opt %s differs (nd=%d, flag=%d)',name,nd,flag);
        end
    end
end
fprintf('1. st3ddir_opt matches st3ddir\n');

%% 2. st3ddir output flags with all four kernel types
% Stokeslets, stresslets, rotlets and doublets together. Asking for fewer
% outputs must not change the ones that are computed, and the outputs that
% were not requested must be empty.
clear srcinfo
ns = 13;
nt = 9;
srcinfo.sources = rand(3,ns);
srcinfo.stoklet = rand(3,ns);
srcinfo.strslet = rand(3,ns);
srcinfo.strsvec = rand(3,ns);
srcinfo.rotlet = rand(3,ns);
srcinfo.rotvec = rand(3,ns);
srcinfo.doublet = rand(3,ns);
srcinfo.doubvec = rand(3,ns);
targ = rand(3,nt);

U3 = st3ddir(srcinfo,targ,3);   % velocity, pressure, gradient
U2 = st3ddir(srcinfo,targ,2);   % velocity, pressure
U1 = st3ddir(srcinfo,targ,1);   % velocity

assert(relative_difference(U1.pottarg,U3.pottarg) < tol, ...
    'Velocity from flag 1 differs from flag 3');
assert(relative_difference(U2.pretarg,U3.pretarg) < tol, ...
    'Pressure from flag 2 differs from flag 3');
assert(isempty(U1.pretarg) && isempty(U1.gradtarg) && isempty(U2.gradtarg), ...
    'Outputs that were not requested must be empty');

% st3ddir_opt handles Stokeslets only and must refuse the other kernels.
try
    st3ddir_opt(srcinfo,targ,1);
    rejected = false;
catch exception
    rejected = strcmp(exception.identifier,'st3ddir_opt:kernel');
end
assert(rejected,'st3ddir_opt must reject stresslets, rotlets and doublets');
fprintf('2. st3ddir output flags agree for all kernels\n');

%% 3. Timing
% Every source count is timed with one density; the first is also timed
% with four. Sources and targets are random points in the unit cube.
ncase = numel(source_counts) + 1;
case_ns = [source_counts source_counts(1)];
case_nt = [target_counts target_counts(1)];
case_nd = [ones(1,numel(source_counts)) 4];

fprintf('\n3. OpenMP threads: %d. Median times in seconds over %d runs; ', ...
    threads,nrepeat);
fprintf('ratios show speedups.\n');
fprintf('%8s %8s %3s %12s %12s %12s %9s %9s %9s %10s\n','ns','nt','nd', ...
    'full','velocity','opt','full/vel','vel/opt','full/opt','rel. diff');

rows = cell(ncase,1);
for k = 1:ncase
    clear srcinfo
    ns = case_ns(k);
    nt = case_nt(k);
    nd = case_nd(k);
    srcinfo.sources = rand(3,ns);
    srcinfo.nd = nd;
    srcinfo.stoklet = squeeze(randn(nd,3,ns));
    targ = rand(3,nt);

    methods = {@() st3ddir(srcinfo,targ,3), ...      % full
               @() st3ddir(srcinfo,targ,1), ...      % velocity
               @() st3ddir_opt(srcinfo,targ,1)};     % opt

    % Run each method once untimed, to load the MEX files and warm the caches.
    U = cellfun(@(method) method(),methods,'UniformOutput',false);

    % Then time each method. The order rotates every repeat (1-2-3, then
    % 2-3-1, then 3-1-2, ...) so that no method always runs first.
    samples = zeros(nrepeat,3);   % seconds, one column per method
    for irepeat = 1:nrepeat
        for m = circshift(1:3,1-irepeat)
            started = tic;
            U{m} = methods{m}();
            samples(irepeat,m) = toc(started);
        end
    end
    seconds = median(samples,1);

    % All three methods must return the same velocity.
    err = max(relative_difference(U{2}.pottarg,U{1}.pottarg), ...
              relative_difference(U{3}.pottarg,U{1}.pottarg));
    assert(err < tol,'Timed velocities differ by %.2e',err);

    rows{k} = struct('ns',ns,'nt',nt,'nd',nd,'threads',threads, ...
        'samples',samples,'seconds',seconds, ...
        'speedup_full_over_velocity',seconds(1)/seconds(2), ...
        'speedup_velocity_over_opt',seconds(2)/seconds(3), ...
        'speedup_full_over_opt',seconds(1)/seconds(3), ...
        'relative_velocity_error',err);
    fprintf('%8d %8d %3d %12.6g %12.6g %12.6g %8.2fx %8.2fx %8.2fx %10.2e\n', ...
        ns,nt,nd,seconds,seconds(1)/seconds(2),seconds(2)/seconds(3), ...
        seconds(1)/seconds(3),err);
end
rows = vertcat(rows{:});


function err = relative_difference(observed, expected)
% Size of the difference between two arrays, relative to the size of
% EXPECTED. Two empty arrays have difference 0.
err = norm(observed(:)-expected(:)) / max(norm(expected(:)),realmin);
end
