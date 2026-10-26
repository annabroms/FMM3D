function comparison = stfmm3d_laplacePerfTest(input)
%STFMM3D_LAPLACEPERFTEST Opt-in free-space velocity timing AND accuracy test.
% No wall, reflected sources, harmonic correction, or automatic test execution.
% comparison = stfmm3d_laplacePerfTest(config)
% Default Ns=[512 4096], Nt/Ns=[.25 1 4], tol=1e-9, seven paired repeats.
% config.save=false disables MAT/CSV output in matlab/results.
% All config fields are listed in the defaults below. No thread changes.
if nargin==0, input=struct(); end
state=rng; cleanup=onCleanup(@() rng(state));
c=struct('sizes',[512 4096],'ratios',[0.25 1 4],'tol',1e-9, ...
    'repetitions',7,'seed',481516,'full_pair_limit',2e7, ...
    'subset_size',256,'block_size',[256 2048],'save',true);
names=fieldnames(input);
for k=1:numel(names)
    assert(isfield(c,names{k}),'Unknown comparison option: %s',names{k});
    c.(names{k})=input.(names{k});
end
validateattributes(c.sizes,{'double'},{'integer','positive','vector','finite'});
validateattributes(c.ratios,{'double'},{'positive','finite','vector'});
validateattributes(c.tol,{'double'},{'real','finite','scalar','>',0,'<',1});
validateattributes(c.repetitions,{'double'},{'integer','positive','scalar'});
validateattributes(c.seed,{'double'},{'integer','nonnegative','scalar','finite'});
validateattributes(c.full_pair_limit,{'double'},{'positive','finite','scalar'});
validateattributes(c.subset_size,{'double'},{'integer','positive','scalar'});
validateattributes(c.block_size,{'double'},{'integer','positive','numel',2});
validateattributes(c.save,{'logical'},{'scalar'});
c.sizes=reshape(c.sizes,1,[]); c.ratios=reshape(c.ratios,1,[]);
root=fileparts(mfilename('fullpath'));
comparison=struct('configuration',c,'timing_schema_version',2, ...
    'problem','free-space Stokeslet velocity, mu=1', ...
    'metadata',dependency_manifest());
backends={'upstream','laplace-velocity'};
rows={}; differences={}; all_samples={}; case_id=0;
for ns=c.sizes
    for ratio=c.ratios
        case_id=case_id+1; nt=max(1,round(ns*ratio));
        seed=c.seed+case_id; rng(seed);
        sources=rand(3,ns)-0.5; targets=rand(3,nt)-0.5; forces=randn(3,ns);
        % All input and coincidence checks precede timed evaluations.
        assert(all(isfinite([sources(:);targets(:);forces(:)])));
        assert(isempty(intersect(sources.',targets.','rows')), ...
            'stfmm3d_laplace:SingularInteraction','A generated source coincides with a target.');
        full=double(ns)*nt<=c.full_pair_limit;
        selected=1:nt;
        if ~full
            [~,order]=sortrows(targets.');
            selected=reshape(order(unique(round(linspace(1,nt,min(c.subset_size,nt))))),1,[]);
        end
        reference=direct_stokes(sources,forces,targets(:,selected), ...
            c.block_size);
        velocities=cell(1,2); errors=cell(1,2);
        % Warm both backends on the exact case, then check accuracy outside timers.
        for k=1:2
            [velocities{k},~]=run_backend(backends{k},sources,forces,targets,c.tol);
            errors{k}=error_metrics(velocities{k}(:,selected),reference.velocity,reference.scale);
            check_accuracy(errors{k},c.tol);
        end
        difference=error_metrics(velocities{2}(:,selected),velocities{1}(:,selected),reference.scale);
        repetitions=cell(2*c.repetitions,1); index=0;
        for repeat=1:c.repetitions
            execution_order=[1 2];
            if mod(repeat,2)==0, execution_order=[2 1]; end
            for position=1:2
                k=execution_order(position);
                [velocity,timing]=run_backend(backends{k},sources,forces,targets,c.tol);
                sample=struct('case_id',case_id,'repetition',repeat,'execution_position',position, ...
                    'ns',ns,'nt',nt,'stokes_backend',backends{k}, ...
                    'backend_call_seconds',timing.backend_call_seconds, ...
                    'packing_seconds',timing.packing_seconds,'assembly_seconds',timing.assembly_seconds, ...
                    'total_seconds',timing.total_seconds, ...
                    'overhead_seconds',timing.total_seconds-timing.backend_call_seconds);
                index=index+1; repetitions{index}=sample;
                % Validate each returned velocity after its timed invocation.
                error=error_metrics(velocity(:,selected),reference.velocity,reference.scale);
                check_accuracy(error,c.tol);
            end
        end
        samples=struct2table(vertcat(repetitions{:})); all_samples{end+1}=samples; %#ok<AGROW>
        summaries=cell(1,2);
        for k=1:2
            chosen=strcmp(samples.stokes_backend,backends{k});
            summaries{k}=summarize_timings(samples(chosen,:));
            row=struct('case_id',case_id,'ns',ns,'nt',nt,'seed',seed, ...
                'requested_tolerance',c.tol,'repetitions',c.repetitions, ...
                'stokes_backend',backends{k},'full_direct',full, ...
                'reference_target_count',numel(selected),'max_absolute_error',errors{k}.max_absolute, ...
                'global_relative_error',errors{k}.global_relative, ...
                'constituent_error',errors{k}.constituent_normalized);
            timing_names=fieldnames(summaries{k});
            for j=1:numel(timing_names), row.(timing_names{j})=summaries{k}.(timing_names{j}); end
            rows{end+1}=row; %#ok<AGROW>
        end
        a=summaries{1}; b=summaries{2};
        differences{end+1}=struct('case_id',case_id,'ns',ns,'nt',nt, ...
            'call_seconds_saved',a.backend_call_seconds-b.backend_call_seconds, ...
            'total_seconds_saved',a.total_seconds-b.total_seconds, ...
            'call_time_ratio',a.backend_call_seconds/b.backend_call_seconds, ...
            'total_time_ratio',a.total_seconds/b.total_seconds, ...
            'backend_max_absolute_difference',difference.max_absolute, ...
            'backend_global_relative_difference',difference.global_relative, ...
            'backend_constituent_difference',difference.constituent_normalized); %#ok<AGROW>
        fprintf('Free-space Ns=%d Nt=%d: upstream %.4gs, laplace-velocity %.4gs, total ratio %.3g\n', ...
            ns,nt,a.total_seconds,b.total_seconds,a.total_seconds/b.total_seconds);
    end
end
comparison.summary=struct2table(vertcat(rows{:}));
comparison.differences=struct2table(vertcat(differences{:}));
comparison.samples=vertcat(all_samples{:});
if c.save
    destination=fullfile(root,'results'); if ~isfolder(destination), mkdir(destination); end
    stem=tempname(destination); [~,name]=fileparts(stem);
    comparison.output_stem=fullfile(destination,['stokes_backends_' name]);
    save([comparison.output_stem '.mat'],'comparison');
    writetable(comparison.summary,[comparison.output_stem '_summary.csv']);
    writetable(comparison.differences,[comparison.output_stem '_differences.csv']);
    writetable(comparison.samples,[comparison.output_stem '_samples.csv']);
end
disp(comparison.summary); disp(comparison.differences);
end

function [velocity,timing]=run_backend(backend,sources,forces,targets,tolerance)
% Isolate one free-space evaluation so packing/reconstruction belong to its total.
if strcmp(backend,'laplace-velocity')
    [velocity,timing]=stfmm3d_laplace(sources,forces,targets,tolerance);
else
    total_timer=tic; timer=tic;
    src=struct('sources',sources,'stoklet',forces);
    timing.packing_seconds=toc(timer);
    timer=tic; value=stfmm3d(tolerance,src,0,targets,1);
    timing.backend_call_seconds=toc(timer);
    timer=tic; velocity=reshape(value.pottarg,3,size(targets,2));
    assert(all(isfinite(velocity(:))),'stfmm3d_laplace:FMMFailure','Nonfinite upstream output.');
    timing.assembly_seconds=toc(timer);
    timing.total_seconds=toc(total_timer);
end
end

function check_accuracy(error,tolerance)
% Accept accuracy using measured error; never require a timing speedup.
assert(error.constituent_normalized<=100*tolerance+1e-12, ...
    'stfmm3d_laplace:AccuracyFailure','Constituent-normalized Stokes velocity error is too large.');
if ~isnan(error.global_relative)
    assert(error.global_relative<=100*tolerance+1e-11, ...
        'stfmm3d_laplace:AccuracyFailure','Global relative Stokes velocity error is too large.');
end
end

function out = direct_stokes(sources,forces,targets,block_size)
% Independently sum blocked free-space Stokeslets and constituent magnitudes.
nt=size(targets,2); ns=size(sources,2);
out=struct('velocity',zeros(3,nt),'scale',zeros(1,nt));
for t0=1:block_size(1):nt
    ti=t0:min(t0+block_size(1)-1,nt);
    for s0=1:block_size(2):ns
        si=s0:min(s0+block_size(2)-1,ns);
        r=cell(1,3);
        for k=1:3, r{k}=targets(k,ti).'-sources(k,si); end
        invr=1./sqrt(r{1}.^2+r{2}.^2+r{3}.^2);
        f=forces(:,si);
        rf=r{1}.*f(1,:)+r{2}.*f(2,:)+r{3}.*f(3,:);
        magnitude=zeros(size(invr));
        for i=1:3
            v=(f(i,:).*invr+r{i}.*rf.*invr.^3)/(8*pi);
            out.velocity(i,ti)=out.velocity(i,ti)+sum(v,2).';
            magnitude=magnitude+v.^2;
        end
        out.scale(ti)=out.scale(ti)+sum(sqrt(magnitude),2).';
    end
end
end


function value = error_metrics(velocity,reference,scale)
% Normalize errors by the reference and by sourcewise contribution magnitudes.
pointwise=vecnorm(velocity-reference);
value.max_absolute=max([0 pointwise]);
value.global_relative=NaN;
if norm(reference,'fro')>100*eps*norm(scale)
    value.global_relative=norm(velocity-reference,'fro')/norm(reference,'fro');
end
normalized=zeros(size(scale)); nonzero=scale>0;
normalized(nonzero)=pointwise(nonzero)./scale(nonzero);
normalized(~nonzero & pointwise>0)=Inf;
value.constituent_normalized=max([0 normalized]);
end

function summary = summarize_timings(samples)
%SUMMARIZE_TIMINGS Validate same-invocation timings and median each raw column.
names=samples.Properties.VariableNames;
if ismember('backend_call_seconds',names)
    calls=samples.backend_call_seconds;
else
    calls=samples.stokes_seconds+samples.laplace_seconds;
    assert(all(abs(calls-samples.fmm_evaluation_seconds)<=64*eps.*max(1,calls)), ...
        'stfmm3d_laplace:InvalidTiming','FMM components do not match combined times.');
end
slack=64*eps.*max(1,samples.total_seconds);
assert(all(calls>=0 & samples.overhead_seconds>=-slack), ...
    'stfmm3d_laplace:InvalidTiming','Negative call time or overhead.');
assert(all(abs(samples.total_seconds-calls-samples.overhead_seconds)<=slack), ...
    'stfmm3d_laplace:InvalidTiming','Timing values must describe the same invocation.');
summary=struct();
for k=1:numel(names)
    if endsWith(names{k},'_seconds')
        summary.(names{k})=median(samples.(names{k}));
    end
end
end

function info = dependency_manifest()
%MANIFEST Record source and binary provenance separately; no source checkout edits.
assert(exist('fmm3d','file')==3 && exist('lfmm3d','file')==2 && exist('stfmm3d','file')==2, ...
    'stfmm3d_laplace:MissingDependency','FMM3D compiled MATLAB interfaces are required.');
info=struct('matlab',version,'computer',computer, ...
    'stokes_path',which('stfmm3d'),'laplace_path',which('lfmm3d'), ...
    'mex_path',which('fmm3d'),'omp_num_threads',getenv('OMP_NUM_THREADS'), ...
    'matlab_max_threads',maxNumCompThreads, ...
    'stokes_error_flag_exposed',false,'tree_reuse',false, ...
    'timing_schema_version',2,'velocity_default_backend','laplace-velocity', ...
    'available_stokes_backends',{{'upstream','laplace-velocity'}});
root=fileparts(fileparts(info.stokes_path));
% Quote filesystem paths before passing them to read-only shell commands.
quote = @(s) ['''' strrep(s,'''','''"''"''') ''''];
[status,revision]=system(['git -C ' quote(root) ' rev-parse HEAD']);
info.source_revision='unknown';
if status==0, info.source_revision=strtrim(revision); end
[status,checksum]=system(['sha256sum ' quote(info.mex_path)]);
info.mex_sha256='unknown';
if status==0, info.mex_sha256=extractBefore(strtrim(checksum),' '); end
[~,info.cpu]=system('lscpu');
[~,info.memory]=system('free -b');
[~,info.os]=system('uname -a');
end
