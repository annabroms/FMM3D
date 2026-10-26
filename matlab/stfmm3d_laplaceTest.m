function report = stfmm3d_laplaceTest()
%STFMM3D_LAPLACETEST Deterministic correctness only; no timing comparison.
% Compare independent direct Stokeslets, upstream stfmm3d, and the new adapter.
% Missing compiled dependencies raise an error, never a passing/skipped result.
assert(exist('fmm3d','file')==3 && exist('lfmm3d','file')==2 && exist('stfmm3d','file')==2, ...
    'stfmm3d_laplace:MissingDependency','Compiled FMM3D MATLAB interfaces are required.');
state=rng; cleanup=onCleanup(@() rng(state));
rng(9217); tolerance=1e-10; max_relative=0; max_constituent=0; cases=0;
stfmm3d_laplace();
for fixture=1:7
    switch fixture
        case {1,2,3}
            sources=[0;0;0]; forces=zeros(3,1); forces(fixture)=1;
            targets=[0.2 0.4 -0.7;0.3 -0.8 0.1;0.9 0.5 -0.2];
        case 4
            sources=rand(3,17)-0.5; forces=randn(3,17); targets=[0.8;-0.1;0.3];
        case 5
            sources=rand(3,19)-0.5; forces=randn(3,19); targets=rand(3,31)-0.5;
        case 6
            % More sources than a leaf: include multipole interactions.
            sources=rand(3,1200)-0.5; forces=randn(3,1200); targets=rand(3,37)-0.5;
        case 7
            sources=[0;0;0]; forces=[1;-2;3]; targets=[0.4;0.2;-0.3];
    end
    [reference,scale]=direct_reference(sources,forces,targets);
    [velocity,timing]=stfmm3d_laplace(sources,forces,targets,tolerance);
    upstream=stfmm3d(tolerance,struct('sources',sources,'stoklet',forces),0,targets,1);
    assert(isequal(size(velocity),[3 size(targets,2)]));
    for value={velocity,reshape(upstream.pottarg,3,[])}
        relative=norm(value{1}-reference,'fro')/norm(reference,'fro');
        constituent=max(vecnorm(value{1}-reference)./scale);
        assert(relative<=100*tolerance+1e-11);
        assert(constituent<=100*tolerance+1e-12);
        max_relative=max(max_relative,relative);
        max_constituent=max(max_constituent,constituent);
    end
    assert(timing.total_seconds>=timing.backend_call_seconds);
    assert(timing.total_seconds>=timing.packing_seconds+timing.backend_call_seconds+timing.assembly_seconds);
    % Shared translation and dimensional scaling preserve the physical kernel.
    shift=[10;-20;30];
    moved=stfmm3d_laplace(sources+shift,forces,targets+shift,tolerance);
    assert(norm(moved-reference,'fro')<=1e-8*norm(reference,'fro'));
    for factor=[1e-9 1e9]
        scaled=stfmm3d_laplace(factor*sources,forces,factor*targets,tolerance);
        assert(norm(factor*scaled-reference,'fro')<=1e-8*norm(reference,'fro'));
    end
    cases=cases+1;
end
empty=zeros(3,0);
[u,t]=stfmm3d_laplace(empty,empty,zeros(3,2));
assert(isequal(u,zeros(3,2)) && t.backend_call_seconds==0);
u=stfmm3d_laplace(zeros(3,1),ones(3,1),empty);
assert(isequal(size(u),[3 0]));
assert_rejected(@() stfmm3d_laplace(zeros(2,1),ones(3,1),zeros(3,1)));
assert_rejected(@() stfmm3d_laplace(zeros(3,2),ones(3,1),zeros(3,1)));
assert_rejected(@() stfmm3d_laplace([NaN;0;0],ones(3,1),zeros(3,1)));
assert_rejected(@() stfmm3d_laplace(zeros(3,1),ones(3,1),zeros(3,1),0));
report=struct('passed',true,'fixtures',cases,'tolerance',tolerance, ...
    'max_global_relative_error',max_relative,'max_constituent_error',max_constituent);
disp(report);
end

function [velocity,scale] = direct_reference(sources,forces,targets)
% Independent explicit Stokeslet tensor gives a small-fixture reference.
velocity=zeros(3,size(targets,2)); scale=zeros(1,size(targets,2));
for m=1:size(targets,2)
    for n=1:size(sources,2)
        r=targets(:,m)-sources(:,n); radius=norm(r);
        tensor=(eye(3)/radius+(r*r.')/radius^3)/(8*pi);
        contribution=tensor*forces(:,n);
        velocity(:,m)=velocity(:,m)+contribution;
        scale(m)=scale(m)+norm(contribution);
    end
end
end

function assert_rejected(evaluation)
% Ensure invalid public inputs fail rather than silently reshape.
rejected=false;
try
    evaluation();
catch
    rejected=true;
end
assert(rejected,'Invalid input was accepted.');
end
