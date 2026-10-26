# Stokeslet target velocity through Laplace

These additions are independent of the Blake half-space wrapper and leave the
installed Fortran sources and MEX unchanged. Add this `matlab` folder to the path.

```matlab
[u, timing] = stfmm3d_laplace(sources, forces, targets, 1e-9);
report = stfmm3d_laplaceTest;       % small deterministic correctness checks
% Opt-in timing AND accuracy comparison (not part of correctness tests):
comparison = stfmm3d_laplacePerfTest;
% comparison = stfmm3d_laplacePerfTest(struct('sizes',[512 2048],'save',false));
```

Inputs are finite real double arrays: `sources` and `forces` are `3 x Ns`,
`targets` is `3 x Nt`. Output is always `3 x Nt`, including singleton/empty
cases. Tolerance defaults to `1e-9`. The kernel is
`(I/r + r*r'/r^3)/(8*pi)` at unit viscosity; divide velocity by viscosity when
needed. This array interface is deliberately limited to target Stokeslet
velocity. It does not accept the general `stfmm3d` source structure, stresslets,
source evaluation, pressure or gradients. Use `stfmm3d` for those fields.

Callers must supply nonsingular, adequately resolved targets. Like FMM3D,
interactions below its domain-relative distance cutoff are omitted; this
standalone helper does not include Blake's geometry proximity guard. The
requested tolerance is not a relative-error guarantee for a cancelling sum.

With `G=1/(4*pi*r)`, one `lfmm3d` call evaluates four charge potentials with
strengths `[forces/2; sum(sources.*forces,1)/2]`, requesting no source fields
and target potential/gradient only (`pg=0`, `pgt=2`, `nd=4`). Reconstruct
`u_i = psi_i - sum_j targets_j * d_i psi_j + d_i psi_4`.
The helper jointly recenters sources and targets before constructing these
moments to reduce origin-dependent cancellation. There are no wall images or
harmonic wall corrections.

`timing.backend_call_seconds` surrounds the `lfmm3d` invocation. Packing and
reconstruction have separate times; `total_seconds` includes input validation,
packing, backend invocation, checks and reconstruction. Overhead is total minus
backend-call time. These are MATLAB interface-call times, not kernel-only
profiles. Neither trees nor expansions are reused.

The opt-in comparison uses identical independent centered-cube sources and
targets for the two backends, warms both, alternates their order for seven
paired repeats, and reports medians of raw per-invocation times. Defaults:
`Ns=[512 4096]`, `Nt/Ns=[0.25 1 4]`, `tol=1e-9`. Direct references use all
sources at all targets for at most `2e7` pairs, otherwise up to 256 spatially
stratified targets. Errors include maximum pointwise absolute, global Frobenius
relative, and absolute error divided by sourcewise contribution magnitudes.
Accuracy checks are outside timed sections; no speedup is required to pass.
MAT/CSV outputs have unique filenames in `matlab/results` and include raw
samples, configuration, thread settings, resolved paths, source revision and
MEX SHA-256. Thread settings are recorded, never changed.

A several-fold speedup is plausible for the current upstream driver: it
internally requests Laplace source and target potentials, gradients and
Hessians even when only target Stokes velocity is requested. The new path
avoids source fields and Hessians. Both paths already use four Laplace
densities; four is not a fixed theoretical speedup. Performance depends on
source/target counts, tolerance, threading and the installed driver version.
