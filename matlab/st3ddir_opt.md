# Optimized direct Stokeslet and stresslet sums (`st3ddir_opt`)

`st3ddir_opt` evaluates free-space Stokeslet and (type I) stresslet velocities
by direct O(ns·nt) summation. It has the same inputs, outputs, kernels and
scaling as `st3ddir`, and uses two optimized Fortran routines in
[`src/Stokes/stokkernels.f`](../src/Stokes/stokkernels.f).

**Results** (measured 2026-11-22, velocity only, N = 10^3 to 2·10^4 sources
and targets; full tables below):

| Build | Stokeslet: `st3ddir_opt` vs `st3ddir` | Stresslet: `st3ddir_opt` vs `st3ddir` |
|---|---:|---:|
| Old flags (generic x86-64/SSE4, no fast-math) | 1.7–2.0x | 2.7–3.2x |
| New flags (AVX2 + fast-math, `stokkernels.f` only) | 4.8–6.3x | 9.1–10.6x |

For N ≥ 1,000, the new flags make `st3ddir_opt` itself 2.2–3.8x faster than
the old flags. The relative difference between `st3ddir_opt` and `st3ddir` is
at most 1.1e-14.

## What was changed

### Fortran kernels in `src/Stokes/stokkernels.f`

Both optimized routines have the argument list of the corresponding original
routine and keep its semantics:

- **Outputs are additive.** Results are added to `pot`, `pre` and `grad`.
- **Only velocity is optimized** (`ifppreg = 1`). For `ifppreg = 2` or `3`,
  pressure and gradient, the routines call the original Fortran routine, so
  those results are bit-for-bit unchanged.
- **Same exclusion rule:** pairs with `r^2 < thresh^2` are skipped and no
  self-term is added.
- **Same kernels, signs and scaling** as the definitions at the top of the
  file.

**`st3ddirectstokg_opt`** (Stokeslet; commit `f0cf435`)

1. Source coordinates and force components are packed once per call into
   separate contiguous arrays.
2. One reciprocal distance `ri = 1/sqrt(r2)` per pair is reused for all three
   velocity components, which replaces repeated divisions.
3. The sum over sources is accumulated in three scalars (`nd = 1`) or three
   length-`nd` arrays (`nd > 1`) instead of updating `pot` for every pair.
   `1/(8 pi)` is applied once per target.
4. For `nd = 1`, `!$omp simd reduction` vectorizes the loop over sources (the
   summation order may change). For `nd > 1`, densities are contiguous and
   vectorized while the geometry is shared between them.
5. OpenMP runs over targets with static scheduling, only when
   `nt*ns >= 4096`. Inside an existing parallel region no new thread team is
   created.

**`st3ddirectstrsg_opt`** (stresslet only; new)

The original stresslet sum is slower than the Stokeslet sum for two reasons.
`st3ddir` calls `st3ddirectstokstrsrotdoubg` → `st3ddirectstokstrsg`, which
always first runs a full Stokeslet loop. For a stresslet-only call that loop
adds a zero Stokeslet field. The stresslet loop itself also reads strided
memory, updates `pot` for every pair, and divides by `r^5` for every pair.

The new routine evaluates
`pot_i += 3/(4 pi) · r_i (r·mu)(r·nu) / r^5` with
`mu = strslet`, `nu = strsvec`:

1. It takes no Stokeslet input, so there is no zero Stokeslet pass.
2. Sources, `mu` and `nu` are packed once per call into contiguous arrays.
3. One `ri = 1/sqrt(r2)` per pair gives `ri^5` by multiplication. The scalar
   `s = (r·mu)(r·nu)·ri^5` multiplies all three components.
4. Accumulation, the `nd = 1` / `nd > 1` vectorization and the OpenMP
   strategy are as in `st3ddirectstokg_opt`. `3/(4 pi)` is applied once per
   target.
5. For `ifppreg = 2` or `3` it calls `st3ddirectstokstrsg` with zero Stokeslet
   strengths.

The original routines are unchanged by this work. Commit `f0cf435` had
already added the `ifppreg` output flag and OpenMP over targets to them.

### MATLAB interface

- **[`st3ddir_opt.m`](st3ddir_opt.m)** accepts `srcinfo.stoklet`,
  `srcinfo.strslet` together with `srcinfo.strsvec`, or both. If both are
  given, the contributions are added. Rotlets and doublets give an error; use
  `st3ddir` for those. The call and outputs are as before:

  ```matlab
  srcinfo.sources = src;          % 3 x ns
  srcinfo.strslet = mu;           % 3 x ns (or nd x 3 x ns, with srcinfo.nd)
  srcinfo.strsvec = nu;           % 3 x ns
  [U, threads] = st3ddir_opt(srcinfo, targ, 1);   % U.pottarg: 3 x nt
  ```

  Flag 1 is the optimized velocity-only path. The default flag is 3, as in
  `st3ddir`. `threads` is the number of OpenMP threads available. In MATLAB
  this follows `maxNumCompThreads`; `OMP_NUM_THREADS` has no effect.
- **[`st3ddir_opt_mex.c`](st3ddir_opt_mex.c)** is the gateway. It now takes
  `(sources, stoklet, strslet, strsvec, targ, flag, nd)`, where an empty
  array omits that kernel, and calls one or both optimized routines. Only
  `st3ddir_opt.m` calls it.

### Build system

The [`makefile`](../makefile) applies extra flags to `stokkernels.o` only:

```make
STOKDIRFLAGS ?= -ffast-math
$(STOK)/stokkernels.o: FFLAGS += $(STOKDIRFLAGS)
```

`stfmm3d` does not call any routine in `stokkernels.f`; its near field uses
the Laplace kernels. These flags therefore change only the direct Stokes sums
(`st3ddir`, `st3ddir_opt`), not the FMM. Architecture flags are not set in the
makefile. They come from `FFLAGS` in `make.inc`, which may also set
`STOKDIRFLAGS` (see below). `make STOKDIRFLAGS=` disables the extra flags.

### Tests

[`test/Stokes/test_stokkernels_opt.f90`](../test/Stokes/test_stokkernels_opt.f90)
(`make test`) checks both optimized routines against the originals:

- nd = 1, 2, 3, 4, 8, 16
- flags 0–3, with bit-for-bit equality for the fallback flags and a relative
  budget of 3e-12 for velocity
- prefilled outputs
- coincident, near-threshold and exact-threshold pairs
- empty sources, all pairs excluded, and zero strengths
- 1 and 8 threads, and calls from inside a parallel region

It passed with a bounds-checked strict build
(`gfortran -O2 -fopenmp -fcheck=all`), in `--serial` mode, and with the
AVX2 + fast-math library objects.

In MATLAB, `st3ddir_opt` was compared with `st3ddir` for stresslets alone and
for Stokeslets plus stresslets, nd = 1, 2, 3, 8 and flags 1–3. The maximum
relative difference was 6.2e-16. A target on a source and a target 1e-16 from
a source were included.

The Python (f2py/meson) build also wraps `stokkernels.f`, so it will expose
`st3ddirectstrsg_opt`. It has not been built or tested.

## How the MEX files were compiled

Both builds used gfortran/gcc 11.4 and MATLAB R2025a on Linux. The general
flags in `make.inc` targeted generic x86-64 with SSE4:

```make
FFLAGS= -fPIC -O3 -funroll-loops -march=x86-64 -mtune=generic -msse4 -fcx-limited-range -std=legacy -w
```

**Old flags.** `stokkernels.f` is compiled like the rest of the library:

```sh
rm -f src/Stokes/stokkernels.o
make matlab STOKDIRFLAGS=
# gfortran -c -fPIC -O3 -funroll-loops -march=x86-64 -mtune=generic -msse4
#   -fcx-limited-range -std=legacy -w -fopenmp ... src/Stokes/stokkernels.f
```

**New flags.** `make.inc` additionally contains

```make
STOKDIRFLAGS = -march=x86-64-v3 -ffast-math
```

and the build is

```sh
rm -f src/Stokes/stokkernels.o
make matlab
# gfortran -c ... -msse4 ... -fopenmp -march=x86-64-v3 -ffast-math ... src/Stokes/stokkernels.f
```

The later `-march=x86-64-v3` overrides the generic `-march=x86-64` for this
file only. All other objects keep the generic flags.

In both builds the MEX gateways are compiled and linked by `make matlab` with

```
mex matlab/st3ddir_opt_mex.c lib-static/libfmm3d.a -compatibleArrayDims -DMWF77_UNDERSCORE1
    "CFLAGS=-std=gnu17 -Wno-implicit-function-declaration -fPIC -O3 -funroll-loops"
    -D_OPENMP -lm -lstdc++ -ldl -lgfortran -lgomp
```

`-ffast-math` is not passed at link time. As a result `crtfastmath.o` is not
linked, and MATLAB's floating-point mode (flush-to-zero) is not changed.

The disassembly confirms the difference. With the old flags both optimized
loops use only SSE registers, with scalar `sqrtsd`/`divsd`. With the new
flags they use 256-bit AVX2 registers with packed `sqrtpd`/`divpd`. The
`sqrt` and division are vectorized only when AVX2 and `-ffast-math` are
combined. In a separate run with the Stokeslet routine (1 thread), AVX2
without fast-math gave 1.19x over the old flags, and fast-math with SSE4 gave
no gain.

**Portability.** Binaries built with `-march=x86-64-v3` (or `-march=native`)
need a CPU with AVX2 at run time. On older or low-end x86 CPUs without it
(roughly pre-2013 Intel, pre-2015 AMD), MATLAB stops with "illegal
instruction". The flag does not compile for non-x86 targets. Set
`STOKDIRFLAGS` in `make.inc` only on machines where the binaries will run.
The makefile default (`-ffast-math` alone) is portable.

## Benchmark setup

Velocity only (flag 1), ns = nt = N:
- Sources and forces/strengths/orientations are standard normal (`randn`) and
  targets are `randn + 0.01`, with seed N.
- `st3ddir` (the original routines through `fmm3d`) and `st3ddir_opt` are both
  timed through their MATLAB wrappers. `st3ddir` was not timed at N = 20,000.
- Times are the minimum over 5 repetitions with rotated execution order. Small
  N are repeated internally.
- There is one MATLAB session per build. Threads are set with
  `maxNumCompThreads` (8 or 1).

The CPU was an 8-core x86-64 processor (16 hardware threads, AVX-512).
Absolute times and speedups depend on the CPU. In repeated sessions, 8-thread
times of the same binary varied by up to about 1.5x; single-thread times were
stable to a few percent.

## Results

Times are in ms. Speedups are `st3ddir` time ÷ `st3ddir_opt` time within the
same build. "Old"/"new" refer to the compiler flags for `stokkernels.f`.

### Stokeslet

| Threads | N | `st3ddir` old | `st3ddir` new | `st3ddir_opt` old | `st3ddir_opt` new | Speedup, old flags | Speedup, new flags | `st3ddir_opt` new vs old flags |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 100 | 0.044 | 0.030 | 0.023 | 0.016 | 1.87 | 1.87 | 1.46 |
| 8 | 1,000 | 2.08 | 1.86 | 1.24 | 0.339 | 1.67 | 5.47 | 3.67 |
| 8 | 5,000 | 30.8 | 28.7 | 15.9 | 5.57 | 1.93 | 5.16 | 2.86 |
| 8 | 10,000 | 145 | 109 | 74.0 | 22.7 | 1.97 | 4.79 | 3.26 |
| 8 | 20,000 | – | – | 244 | 109 | – | – | 2.25 |
| 1 | 1,000 | 7.53 | 6.58 | 3.99 | 1.06 | 1.89 | 6.21 | 3.76 |
| 1 | 10,000 | 753 | 679 | 391 | 107 | 1.93 | 6.34 | 3.65 |

### Stresslet

| Threads | N | `st3ddir` old | `st3ddir` new | `st3ddir_opt` old | `st3ddir_opt` new | Speedup, old flags | Speedup, new flags | `st3ddir_opt` new vs old flags |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 100 | 0.060 | 0.058 | 0.032 | 0.022 | 1.91 | 2.64 | 1.44 |
| 8 | 1,000 | 3.57 | 3.70 | 1.32 | 0.407 | 2.70 | 9.09 | 3.25 |
| 8 | 5,000 | 57.4 | 49.5 | 17.7 | 5.42 | 3.24 | 9.13 | 3.28 |
| 8 | 10,000 | 216 | 201 | 70.5 | 21.0 | 3.07 | 9.58 | 3.36 |
| 8 | 20,000 | – | – | 295 | 126 | – | – | 2.34 |
| 1 | 1,000 | 12.7 | 12.6 | 4.02 | 1.19 | 3.17 | 10.59 | 3.38 |
| 1 | 10,000 | 1,280 | 1,220 | 412 | 124 | 3.10 | 9.83 | 3.32 |

### Accuracy

Relative difference `norm(u_opt - u_st3ddir)/norm(u_st3ddir)` over all
targets, for the same inputs (8 threads):

| N | Stokeslet, old flags | Stokeslet, new flags | Stresslet, old flags | Stresslet, new flags |
|---:|---:|---:|---:|---:|
| 100 | 6.9e-16 | 6.4e-16 | 6.7e-16 | 5.9e-16 |
| 1,000 | 1.9e-15 | 1.5e-15 | 1.7e-15 | 5.6e-16 |
| 5,000 | 6.2e-15 | 5.2e-15 | 3.3e-15 | 2.8e-15 |
| 10,000 | 7.0e-15 | 6.3e-15 | 4.3e-15 | 3.2e-15 |
| 20,000 | 1.1e-14 | 8.2e-15 | 2.8e-15 | 2.3e-15 |

### Observations

- **With the old flags,** the optimized routines are 1.7–2.0x (Stokeslet) and
  2.7–3.2x (stresslet) faster than `st3ddir`. Part of the stresslet gain comes
  from skipping the zero Stokeslet pass.
- **AVX2 + fast-math on `stokkernels.f`** makes `st3ddir_opt` 2.2–3.8x faster
  than with the old flags for N ≥ 1,000. The original `st3ddir` changes by
  only 1.0–1.3x for N ≥ 1,000, because its loops do not vectorize. With these
  flags, `st3ddir_opt` is 5–6x (Stokeslet) and 9–11x (stresslet) faster than
  `st3ddir`.
- **For N ≈ 100,** a call to `st3ddir_opt` takes 16–32 µs, much of it the
  input checks in `st3ddir_opt.m`. It is still 1.9–2.6x faster than
  `st3ddir`.
- **Accuracy:** fast-math does not reduce accuracy measurably here. The
  differences are at the level of rounding and are no larger with the new
  flags.
