To install in the system:
pip install -e . 


On compilation you can import the package fmm3dpy
which has 2 main subroutines

- hfmm3d: for computing N-body Helmholtz interactions
- lfmm3d: for computing N-body Laplace interactions

Examples for both of these are demonstrated in
``hfmmexample.py`` and ``lfmmexample.py``

## Direct Stokes evaluation

Direct Stokes evaluation supports OpenMP. The Meson Python build enables it
automatically when the Fortran compiler supports it; the Makefile build uses
its configured OpenMP flags. Set `OMP_NUM_THREADS` before importing Python:

```bash
OMP_NUM_THREADS=8 python your_script.py
```

`st3ddir(..., ifppreg=0, ifppregtarg=1)` evaluates target velocity only.
Independent targets use static scheduling and retain source summation order.
Stokeslet, stresslet, rotlet and doublet paths support output flags 0--3.
Interactions with fewer than 4096 source-target pairs remain serial. Calls
inside an existing OpenMP region do not create nested worker teams. Builds
without OpenMP remain usable as serial code. Rebuild the Python bindings
after changing the Fortran source.

## Optimized Stokeslet velocity

The additional `st3ddirectstokg_opt` symbol optimizes Stokeslet velocity for
one or multiple densities. After rebuilding the bindings:

```python
from fmm3dpy import stfmm3d_fortran as direct

# sources: (3, ns), strengths: (nd, 3, ns), targets: (3, nt)
# strengths = force / viscosity; the kernel includes 1/(8*pi).
velocity, pressure, gradient = direct.st3ddirectstokg_opt(
    sources, strengths, targets, 1, 1e-16)
```

The signature and additive Fortran output semantics match `st3ddirectstokg`.
Output flag 1 uses the optimized velocity path for every positive density
count. Flags 0, 2 and 3 call the original implementation. The f2py interface
allocates pressure/gradient arrays even for flag 1; those arrays remain zero.
The single-density path uses component packing and a SIMD source reduction;
the multiple-density path shares distances and vectorizes across densities.
Packing and worker-private accumulators are allocated internally. Floating
point sums can differ by rounding. Both target paths use OpenMP, without
nested worker teams. `st3ddir` and the FMM retain their existing dispatch;
using `_opt` is explicit.

The regression test `test/Stokes/test_stokkernels_opt.f90` runs with
`make test` unless `OMP=OFF`. It checks multiple densities, all flags,
prefilled outputs, thresholds and parallel execution.
