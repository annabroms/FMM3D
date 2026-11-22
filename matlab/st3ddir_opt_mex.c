/* Stokeslet/stresslet MEX gateway for st3ddir_opt. Built from the same
 * object as fmm3d. Inputs: sources, stoklet, strslet, strsvec, targ, flag, nd.
 * Pass [] for stoklet, or for strslet and strsvec, to omit that kernel. */
#include "mex.h"
#include <stdint.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
extern void st3ddirectstokg_opt_(int64_t*,double*,double*,int64_t*,
    double*,int64_t*,int64_t*,double*,double*,double*,double*);
extern void st3ddirectstrsg_opt_(int64_t*,double*,double*,double*,
    int64_t*,double*,int64_t*,int64_t*,double*,double*,double*,double*);
static int strength_ok(const mxArray *a,int64_t nd,int64_t ns)
{
    return mxGetM(a)==(mwSize)(3*nd) && mxGetN(a)==(mwSize)ns;
}
void mexFunction(int nlhs,mxArray **out,int nrhs,const mxArray **in)
{
    if(nrhs!=7 || nlhs<3 || nlhs>4)
        mexErrMsgIdAndTxt("st3ddir_opt:arguments","Require seven inputs and three/four outputs.");
    for(int k=0;k<7;++k)
        if(!mxIsDouble(in[k]) || mxIsComplex(in[k]) || mxIsSparse(in[k]) || mxGetNumberOfDimensions(in[k])!=2)
            mexErrMsgIdAndTxt("st3ddir_opt:type","Inputs must be full real double arrays.");
    if(mxGetNumberOfElements(in[5])!=1 || mxGetNumberOfElements(in[6])!=1)
        mexErrMsgIdAndTxt("st3ddir_opt:scalar","Flag and density count must be scalars.");
    double flag=mxGetScalar(in[5]), densities=mxGetScalar(in[6]);
    if(!isfinite(flag) || flag!=floor(flag) || flag<0 || flag>3 ||
       !isfinite(densities) || densities!=floor(densities) || densities<1 || densities>INT32_MAX/9)
        mexErrMsgIdAndTxt("st3ddir_opt:value","Invalid flag or density count.");
    int64_t nd=(int64_t)densities, ns=(int64_t)mxGetN(in[0]);
    int64_t nt=(int64_t)mxGetN(in[4]), ifppreg=(int64_t)flag;
    int ifstok=!mxIsEmpty(in[1]), ifstrs=!mxIsEmpty(in[2]);
    if(mxGetM(in[0])!=3 || mxGetM(in[4])!=3)
        mexErrMsgIdAndTxt("st3ddir_opt:shape","Invalid coordinate shape.");
    if((ifstok && !strength_ok(in[1],nd,ns)) ||
       (ifstrs && (!strength_ok(in[2],nd,ns) || !strength_ok(in[3],nd,ns))) ||
       (!ifstrs && !mxIsEmpty(in[3])))
        mexErrMsgIdAndTxt("st3ddir_opt:shape","Invalid strength shape.");
    out[0]=mxCreateDoubleMatrix((mwSize)(3*nd),(mwSize)nt,mxREAL);
    out[1]=mxCreateDoubleMatrix((mwSize)nd,(mwSize)nt,mxREAL);
    out[2]=mxCreateDoubleMatrix((mwSize)(9*nd),(mwSize)nt,mxREAL);
    double threshold=1e-15;
    /* Both routines add to the outputs, so the kernels can be combined. */
    if(ifstok)
        st3ddirectstokg_opt_(&nd,mxGetPr(in[0]),mxGetPr(in[1]),&ns,
            mxGetPr(in[4]),&nt,&ifppreg,mxGetPr(out[0]),mxGetPr(out[1]),
            mxGetPr(out[2]),&threshold);
    if(ifstrs)
        st3ddirectstrsg_opt_(&nd,mxGetPr(in[0]),mxGetPr(in[2]),mxGetPr(in[3]),
            &ns,mxGetPr(in[4]),&nt,&ifppreg,mxGetPr(out[0]),mxGetPr(out[1]),
            mxGetPr(out[2]),&threshold);
    if(nlhs==4) {
        int threads=1;
#ifdef _OPENMP
        threads=omp_get_max_threads();
#endif
        out[3]=mxCreateDoubleScalar(threads);
    }
}
