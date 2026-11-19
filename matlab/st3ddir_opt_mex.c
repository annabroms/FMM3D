/* Stokeslet-only MEX gateway. Both direct versions are built from one object. */
#include "mex.h"
#include <stdint.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
extern void st3ddirectstokg_opt_(int64_t*,double*,double*,int64_t*,
    double*,int64_t*,int64_t*,double*,double*,double*,double*);
void mexFunction(int nlhs,mxArray **out,int nrhs,const mxArray **in)
{
    if(nrhs!=5 || nlhs<3 || nlhs>4)
        mexErrMsgIdAndTxt("st3ddir_opt:arguments","Require five inputs and three/four outputs.");
    for(int k=0;k<5;++k)
        if(!mxIsDouble(in[k]) || mxIsComplex(in[k]) || mxIsSparse(in[k]) || mxGetNumberOfDimensions(in[k])!=2)
            mexErrMsgIdAndTxt("st3ddir_opt:type","Inputs must be full real double arrays.");
    if(mxGetNumberOfElements(in[3])!=1 || mxGetNumberOfElements(in[4])!=1)
        mexErrMsgIdAndTxt("st3ddir_opt:scalar","Flag and density count must be scalars.");
    double flag=mxGetScalar(in[3]), densities=mxGetScalar(in[4]);
    if(!isfinite(flag) || flag!=floor(flag) || flag<0 || flag>3 ||
       !isfinite(densities) || densities!=floor(densities) || densities<1 || densities>INT32_MAX/9)
        mexErrMsgIdAndTxt("st3ddir_opt:value","Invalid flag or density count.");
    int64_t nd=(int64_t)densities, ns=(int64_t)mxGetN(in[0]);
    int64_t nt=(int64_t)mxGetN(in[2]), ifppreg=(int64_t)flag;
    if(mxGetM(in[0])!=3 || mxGetM(in[2])!=3 ||
       mxGetM(in[1])!=(mwSize)(3*nd) || mxGetN(in[1])!=(mwSize)ns)
        mexErrMsgIdAndTxt("st3ddir_opt:shape","Invalid coordinate or strength shape.");
    out[0]=mxCreateDoubleMatrix((mwSize)(3*nd),(mwSize)nt,mxREAL);
    out[1]=mxCreateDoubleMatrix((mwSize)nd,(mwSize)nt,mxREAL);
    out[2]=mxCreateDoubleMatrix((mwSize)(9*nd),(mwSize)nt,mxREAL);
    double threshold=1e-15;
    st3ddirectstokg_opt_(&nd,mxGetPr(in[0]),mxGetPr(in[1]),&ns,
        mxGetPr(in[2]),&nt,&ifppreg,mxGetPr(out[0]),mxGetPr(out[1]),
        mxGetPr(out[2]),&threshold);
    if(nlhs==4) {
        int threads=1;
#ifdef _OPENMP
        threads=omp_get_max_threads();
#endif
        out[3]=mxCreateDoubleScalar(threads);
    }
}
