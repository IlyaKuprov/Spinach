/* gausscon.cpp
 *
 * Gaussian-triangle convolution MEX utility for Spinach
 *
 * Syntax:
 *
 *                    y=gausscon(offs,ampl,fwhm,x)
 */

#include "mex.h"
#include <cmath>
#include <cstdint>
#include <limits>

struct LinePlan
{
    int mode;
    double sigma;
    double sigma2;
    double inv_sqrt2_sigma;
    double centre;
    double edge0;
    double edge1;
    double offs0;
    double offs1;
    double offs2;
    double ampl;
    double norm_gauss;
    double coef0;
    double coef1;
};

static inline double real_value(const mxArray *array,mwIndex idx);
static inline double gauss_value(double diff,const LinePlan& plan);
static inline double prob_int(double edge0,double edge1,double x,const LinePlan& plan);
static inline double line_value0(double x,const LinePlan& plan);
static inline double line_value1(double x,const LinePlan& plan);
static inline double line_value2(double x,const LinePlan& plan);
static inline double line_value3(double x,const LinePlan& plan);
template<typename XType,typename YType>
static void fill_direct(const XType *x_data,YType *y_data,mwSize n_elem,const LinePlan& plan);
template<typename YType>
static void fill_real(const mxArray *x_arg,YType *y_data,mwSize n_elem,const LinePlan& plan);
static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);

template<typename XType,typename YType>
static void fill_direct(const XType *x_data,YType *y_data,mwSize n_elem,const LinePlan& plan)
{
    switch (plan.mode)
    {
        case 0:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value0((double)x_data[n],plan);
            return;

        case 1:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value1((double)x_data[n],plan);
            return;

        case 2:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value2((double)x_data[n],plan);
            return;

        default:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value3((double)x_data[n],plan);
            return;
    }
}

template<typename YType>
static void fill_real(const mxArray *x_arg,YType *y_data,mwSize n_elem,const LinePlan& plan)
{
    switch (plan.mode)
    {
        case 0:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value0(real_value(x_arg,n),plan);
            return;

        case 1:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value1(real_value(x_arg,n),plan);
            return;

        case 2:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value2(real_value(x_arg,n),plan);
            return;

        default:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value3(real_value(x_arg,n),plan);
            return;
    }
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

    // Get input arguments
    const mxArray *offs_arg=prhs[0];
    const mxArray *x_arg=prhs[3];
    const mwSize n_offs=mxGetNumberOfElements(offs_arg);
    const double ampl=mxGetScalar(prhs[1]);
    const double fwhm=mxGetScalar(prhs[2]);
    const double pi=3.141592653589793238462643383279502884;
    const double ln2=0.693147180559945309417232121458176568;
    const double sqrt2=1.414213562373095048801688724209698079;
    const double sigma=fwhm/(2.0*std::sqrt(2.0*ln2));
    const double sigma2=sigma*sigma;

    // Read and sort offsets explicitly
    double offs[3]={0.0,0.0,0.0};

    for (mwIndex n=0;n<(mwIndex)n_offs;n++)
        offs[n]=real_value(offs_arg,n);

    if (n_offs==3)
    {
        if (offs[1]<offs[0])
        {
            const double tmp=offs[0];
            offs[0]=offs[1];
            offs[1]=tmp;
        }

        if (offs[2]<offs[1])
        {
            const double tmp=offs[1];
            offs[1]=offs[2];
            offs[2]=tmp;
        }

        if (offs[1]<offs[0])
        {
            const double tmp=offs[0];
            offs[0]=offs[1];
            offs[1]=tmp;
        }
    }

    // Precompute branch constants
    LinePlan plan;
    plan.mode=0;
    plan.sigma=sigma;
    plan.sigma2=sigma2;
    plan.inv_sqrt2_sigma=1.0/(sqrt2*sigma);
    plan.centre=offs[0];
    plan.edge0=0.0;
    plan.edge1=0.0;
    plan.offs0=offs[0];
    plan.offs1=offs[1];
    plan.offs2=offs[2];
    plan.ampl=ampl;
    plan.norm_gauss=1.0/(sigma*std::sqrt(2.0*pi));
    plan.coef0=0.0;
    plan.coef1=0.0;

    if (n_offs==3)
    {
        double sim_tol=std::sqrt(std::numeric_limits<double>::epsilon())*
                       std::sqrt(offs[0]*offs[0]+offs[1]*offs[1]+offs[2]*offs[2]);

        if (sim_tol==0.0)
            sim_tol=std::numeric_limits<double>::epsilon();

        if (offs[2]-offs[0]<=sim_tol)
        {
            plan.mode=0;
            plan.centre=(offs[0]+offs[1]+offs[2])/3.0;
        }
        else if (offs[1]-offs[0]<=sim_tol)
        {
            plan.mode=1;
            plan.edge0=0.5*(offs[0]+offs[1]);
            plan.edge1=offs[2];
            const double diff=plan.edge1-plan.edge0;
            plan.coef0=2.0*ampl/(diff*diff);
        }
        else if (offs[2]-offs[1]<=sim_tol)
        {
            plan.mode=2;
            plan.edge0=offs[0];
            plan.edge1=0.5*(offs[1]+offs[2]);
            const double diff=plan.edge1-plan.edge0;
            plan.coef0=2.0*ampl/(diff*diff);
        }
        else
        {
            plan.mode=3;
            plan.coef0=2.0*ampl/((offs[1]-offs[0])*(offs[2]-offs[0]));
            plan.coef1=2.0*ampl/((offs[2]-offs[1])*(offs[2]-offs[0]));
        }
    }

    // Get input array dimensions
    const mwSize n_dims=mxGetNumberOfDimensions(x_arg);
    const mwSize *dims=mxGetDimensions(x_arg);
    const mwSize n_elem=mxGetNumberOfElements(x_arg);

    // Allocate output array
    const bool single_out=(mxGetClassID(x_arg)==mxSINGLE_CLASS);
    mxArray *y_out=mxCreateUninitNumericArray(n_dims,const_cast<mwSize*>(dims),
                                              single_out?mxSINGLE_CLASS:mxDOUBLE_CLASS,mxREAL);

    if (y_out==nullptr)
        mexErrMsgIdAndTxt("Spinach:gausscon:alloc","Failed to allocate the output array.");

    // Return empty arrays without resolving a data pointer
    if (n_elem==0)
    {
        plhs[0]=y_out;
        return;
    }

    // Fill single-precision output array
    if (single_out)
    {
        float *y_data=mxGetSingles(y_out);
        const float *x_data=mxGetSingles(x_arg);

        if ((y_data==nullptr)||(x_data==nullptr))
        {
            mxDestroyArray(y_out);
            mexErrMsgIdAndTxt("Spinach:gausscon:dataPtr","Failed to access array data.");
        }

        // Evaluate the selected branch in a single pass
        fill_direct(x_data,y_data,n_elem,plan);
    }

    // Fill double-precision output array
    else
    {
        double *y_data=mxGetDoubles(y_out);

        if (y_data==nullptr)
        {
            mxDestroyArray(y_out);
            mexErrMsgIdAndTxt("Spinach:gausscon:dataOut","Failed to access the output array data.");
        }

        if (mxGetClassID(x_arg)==mxDOUBLE_CLASS)
        {
            const double *x_data=mxGetDoubles(x_arg);

            if (x_data==nullptr)
            {
                mxDestroyArray(y_out);
                mexErrMsgIdAndTxt("Spinach:gausscon:dataPtr","Failed to access input array data.");
            }

            // Evaluate the selected branch in a single pass
            fill_direct(x_data,y_data,n_elem,plan);
        }
        else
        {
            // Evaluate non-floating-point input through scalar conversion
            fill_real(x_arg,y_data,n_elem,plan);
        }
    }

    // Return output argument
    plhs[0]=y_out;
}

static inline double real_value(const mxArray *array,mwIndex idx)
{
    switch (mxGetClassID(array))
    {
        case mxDOUBLE_CLASS:
            return mxGetDoubles(array)[idx];

        case mxSINGLE_CLASS:
            return (double)mxGetSingles(array)[idx];

        case mxINT8_CLASS:
            return (double)((int8_t*)mxGetData(array))[idx];

        case mxUINT8_CLASS:
            return (double)((uint8_t*)mxGetData(array))[idx];

        case mxINT16_CLASS:
            return (double)((int16_t*)mxGetData(array))[idx];

        case mxUINT16_CLASS:
            return (double)((uint16_t*)mxGetData(array))[idx];

        case mxINT32_CLASS:
            return (double)((int32_t*)mxGetData(array))[idx];

        case mxUINT32_CLASS:
            return (double)((uint32_t*)mxGetData(array))[idx];

        case mxINT64_CLASS:
            return (double)((int64_t*)mxGetData(array))[idx];

        case mxUINT64_CLASS:
            return (double)((uint64_t*)mxGetData(array))[idx];

        default:
            mexErrMsgIdAndTxt("Spinach:gausscon:type","Inputs must be real numeric arrays.");
            return 0.0;
    }
}

static inline double gauss_value(double diff,const LinePlan& plan)
{
    const double arg=diff/plan.sigma;
    return plan.norm_gauss*std::exp(-0.5*arg*arg);
}

static inline double prob_int(double edge0,double edge1,double x,const LinePlan& plan)
{
    return 0.5*(std::erf((edge1-x)*plan.inv_sqrt2_sigma)-
                std::erf((edge0-x)*plan.inv_sqrt2_sigma));
}

static inline double line_value0(double x,const LinePlan& plan)
{
    return plan.ampl*gauss_value(x-plan.centre,plan);
}

static inline double line_value1(double x,const LinePlan& plan)
{
    const double g0=gauss_value(x-plan.edge0,plan);
    const double g1=gauss_value(x-plan.edge1,plan);
    const double prob=prob_int(plan.edge0,plan.edge1,x,plan);

    return plan.coef0*((plan.edge1-x)*prob+plan.sigma2*(g1-g0));
}

static inline double line_value2(double x,const LinePlan& plan)
{
    const double g0=gauss_value(x-plan.edge0,plan);
    const double g1=gauss_value(x-plan.edge1,plan);
    const double prob=prob_int(plan.edge0,plan.edge1,x,plan);

    return plan.coef0*((x-plan.edge0)*prob-plan.sigma2*(g1-g0));
}

static inline double line_value3(double x,const LinePlan& plan)
{
    const double g0=gauss_value(x-plan.offs0,plan);
    const double g1=gauss_value(x-plan.offs1,plan);
    const double g2=gauss_value(x-plan.offs2,plan);
    const double prob01=prob_int(plan.offs0,plan.offs1,x,plan);
    const double prob12=prob_int(plan.offs1,plan.offs2,x,plan);
    const double left_int=(x-plan.offs0)*prob01-plan.sigma2*(g1-g0);
    const double right_int=(plan.offs2-x)*prob12+plan.sigma2*(g2-g1);

    return plan.coef0*left_int+plan.coef1*right_int;
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=4)
        mexErrMsgIdAndTxt("Spinach:gausscon:nrhs","Four inputs are required: offs, ampl, fwhm, x.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:gausscon:nlhs","Exactly one output is required: y=gausscon(offs,ampl,fwhm,x).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:gausscon:prhsNull","Input argument array pointer is null.");

    const mxArray *offs=prhs[0];
    const mxArray *ampl=prhs[1];
    const mxArray *fwhm=prhs[2];
    const mxArray *x=prhs[3];

    if ((offs==nullptr)||(ampl==nullptr)||(fwhm==nullptr)||(x==nullptr))
        mexErrMsgIdAndTxt("Spinach:gausscon:argNull","Input argument pointer is null.");

    if ((!mxIsNumeric(offs))||mxIsSparse(offs)||mxIsComplex(offs)||mxIsLogical(offs))
        mexErrMsgIdAndTxt("Spinach:gausscon:offsType","offs must be a full real numeric array.");

    const mwSize n_offs=mxGetNumberOfElements(offs);

    if ((n_offs!=1)&&(n_offs!=3))
        mexErrMsgIdAndTxt("Spinach:gausscon:offsSize","offs must have one or three elements.");

    for (mwIndex n=0;n<(mwIndex)n_offs;n++)
        if (!std::isfinite(real_value(offs,n)))
            mexErrMsgIdAndTxt("Spinach:gausscon:offsFinite","offs must contain finite real numbers.");

    if ((!mxIsNumeric(ampl))||mxIsSparse(ampl)||mxIsComplex(ampl)||mxIsLogical(ampl)||
        (mxGetNumberOfElements(ampl)!=1))
        mexErrMsgIdAndTxt("Spinach:gausscon:amplType","ampl must be a finite real scalar.");

    const double ampl_val=mxGetScalar(ampl);

    if (!std::isfinite(ampl_val))
        mexErrMsgIdAndTxt("Spinach:gausscon:amplFinite","ampl must be a finite real scalar.");

    if ((!mxIsNumeric(fwhm))||mxIsSparse(fwhm)||mxIsComplex(fwhm)||mxIsLogical(fwhm)||
        (mxGetNumberOfElements(fwhm)!=1))
        mexErrMsgIdAndTxt("Spinach:gausscon:fwhmType","fwhm must be a finite positive real scalar.");

    const double fwhm_val=mxGetScalar(fwhm);

    if ((!std::isfinite(fwhm_val))||(!(fwhm_val>0.0)))
        mexErrMsgIdAndTxt("Spinach:gausscon:fwhmFinite","fwhm must be a finite positive real scalar.");

    if ((!mxIsNumeric(x))||mxIsSparse(x)||mxIsComplex(x)||mxIsLogical(x))
        mexErrMsgIdAndTxt("Spinach:gausscon:xType","x must be a full real numeric array.");

    const mwSize n_elem=mxGetNumberOfElements(x);

    for (mwIndex n=0;n<n_elem;n++)
        if (!std::isfinite(real_value(x,n)))
            mexErrMsgIdAndTxt("Spinach:gausscon:xFinite","x must contain finite real numbers.");
}
