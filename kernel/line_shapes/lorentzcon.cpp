/* lorentzcon.cpp
 *
 * Lorentzian-triangle convolution MEX utility for Spinach
 *
 * Syntax:
 *
 *                    y=lorentzcon(offs,ampl,fwhm,x)
 */

#include "mex.h"
#include <cmath>
#include <cstdint>
#include <limits>

struct LinePlan
{
    int mode;
    double gam;
    double gam2;
    double centre;
    double edge0;
    double edge1;
    double offs0;
    double offs1;
    double offs2;
    double scale_lor;
    double scale_seg;
    double coef0;
    double coef1;
    double coef2;
    double coef_log0;
    double coef_log1;
};

static inline double real_value(const mxArray *array,mwIndex idx);
static inline double line_value0(double x,const LinePlan& plan);
static inline double line_value1(double x,const LinePlan& plan);
static inline double line_value2(double x,const LinePlan& plan);
static inline double line_value3(double x,const LinePlan& plan);
static inline double line_value4(double x,const LinePlan& plan);
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

        case 3:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value3((double)x_data[n],plan);
            return;

        default:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value4((double)x_data[n],plan);
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

        case 3:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value3(real_value(x_arg,n),plan);
            return;

        default:
            for (mwIndex n=0;n<n_elem;n++)
                y_data[n]=(YType)line_value4(real_value(x_arg,n),plan);
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
    const double gam=0.5*mxGetScalar(prhs[2]);
    const double pi=3.141592653589793238462643383279502884;
    const double two_ampl_pi=2.0*ampl/pi;
    const double ampl_gam_pi=ampl*gam/pi;

    // Read and sort offsets explicitly
    double offs[3]={0.0,0.0,0.0};

    for (mwIndex n=0;n<(mwIndex)n_offs;n++)
        offs[n]=real_value(offs_arg,n);

    if ((n_offs>=2)&&(offs[1]<offs[0]))
    {
        const double tmp=offs[0];
        offs[0]=offs[1];
        offs[1]=tmp;
    }

    if (n_offs==3)
    {
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
    plan.gam=gam;
    plan.gam2=gam*gam;
    plan.centre=offs[0];
    plan.edge0=0.0;
    plan.edge1=0.0;
    plan.offs0=offs[0];
    plan.offs1=offs[1];
    plan.offs2=offs[2];
    plan.scale_lor=ampl/(pi*gam);
    plan.scale_seg=0.0;
    plan.coef0=0.0;
    plan.coef1=0.0;
    plan.coef2=0.0;
    plan.coef_log0=0.0;
    plan.coef_log1=0.0;

    if (n_offs==2)
    {
        const double norm_offs=std::sqrt(offs[0]*offs[0]+offs[1]*offs[1]);
        const double sim_tol=std::sqrt(std::numeric_limits<double>::epsilon())*norm_offs;

        if (offs[1]-offs[0]<=sim_tol)
        {
            plan.mode=0;
            plan.centre=0.5*(offs[0]+offs[1]);
        }
        else
        {
            plan.mode=1;
            plan.edge0=offs[0];
            plan.edge1=offs[1];
            plan.scale_seg=ampl/(pi*(offs[1]-offs[0]));
        }
    }

    if (n_offs==3)
    {
        const double norm_offs=std::sqrt(offs[0]*offs[0]+offs[1]*offs[1]+offs[2]*offs[2]);
        const double sim_tol=std::sqrt(std::numeric_limits<double>::epsilon())*norm_offs;

        if (offs[2]-offs[0]<sim_tol)
        {
            plan.mode=0;
            plan.centre=(offs[0]+offs[1]+offs[2])/3.0;
        }
        else if (offs[1]-offs[0]<sim_tol)
        {
            plan.mode=2;
            plan.edge0=0.5*(offs[0]+offs[1]);
            plan.edge1=offs[2];
            const double diff=plan.edge0-plan.edge1;
            const double den=diff*diff;
            plan.coef0=two_ampl_pi/den;
            plan.coef_log0=ampl_gam_pi/den;
        }
        else if (offs[2]-offs[1]<sim_tol)
        {
            plan.mode=3;
            plan.edge0=offs[0];
            plan.edge1=0.5*(offs[1]+offs[2]);
            const double diff=plan.edge0-plan.edge1;
            const double den=diff*diff;
            plan.coef0=two_ampl_pi/den;
            plan.coef_log0=ampl_gam_pi/den;
        }
        else
        {
            plan.mode=4;
            const double den01=offs[0]-offs[1];
            const double den02=offs[0]-offs[2];
            const double den12=offs[1]-offs[2];
            const double den21=offs[2]-offs[1];
            plan.coef0=two_ampl_pi/(den01*den02);
            plan.coef1=two_ampl_pi/(den01*den12);
            plan.coef2=two_ampl_pi/(den02*den21);
            plan.coef_log0=ampl_gam_pi/(den01*den02);
            plan.coef_log1=ampl_gam_pi/(den02*den21);
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
        mexErrMsgIdAndTxt("Spinach:lorentzcon:alloc","Failed to allocate the output array.");

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
            mexErrMsgIdAndTxt("Spinach:lorentzcon:dataPtr","Failed to access array data.");
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
            mexErrMsgIdAndTxt("Spinach:lorentzcon:dataOut","Failed to access the output array data.");
        }

        if (mxGetClassID(x_arg)==mxDOUBLE_CLASS)
        {
            const double *x_data=mxGetDoubles(x_arg);

            if (x_data==nullptr)
            {
                mxDestroyArray(y_out);
                mexErrMsgIdAndTxt("Spinach:lorentzcon:dataPtr","Failed to access input array data.");
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
            mexErrMsgIdAndTxt("Spinach:lorentzcon:type","Inputs must be real numeric arrays.");
            return 0.0;
    }
}

static inline double line_value0(double x,const LinePlan& plan)
{
    const double arg=(x-plan.centre)/plan.gam;
    return plan.scale_lor/(1.0+arg*arg);
}

static inline double line_value1(double x,const LinePlan& plan)
{
    return plan.scale_seg*(std::atan2(plan.edge1-x,plan.gam)-
                           std::atan2(plan.edge0-x,plan.gam));
}

static inline double line_value2(double x,const LinePlan& plan)
{
    const double dx0=plan.edge0-x;
    const double dx1=plan.edge1-x;
    const double r0=dx0*dx0+plan.gam2;
    const double r1=dx1*dx1+plan.gam2;
    const double ang0=std::atan2(dx0,plan.gam);
    const double ang1=std::atan2(dx1,plan.gam);

    return plan.coef0*dx1*(ang1-ang0)+
           plan.coef_log0*std::log(r0/r1);
}

static inline double line_value3(double x,const LinePlan& plan)
{
    const double dx0=plan.edge0-x;
    const double dx1=plan.edge1-x;
    const double r0=dx0*dx0+plan.gam2;
    const double r1=dx1*dx1+plan.gam2;
    const double ang0=std::atan2(dx0,plan.gam);
    const double ang1=std::atan2(dx1,plan.gam);

    return plan.coef0*dx0*(ang0-ang1)+
           plan.coef_log0*std::log(r1/r0);
}

static inline double line_value4(double x,const LinePlan& plan)
{
    const double dx0=plan.offs0-x;
    const double dx1=plan.offs1-x;
    const double dx2=plan.offs2-x;
    const double r0=dx0*dx0+plan.gam2;
    const double r1=dx1*dx1+plan.gam2;
    const double r2=dx2*dx2+plan.gam2;

    return plan.coef0*dx0*std::atan2(dx0,plan.gam)-
           plan.coef1*dx1*std::atan2(dx1,plan.gam)-
           plan.coef2*dx2*std::atan2(dx2,plan.gam)+
           plan.coef_log0*std::log(r1/r0)+
           plan.coef_log1*std::log(r2/r1);
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=4)
        mexErrMsgIdAndTxt("Spinach:lorentzcon:nrhs","Four inputs are required: offs, ampl, fwhm, x.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:lorentzcon:nlhs","Exactly one output is required: y=lorentzcon(offs,ampl,fwhm,x).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:lorentzcon:prhsNull","Input argument array pointer is null.");

    const mxArray *offs=prhs[0];
    const mxArray *ampl=prhs[1];
    const mxArray *fwhm=prhs[2];
    const mxArray *x=prhs[3];

    if ((offs==nullptr)||(ampl==nullptr)||(fwhm==nullptr)||(x==nullptr))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:argNull","Input argument pointer is null.");

    if ((!mxIsNumeric(offs))||mxIsSparse(offs)||mxIsComplex(offs)||mxIsLogical(offs))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:offsType","offs must be a full real numeric array.");

    const mwSize n_offs=mxGetNumberOfElements(offs);

    if ((n_offs<1)||(n_offs>3))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:offsSize","offs must have one, two, or three elements.");

    for (mwIndex n=0;n<(mwIndex)n_offs;n++)
        if (!std::isfinite(real_value(offs,n)))
            mexErrMsgIdAndTxt("Spinach:lorentzcon:offsFinite","offs must contain finite real numbers.");

    if ((!mxIsNumeric(ampl))||mxIsSparse(ampl)||mxIsComplex(ampl)||mxIsLogical(ampl)||
        (mxGetNumberOfElements(ampl)!=1))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:amplType","ampl must be a finite real scalar.");

    const double ampl_val=mxGetScalar(ampl);

    if (!std::isfinite(ampl_val))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:amplFinite","ampl must be a finite real scalar.");

    if ((!mxIsNumeric(fwhm))||mxIsSparse(fwhm)||mxIsComplex(fwhm)||mxIsLogical(fwhm)||
        (mxGetNumberOfElements(fwhm)!=1))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:fwhmType","fwhm must be a finite positive real scalar.");

    const double fwhm_val=mxGetScalar(fwhm);

    if ((!std::isfinite(fwhm_val))||(!(fwhm_val>0.0)))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:fwhmFinite","fwhm must be a finite positive real scalar.");

    if ((!mxIsNumeric(x))||mxIsSparse(x)||mxIsComplex(x)||mxIsLogical(x))
        mexErrMsgIdAndTxt("Spinach:lorentzcon:xType","x must be a full real numeric array.");

    const mwSize n_elem=mxGetNumberOfElements(x);

    for (mwIndex n=0;n<n_elem;n++)
        if (!std::isfinite(real_value(x,n)))
            mexErrMsgIdAndTxt("Spinach:lorentzcon:xFinite","x must contain finite real numbers.");
}
