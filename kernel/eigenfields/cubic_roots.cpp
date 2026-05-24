/* cubic_roots.cpp
 *
 * Cubic-polynomial real-root filter MEX utility for Spinach
 *
 * Syntax:
 *
 *                    root_list=cubic_roots(poly_coeffs,root_tol)
 *
 * Returns sorted real roots in the interval [0,1] of
 *
 *                    y=a*x^3+b*x^2+c*x+d
 *
 * where poly_coeffs=[a b c d]. Leading numerical zeros are dropped
 * relative to the largest coefficient.
 */

#include "mex.h"
#include <cmath>
#include <limits>

static const double pi=3.141592653589793238462643383279502884;

static inline double clamp_value(const double x,const double lo,const double hi);
static inline void append_root(double *roots,mwSize& n_roots,double root,const double root_tol);
static inline void sort_roots(double *roots,const mwSize n_roots);
static mwSize solve_linear(const double c,const double d,double *roots,const double root_tol);
static mwSize solve_quadratic(const double b2,const double c2,const double d2,double *roots,const double root_tol);
static mwSize solve_cubic(const double a3,const double b3,const double c3,const double d3,double *roots,const double root_tol);
static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);

static inline double clamp_value(const double x,const double lo,const double hi)
{
    if (x<lo)
        return lo;

    if (x>hi)
        return hi;

    return x;
}

static inline void append_root(double *roots,mwSize& n_roots,double root,const double root_tol)
{
    if (!std::isfinite(root))
        return;

    if ((root<-root_tol)||(root>1.0+root_tol))
        return;

    root=clamp_value(root,0.0,1.0);

    for (mwSize n=0;n<n_roots;n++)
        if (std::fabs(root-roots[n])<=root_tol)
            return;

    roots[n_roots]=root;
    n_roots++;
}

static inline void sort_roots(double *roots,const mwSize n_roots)
{
    for (mwSize n=1;n<n_roots;n++)
    {
        const double root=roots[n];
        mwSize k=n;

        while ((k>0)&&(roots[k-1]>root))
        {
            roots[k]=roots[k-1];
            k--;
        }

        roots[k]=root;
    }
}

static mwSize solve_linear(const double c,const double d,double *roots,const double root_tol)
{
    mwSize n_roots=0;

    if (std::fabs(c)>root_tol)
        append_root(roots,n_roots,-d/c,root_tol);

    return n_roots;
}

static mwSize solve_quadratic(const double b2,const double c2,const double d2,double *roots,const double root_tol)
{
    if (std::fabs(b2)<=root_tol)
        return solve_linear(c2,d2,roots,root_tol);

    mwSize n_roots=0;
    double discr=c2*c2-4.0*b2*d2;
    const double discr_tol=64.0*std::numeric_limits<double>::epsilon()*(1.0+c2*c2+std::fabs(4.0*b2*d2));

    if (discr<-discr_tol)
        return 0;

    if (discr<0.0)
        discr=0.0;

    const double sqrt_discr=std::sqrt(discr);

    if (sqrt_discr==0.0)
    {
        append_root(roots,n_roots,-0.5*c2/b2,root_tol);
    }
    else
    {
        const double q=-0.5*(c2+std::copysign(sqrt_discr,c2));
        append_root(roots,n_roots,q/b2,root_tol);

        if (q!=0.0)
            append_root(roots,n_roots,d2/q,root_tol);
    }

    sort_roots(roots,n_roots);
    return n_roots;
}

static mwSize solve_cubic(const double a3,const double b3,const double c3,const double d3,double *roots,const double root_tol)
{
    const double coeff_scale=std::fmax(std::fmax(std::fabs(a3),std::fabs(b3)),
                                      std::fmax(std::fabs(c3),std::fabs(d3)));

    if (coeff_scale==0.0)
        return 0;

    const double a=a3/coeff_scale;
    const double b=b3/coeff_scale;
    const double c=c3/coeff_scale;
    const double d=d3/coeff_scale;

    if (std::fabs(a)<=root_tol)
        return solve_quadratic(b,c,d,roots,root_tol);

    mwSize n_roots=0;
    const double aa=b/a;
    const double bb=c/a;
    const double cc=d/a;
    const double aa2=aa*aa;
    const double p=bb-aa2/3.0;
    const double q=(2.0*aa*aa2)/27.0-(aa*bb)/3.0+cc;
    const double q_half=0.5*q;
    const double p_third=p/3.0;
    const double discr=q_half*q_half+p_third*p_third*p_third;
    const double shift=aa/3.0;

    if (discr>0.0)
    {
        const double sqrt_discr=std::sqrt(discr);
        const double u=std::cbrt(-q_half+sqrt_discr);
        const double v=std::cbrt(-q_half-sqrt_discr);
        append_root(roots,n_roots,u+v-shift,root_tol);
    }
    else if (discr==0.0)
    {
        const double u=std::cbrt(-q_half);
        append_root(roots,n_roots,2.0*u-shift,root_tol);
        append_root(roots,n_roots,-u-shift,root_tol);
    }
    else
    {
        if (p>=0.0)
        {
            const double sqrt_discr=std::sqrt(std::fmax(discr,0.0));
            const double u=std::cbrt(-q_half+sqrt_discr);
            const double v=std::cbrt(-q_half-sqrt_discr);
            append_root(roots,n_roots,u+v-shift,root_tol);
        }
        else
        {
            const double radius=2.0*std::sqrt(-p_third);
            const double arg=clamp_value(-q_half/std::sqrt(-p_third*p_third*p_third),-1.0,1.0);
            const double angle=std::acos(arg)/3.0;

            append_root(roots,n_roots,radius*std::cos(angle)-shift,root_tol);
            append_root(roots,n_roots,radius*std::cos(angle+2.0*pi/3.0)-shift,root_tol);
            append_root(roots,n_roots,radius*std::cos(angle-2.0*pi/3.0)-shift,root_tol);
        }
    }

    // Add repeated roots detected through stationary points
    double turn_roots[3]={0.0,0.0,0.0};
    const mwSize n_turns=solve_quadratic(3.0*a,2.0*b,c,turn_roots,root_tol);

    const double repeat_tol=std::fmax(128.0*std::numeric_limits<double>::epsilon(),root_tol*root_tol);

    for (mwSize n=0;n<n_turns;n++)
    {
        const double x=turn_roots[n];
        const double y=((a*x+b)*x+c)*x+d;

        if (std::fabs(y)<=repeat_tol)
            append_root(roots,n_roots,x,root_tol);
    }

    sort_roots(roots,n_roots);
    return n_roots;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

    // Get input arguments
    const double *poly_coeffs=mxGetDoubles(prhs[0]);
    const double root_tol=mxGetScalar(prhs[1]);

    // Solve the cubic equation
    double roots[3]={0.0,0.0,0.0};
    const mwSize n_roots=solve_cubic(poly_coeffs[0],poly_coeffs[1],
                                     poly_coeffs[2],poly_coeffs[3],
                                     roots,root_tol);

    // Allocate output vector
    mxArray *root_out=mxCreateDoubleMatrix(1,n_roots,mxREAL);

    if (root_out==nullptr)
        mexErrMsgIdAndTxt("Spinach:cubic_roots:alloc","Failed to allocate the output array.");

    // Copy roots into the output vector
    double *root_data=mxGetDoubles(root_out);

    if ((n_roots>0)&&(root_data==nullptr))
    {
        mxDestroyArray(root_out);
        mexErrMsgIdAndTxt("Spinach:cubic_roots:dataOut","Failed to access the output array data.");
    }

    for (mwSize n=0;n<n_roots;n++)
        root_data[n]=roots[n];

    // Return output argument
    plhs[0]=root_out;
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=2)
        mexErrMsgIdAndTxt("Spinach:cubic_roots:nrhs","Two inputs are required: poly_coeffs, root_tol.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:cubic_roots:nlhs","Exactly one output is required: root_list=cubic_roots(poly_coeffs,root_tol).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:cubic_roots:prhsNull","Input argument array pointer is null.");

    const mxArray *poly_coeffs=prhs[0];
    const mxArray *root_tol=prhs[1];

    if ((poly_coeffs==nullptr)||(root_tol==nullptr))
        mexErrMsgIdAndTxt("Spinach:cubic_roots:argNull","Input argument pointer is null.");

    if ((!mxIsDouble(poly_coeffs))||mxIsComplex(poly_coeffs)||mxIsSparse(poly_coeffs)||mxIsLogical(poly_coeffs)||
        (mxGetNumberOfElements(poly_coeffs)!=4))
        mexErrMsgIdAndTxt("Spinach:cubic_roots:coeffType","poly_coeffs must be a full real double vector with four elements.");

    if ((!mxIsDouble(root_tol))||mxIsComplex(root_tol)||mxIsSparse(root_tol)||mxIsLogical(root_tol)||
        (mxGetNumberOfElements(root_tol)!=1))
        mexErrMsgIdAndTxt("Spinach:cubic_roots:tolType","root_tol must be a full real double scalar.");

    const double *coeff_data=mxGetDoubles(poly_coeffs);

    if (coeff_data==nullptr)
        mexErrMsgIdAndTxt("Spinach:cubic_roots:dataIn","Failed to access coefficient data.");

    for (mwIndex n=0;n<4;n++)
        if (!std::isfinite(coeff_data[n]))
            mexErrMsgIdAndTxt("Spinach:cubic_roots:coeffFinite","poly_coeffs must contain finite real numbers.");

    const double tol=mxGetScalar(root_tol);

    if ((!std::isfinite(tol))||(tol<=0.0))
        mexErrMsgIdAndTxt("Spinach:cubic_roots:tolValue","root_tol must be a finite positive real scalar.");
}
