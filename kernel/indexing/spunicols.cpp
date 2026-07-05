/* spunicols.cpp
 *
 * Sparse real double unique-column MEX utility for Spinach
 *
 * Syntax:
 *
 *                    A=spunicols(A)
 *
 * Performs the equivalent of:
 *
 *                    A=unique(A.','rows').'
 *
 * for sparse real double matrices
 */

#include "mex.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <vector>

struct sparse_col_less
{
    const mwIndex *jc;
    const mwIndex *ir;
    const double *values;

    bool operator()(const mwIndex col_a,const mwIndex col_b) const;
};

static inline int compare_values(const double val_a,const double val_b);
static bool columns_equal(const mwIndex col_a,const mwIndex col_b,
                          const mwIndex *jc,const mwIndex *ir,
                          const double *values);
static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

    // Get input matrix descriptor and dimensions
    const mxArray *A=prhs[0];
    const mwSize n_rows=mxGetM(A);
    const mwSize n_cols=mxGetN(A);

    // Get sparse matrix structural arrays
    const mwIndex *jc=mxGetJc(A);
    const mwIndex *ir=mxGetIr(A);
    const double *pr=mxGetDoubles(A);

    try
    {
        // Check sparse structure before comparator use
        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
            if (jc[col]>jc[col+1])
                mexErrMsgIdAndTxt("Spinach:spunicols:badJc","Sparse column pointers are not monotonic.");

        if (jc[n_cols]>mxGetNzmax(A))
            mexErrMsgIdAndTxt("Spinach:spunicols:badJcEnd","Sparse column pointer exceeds allocated storage.");

        for (mwIndex ptr=0;ptr<jc[n_cols];ptr++)
            if (ir[ptr]>=(mwIndex)n_rows)
                mexErrMsgIdAndTxt("Spinach:spunicols:badIr","Sparse row index exceeds matrix dimensions.");

        // Initialise zero-based column permutation
        std::vector<mwIndex> order(n_cols);
        std::iota(order.begin(),order.end(),(mwIndex)0);

        // Sort columns lexicographically using sparse column descriptors
        sparse_col_less col_less;
        col_less.jc=jc;
        col_less.ir=ir;
        col_less.values=pr;
        std::sort(order.begin(),order.end(),col_less);

        // Keep the first entry in each equal-column run
        std::vector<mwIndex> keep;
        keep.reserve(n_cols);
        for (mwIndex n=0;n<(mwIndex)n_cols;n++)
        {
            if ((n==0)||(!columns_equal(order[n-1],order[n],jc,ir,pr)))
                keep.push_back(order[n]);
        }

        // Count numerically non-zero entries in kept columns
        mwSize nnz_out=0;
        for (mwIndex n=0;n<(mwIndex)keep.size();n++)
        {
            const mwIndex col=keep[n];
            for (mwIndex ptr=jc[col];ptr<jc[col+1];ptr++)
                if (pr[ptr]!=0.0)
                    nnz_out++;
        }

        // Allocate output sparse matrix
        plhs[0]=mxCreateSparse(n_rows,(mwSize)keep.size(),nnz_out,mxREAL);

        if (plhs[0]==nullptr)
            mexErrMsgIdAndTxt("Spinach:spunicols:allocOut","Failed to allocate output sparse matrix.");

        // Get output sparse matrix structural arrays
        mwIndex *jc_out=mxGetJc(plhs[0]);
        mwIndex *ir_out=mxGetIr(plhs[0]);
        double *pr_out=mxGetDoubles(plhs[0]);

        if (jc_out==nullptr)
            mexErrMsgIdAndTxt("Spinach:spunicols:jcOut","Failed to access output column pointers.");

        if ((nnz_out!=0)&&((ir_out==nullptr)||(pr_out==nullptr)))
            mexErrMsgIdAndTxt("Spinach:spunicols:dataOut","Failed to access output sparse data.");

        // Copy kept columns into the output matrix
        mwIndex out_ptr=0;
        jc_out[0]=0;
        for (mwIndex n=0;n<(mwIndex)keep.size();n++)
        {
            const mwIndex col=keep[n];
            for (mwIndex ptr=jc[col];ptr<jc[col+1];ptr++)
            {
                if (pr[ptr]!=0.0)
                {
                    ir_out[out_ptr]=ir[ptr];
                    pr_out[out_ptr]=pr[ptr];
                    out_ptr++;
                }
            }
            jc_out[n+1]=out_ptr;
        }
    }
    catch (const std::bad_alloc&)
    {
        mexErrMsgIdAndTxt("Spinach:spunicols:badAlloc","Memory allocation failed.");
    }
    catch (const std::exception& err)
    {
        mexErrMsgIdAndTxt("Spinach:spunicols:cppException","C++ exception: %s",err.what());
    }
}

bool sparse_col_less::operator()(const mwIndex col_a,const mwIndex col_b) const
{
    mwIndex ptr_a=jc[col_a];
    mwIndex ptr_b=jc[col_b];
    const mwIndex end_a=jc[col_a+1];
    const mwIndex end_b=jc[col_b+1];
    const mwIndex no_row=std::numeric_limits<mwIndex>::max();

    while ((ptr_a<end_a)||(ptr_b<end_b))
    {
        while ((ptr_a<end_a)&&(values[ptr_a]==0.0))
            ptr_a++;

        while ((ptr_b<end_b)&&(values[ptr_b]==0.0))
            ptr_b++;

        if ((ptr_a>=end_a)&&(ptr_b>=end_b))
            break;

        const mwIndex row_a=(ptr_a<end_a)?ir[ptr_a]:no_row;
        const mwIndex row_b=(ptr_b<end_b)?ir[ptr_b]:no_row;
        double val_a=0.0;
        double val_b=0.0;

        if (row_a==row_b)
        {
            val_a=values[ptr_a++];
            val_b=values[ptr_b++];
        }
        else if (row_a<row_b)
        {
            val_a=values[ptr_a++];
        }
        else
        {
            val_b=values[ptr_b++];
        }

        const int cmp=compare_values(val_a,val_b);

        if (cmp<0)
            return true;

        if (cmp>0)
            return false;
    }

    return col_a<col_b;
}

static inline int compare_values(const double val_a,const double val_b)
{
    const bool nan_a=std::isnan(val_a);
    const bool nan_b=std::isnan(val_b);

    if (nan_a||nan_b)
    {
        if (nan_a&&nan_b)
            return 0;

        return nan_a?1:-1;
    }

    if (val_a<val_b)
        return -1;

    if (val_a>val_b)
        return 1;

    return 0;
}

static bool columns_equal(const mwIndex col_a,const mwIndex col_b,
                          const mwIndex *jc,const mwIndex *ir,
                          const double *values)
{
    mwIndex ptr_a=jc[col_a];
    mwIndex ptr_b=jc[col_b];
    const mwIndex end_a=jc[col_a+1];
    const mwIndex end_b=jc[col_b+1];

    while ((ptr_a<end_a)||(ptr_b<end_b))
    {
        while ((ptr_a<end_a)&&(values[ptr_a]==0.0))
            ptr_a++;

        while ((ptr_b<end_b)&&(values[ptr_b]==0.0))
            ptr_b++;

        if ((ptr_a>=end_a)&&(ptr_b>=end_b))
            return true;

        if ((ptr_a>=end_a)||(ptr_b>=end_b))
            return false;

        if (ir[ptr_a]!=ir[ptr_b])
            return false;

        if (values[ptr_a]!=values[ptr_b])
            return false;

        ptr_a++;
        ptr_b++;
    }

    return true;
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=1)
        mexErrMsgIdAndTxt("Spinach:spunicols:nrhs","One input is required: A.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:spunicols:nlhs","Exactly one output is required: A=spunicols(A).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:spunicols:prhsNull","Input argument array pointer is null.");

    const mxArray *A=prhs[0];

    if (A==nullptr)
        mexErrMsgIdAndTxt("Spinach:spunicols:argNull","Input argument A is null.");

    if (!mxIsNumeric(A))
        mexErrMsgIdAndTxt("Spinach:spunicols:notNumeric","A must be numeric.");

    if (!mxIsSparse(A))
        mexErrMsgIdAndTxt("Spinach:spunicols:notSparse","A must be sparse.");

    if (mxIsComplex(A))
        mexErrMsgIdAndTxt("Spinach:spunicols:complex","A must be real.");

    if (mxIsLogical(A))
        mexErrMsgIdAndTxt("Spinach:spunicols:logical","A must not be logical.");

    if (mxGetClassID(A)!=mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("Spinach:spunicols:notDouble","A must be double precision.");

    if (mxGetNumberOfDimensions(A)!=2)
        mexErrMsgIdAndTxt("Spinach:spunicols:not2d","A must be a two-dimensional matrix.");

    const mwSize n_rows=mxGetM(A);
    const mwSize n_cols=mxGetN(A);
    const mwSize max_elem=std::numeric_limits<mwSize>::max();

    if ((n_cols!=0)&&(n_rows>(max_elem/n_cols)))
        mexErrMsgIdAndTxt("Spinach:spunicols:elemLim","A dimensions overflow element indexing.");

    if (mxGetJc(A)==nullptr)
        mexErrMsgIdAndTxt("Spinach:spunicols:jcAccess","A column pointer is null.");

    if ((mxGetNzmax(A)!=0)&&(mxGetIr(A)==nullptr))
        mexErrMsgIdAndTxt("Spinach:spunicols:irAccess","A row index pointer is null.");

    if ((mxGetNzmax(A)!=0)&&(mxGetDoubles(A)==nullptr))
        mexErrMsgIdAndTxt("Spinach:spunicols:dataAccess","A data pointer is null.");
}
