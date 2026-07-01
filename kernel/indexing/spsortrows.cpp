/* spsortrows.cpp
 *
 * Sparse real double sortrows index MEX utility for Spinach
 *
 * Syntax:
 *
 *                    idx=spsortrows(A)
 *
 * Performs:
 *
 *                    [~,idx]=sortrows(A)
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

#if defined(_OPENMP)
#include <omp.h>
#endif

struct sparse_row_less
{
    const mwIndex *row_ptr;
    const mwIndex *col_idx;
    const double *values;

    bool operator()(const mwIndex row_a,const mwIndex row_b) const;
};

static inline int compare_values(const double val_a,const double val_b);
static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

#ifdef _OPENMP

    // Do what Matlab says
    omp_set_dynamic(0);

#endif

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
        // Count numerically non-zero entries in every row
        std::vector<mwIndex> row_ptr(n_rows+1,0);
        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            for (mwIndex ptr=jc[col];ptr<jc[col+1];ptr++)
            {
                const mwIndex row=ir[ptr];

                if (row>=(mwIndex)n_rows)
                    mexErrMsgIdAndTxt("Spinach:spsortrows:badIr","Sparse row index exceeds matrix dimensions.");

                if (pr[ptr]!=0.0)
                    row_ptr[row+1]++;
            }
        }

        // Prefix-sum row counts into CSR row pointers
        for (mwIndex row=0;row<(mwIndex)n_rows;row++)
            row_ptr[row+1]+=row_ptr[row];

        // Allocate CSR column and value arrays
        const mwIndex nnz=row_ptr[n_rows];
        std::vector<mwIndex> row_next(row_ptr);
        std::vector<mwIndex> col_idx(nnz);
        std::vector<double> values(nnz);

        // Fill CSR arrays in column order
        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            for (mwIndex ptr=jc[col];ptr<jc[col+1];ptr++)
            {
                if (pr[ptr]!=0.0)
                {
                    const mwIndex row=ir[ptr];
                    const mwIndex out_ptr=row_next[row]++;
                    col_idx[out_ptr]=col;
                    values[out_ptr]=pr[ptr];
                }
            }
        }

        // Initialise zero-based row permutation
        std::vector<mwIndex> order(n_rows);
        std::iota(order.begin(),order.end(),(mwIndex)0);

        // Sort rows lexicographically using sparse row descriptors
        sparse_row_less row_less;
        row_less.row_ptr=row_ptr.data();
        row_less.col_idx=col_idx.data();
        row_less.values=values.data();

#if defined(_OPENMP)

        // Get OpenMP team size
        const int max_threads=omp_get_max_threads();

        if ((max_threads>1)&&(n_rows>1))
        {
            // Cap the team at one thread per row
            int n_threads=max_threads;

            if ((mwSize)n_threads>n_rows)
                n_threads=(int)n_rows;

            // Split the row permutation into balanced blocks
            std::vector<mwIndex> run_start((mwIndex)n_threads+1);
            const mwIndex chunk=(mwIndex)(n_rows/(mwSize)n_threads);
            const mwIndex rem=(mwIndex)(n_rows%(mwSize)n_threads);
            mwIndex run_pos=0;

            for (int n=0;n<n_threads;n++)
            {
                run_start[(mwIndex)n]=run_pos;
                run_pos+=chunk+((n<(int)rem)?1:0);
            }

            run_start[(mwIndex)n_threads]=run_pos;

            // Sort local blocks in parallel
#pragma omp parallel for schedule(static) num_threads(n_threads)
            for (int n=0;n<n_threads;n++)
                std::sort(order.begin()+run_start[(mwIndex)n],
                          order.begin()+run_start[(mwIndex)n+1],row_less);

            // Allocate merge workspace
            std::vector<mwIndex> buffer(n_rows);
            std::vector<mwIndex> *src=&order;
            std::vector<mwIndex> *dst=&buffer;

            // Merge sorted blocks using a reduced GNU-style sort/merge schedule
            for (int run_len=1;run_len<n_threads;run_len*=2)
            {
                const int pair_count=(n_threads+2*run_len-1)/(2*run_len);

#pragma omp parallel for schedule(static) num_threads(n_threads)
                for (int pair=0;pair<pair_count;pair++)
                {
                    const int left_run=pair*2*run_len;
                    const int mid_run=std::min(left_run+run_len,n_threads);
                    const int right_run=std::min(left_run+2*run_len,n_threads);
                    const mwIndex left=run_start[(mwIndex)left_run];
                    const mwIndex mid=run_start[(mwIndex)mid_run];
                    const mwIndex right=run_start[(mwIndex)right_run];

                    if (mid<right)
                        std::merge(src->begin()+left,src->begin()+mid,
                                   src->begin()+mid,src->begin()+right,
                                   dst->begin()+left,row_less);
                    else
                        std::copy(src->begin()+left,src->begin()+right,
                                  dst->begin()+left);
                }

                std::swap(src,dst);
            }

            // Copy the final pass back if it ended in the workspace
            if (src!=&order)
                std::copy(src->begin(),src->end(),order.begin());
        }
        else
        {
            // Use standard introsort when there is no useful OpenMP team
            std::sort(order.begin(),order.end(),row_less);
        }

#else

        // Use standard introsort otherwise
        std::sort(order.begin(),order.end(),row_less);

#endif

        // Allocate output permutation vector
        plhs[0]=mxCreateDoubleMatrix(n_rows,1,mxREAL);

        if (plhs[0]==nullptr)
            mexErrMsgIdAndTxt("Spinach:spsortrows:allocOut","Failed to allocate output index vector.");

        // Fill output permutation vector
        double *idx=mxGetDoubles(plhs[0]);

        if ((n_rows!=0)&&(idx==nullptr))
            mexErrMsgIdAndTxt("Spinach:spsortrows:dataOut","Failed to access output index data.");

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (mwSignedIndex n=0;n<(mwSignedIndex)n_rows;n++)
            idx[(mwIndex)n]=(double)(order[(mwIndex)n]+1);
    }
    catch (const std::bad_alloc&)
    {
        mexErrMsgIdAndTxt("Spinach:spsortrows:badAlloc","Memory allocation failed.");
    }
    catch (const std::exception& err)
    {
        mexErrMsgIdAndTxt("Spinach:spsortrows:cppException","C++ exception: %s",err.what());
    }
}

bool sparse_row_less::operator()(const mwIndex row_a,const mwIndex row_b) const
{
    mwIndex ptr_a=row_ptr[row_a];
    mwIndex ptr_b=row_ptr[row_b];
    const mwIndex end_a=row_ptr[row_a+1];
    const mwIndex end_b=row_ptr[row_b+1];
    const mwIndex no_col=std::numeric_limits<mwIndex>::max();

    while ((ptr_a<end_a)||(ptr_b<end_b))
    {
        const mwIndex col_a=(ptr_a<end_a)?col_idx[ptr_a]:no_col;
        const mwIndex col_b=(ptr_b<end_b)?col_idx[ptr_b]:no_col;
        double val_a=0.0;
        double val_b=0.0;

        if (col_a==col_b)
        {
            val_a=values[ptr_a++];
            val_b=values[ptr_b++];
        }
        else if (col_a<col_b)
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

    return row_a<row_b;
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

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=1)
        mexErrMsgIdAndTxt("Spinach:spsortrows:nrhs","One input is required: A.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:spsortrows:nlhs","Exactly one output is required: idx=spsortrows(A).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:spsortrows:prhsNull","Input argument array pointer is null.");

    const mxArray *A=prhs[0];

    if (A==nullptr)
        mexErrMsgIdAndTxt("Spinach:spsortrows:argNull","Input argument A is null.");

    if (!mxIsNumeric(A))
        mexErrMsgIdAndTxt("Spinach:spsortrows:notNumeric","A must be numeric.");

    if (!mxIsSparse(A))
        mexErrMsgIdAndTxt("Spinach:spsortrows:notSparse","A must be sparse.");

    if (mxIsComplex(A))
        mexErrMsgIdAndTxt("Spinach:spsortrows:complex","A must be real.");

    if (mxIsLogical(A))
        mexErrMsgIdAndTxt("Spinach:spsortrows:logical","A must not be logical.");

    if (mxGetClassID(A)!=mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("Spinach:spsortrows:notDouble","A must be double precision.");

    if (mxGetNumberOfDimensions(A)!=2)
        mexErrMsgIdAndTxt("Spinach:spsortrows:not2d","A must be a two-dimensional matrix.");

    const mwSize n_rows=mxGetM(A);
    const mwSize n_cols=mxGetN(A);
    const mwSize row_lim=(mwSize)std::numeric_limits<mwSignedIndex>::max();

    if (n_rows>row_lim)
        mexErrMsgIdAndTxt("Spinach:spsortrows:rowLim","A has too many rows for index conversion.");

    const mwSize max_elem=std::numeric_limits<mwSize>::max();

    if ((n_cols!=0)&&(n_rows>(max_elem/n_cols)))
        mexErrMsgIdAndTxt("Spinach:spsortrows:elemLim","A dimensions overflow element indexing.");

    if (mxGetJc(A)==nullptr)
        mexErrMsgIdAndTxt("Spinach:spsortrows:jcAccess","A column pointer is null.");

    if ((mxGetNzmax(A)!=0)&&(mxGetIr(A)==nullptr))
        mexErrMsgIdAndTxt("Spinach:spsortrows:irAccess","A row index pointer is null.");

    if ((mxGetNzmax(A)!=0)&&(mxGetDoubles(A)==nullptr))
        mexErrMsgIdAndTxt("Spinach:spsortrows:dataAccess","A data pointer is null.");
}
