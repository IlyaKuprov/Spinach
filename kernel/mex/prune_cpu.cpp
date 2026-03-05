/* prune_cpu.cpp
 *
 * Sparse clean-up MEX utility for Spinach
 *
 * Syntax:
 *
 *                    Aout=prune_cpu(A,nonzero_tol)
 *
 * Performs:
 *
 *                    Aout=nonzero_tol*round(A/nonzero_tol)
 *
 * and removes stored zero entries
 */

#include "mex.h"
#include <cmath>

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=2)
        mexErrMsgIdAndTxt("Spinach:prune_cpu:nrhs","Two inputs required: A, nonzero_tol.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:prune_cpu:nlhs","Exactly one output is required: Aout=prune_cpu(A,nonzero_tol).");

    const mxArray *a=prhs[0];
    const mxArray *nonzero_tol=prhs[1];

    if (!mxIsSparse(a))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:notSparse","A must be sparse.");

    if (!mxIsDouble(a)||mxIsLogical(a))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:type","A must be a sparse double array.");

    if (!mxIsDouble(nonzero_tol)||mxIsComplex(nonzero_tol)||(mxGetNumberOfElements(nonzero_tol)!=1))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:tolType","nonzero_tol must be a real double scalar.");

    const double tol=mxGetScalar(nonzero_tol);

    if (!(tol>0.0)||mxIsNaN(tol)||mxIsInf(tol))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:tolVal","nonzero_tol must be a finite positive real scalar.");
}

static inline double quantise_value(const double x,const double inv_tol,const double tol)
{
    return tol*std::round(x*inv_tol);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

    // Get input arguments
    const mxArray *a=prhs[0];
    const double tol=mxGetScalar(prhs[1]);
    const double inv_tol=1.0/tol;

    // Get matrix dimensions
    const mwSize n_rows=mxGetM(a);
    const mwSize n_cols=mxGetN(a);

    // Get sparse matrix structural arrays
    const mwIndex *jc=mxGetJc(a);
    const mwIndex *ir=mxGetIr(a);

    // Determine matrix complexity
    const bool is_cplx=mxIsComplex(a);

    // Allocate per-column entry counters
    mwIndex *col_counts=(mwIndex*)mxCalloc(n_cols,sizeof(mwIndex));

    // Count retained entries in each column
    if (is_cplx)
    {
        const mxComplexDouble *px=mxGetComplexDoubles(a);

        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            const mwIndex start=jc[col];
            const mwIndex finish=jc[col+1];
            mwIndex count=0;

            for (mwIndex k=start;k<finish;k++)
            {
                const double real_part=quantise_value(px[k].real,inv_tol,tol);
                const double imag_part=quantise_value(px[k].imag,inv_tol,tol);

                if ((real_part!=0.0)||(imag_part!=0.0))
                    count++;
            }

            col_counts[col]=count;
        }
    }
    else
    {
        const double *pr=mxGetDoubles(a);

        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            const mwIndex start=jc[col];
            const mwIndex finish=jc[col+1];
            mwIndex count=0;

            for (mwIndex k=start;k<finish;k++)
            {
                const double real_part=quantise_value(pr[k],inv_tol,tol);

                if (real_part!=0.0)
                    count++;
            }

            col_counts[col]=count;
        }
    }

    // Allocate temporary output column pointers
    mwIndex *jc_out_tmp=(mwIndex*)mxCalloc(n_cols+1,sizeof(mwIndex));

    // Build output column pointers by prefix sum
    jc_out_tmp[0]=0;

    for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        jc_out_tmp[col+1]=jc_out_tmp[col]+col_counts[col];

    const mwIndex nnz_out=jc_out_tmp[n_cols];

    // Allocate output sparse matrix
    mxArray *a_out=mxCreateSparse(n_rows,n_cols,nnz_out,is_cplx?mxCOMPLEX:mxREAL);

    // Get output sparse matrix data pointers
    mwIndex *jc_out=mxGetJc(a_out);
    mwIndex *ir_out=mxGetIr(a_out);

    // Copy output column pointer array
    for (mwIndex col=0;col<(mwIndex)(n_cols+1);col++)
        jc_out[col]=jc_out_tmp[col];

    // Fill output sparse matrix with retained entries
    if (is_cplx)
    {
        const mxComplexDouble *px_in=mxGetComplexDoubles(a);
        mxComplexDouble *px_out=mxGetComplexDoubles(a_out);

        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            mwIndex write_ptr=jc_out[col];
            const mwIndex start=jc[col];
            const mwIndex finish=jc[col+1];

            for (mwIndex k=start;k<finish;k++)
            {
                const double real_part=quantise_value(px_in[k].real,inv_tol,tol);
                const double imag_part=quantise_value(px_in[k].imag,inv_tol,tol);

                if ((real_part!=0.0)||(imag_part!=0.0))
                {
                    ir_out[write_ptr]=ir[k];
                    px_out[write_ptr].real=real_part;
                    px_out[write_ptr].imag=imag_part;
                    write_ptr++;
                }
            }

            if (write_ptr!=jc_out[col+1])
                mexErrMsgIdAndTxt("Spinach:prune_cpu:internal","Internal error: column write count mismatch.");
        }
    }
    else
    {
        const double *pr_in=mxGetDoubles(a);
        double *pr_out=mxGetDoubles(a_out);

        for (mwIndex col=0;col<(mwIndex)n_cols;col++)
        {
            mwIndex write_ptr=jc_out[col];
            const mwIndex start=jc[col];
            const mwIndex finish=jc[col+1];

            for (mwIndex k=start;k<finish;k++)
            {
                const double real_part=quantise_value(pr_in[k],inv_tol,tol);

                if (real_part!=0.0)
                {
                    ir_out[write_ptr]=ir[k];
                    pr_out[write_ptr]=real_part;
                    write_ptr++;
                }
            }

            if (write_ptr!=jc_out[col+1])
                mexErrMsgIdAndTxt("Spinach:prune_cpu:internal","Internal error: column write count mismatch.");
        }
    }

    // Release temporary arrays
    mxFree(jc_out_tmp);
    mxFree(col_counts);

    // Return output argument
    plhs[0]=a_out;
}
