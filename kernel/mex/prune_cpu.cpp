/* prune_cpu.cpp
 *
 * Memory-friendly sparse clean-up for Spinach.
 *
 *   Aout = prune_cpu(A, nonzero_tol)
 *
 * Performs: Aout = tol*round(A/tol) and removes stored zeros.
 *
 * Design goal:
 *   Avoid the multiple full sparse temporaries produced by the generic MATLAB
 *   expression tol*round((1/tol)*A). This MEX allocates ONLY the output sparse
 *   matrix (plus O(ncols) bookkeeping) and streams over the input.
 *
 * Note on "in-place":
 *   Modern MATLAB releases actively block older undocumented APIs that were once
 *   used to detach shared arrays for in-place edits (e.g., mxUnshareArray).
 *   Modifying prhs[0] directly is unsupported and can crash MATLAB.
 *   Therefore this implementation is "single-extra-copy" (output only), not
 *   true in-place within the caller's sparse buffers.
 */

#include "mex.h"
#include <cmath>

static void grumble(int nlhs, int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("Spinach:prune_cpu:nrhs", "Two inputs required: A, nonzero_tol.");
    if (nlhs > 1)  mexErrMsgIdAndTxt("Spinach:prune_cpu:nlhs", "One output at most.");

    const mxArray *A = prhs[0];
    const mxArray *tol = prhs[1];

    if (!mxIsSparse(A))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:notSparse", "A must be sparse.");
    if (!mxIsDouble(A) || mxIsLogical(A))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:type", "A must be a sparse double array.");

    if (!mxIsDouble(tol) || mxIsComplex(tol) || mxGetNumberOfElements(tol) != 1)
        mexErrMsgIdAndTxt("Spinach:prune_cpu:tolType", "nonzero_tol must be a real double scalar.");

    const double nonzero_tol = mxGetScalar(tol);
    if (!(nonzero_tol > 0.0) || mxIsNaN(nonzero_tol) || mxIsInf(nonzero_tol))
        mexErrMsgIdAndTxt("Spinach:prune_cpu:tolVal", "nonzero_tol must be a finite positive real scalar.");
}

static inline double q(const double x, const double invtol, const double tol)
{
    return tol * std::round(x * invtol);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    grumble(nlhs, nrhs, prhs);

    const mxArray *A = prhs[0];
    const double tol = mxGetScalar(prhs[1]);
    const double invtol = 1.0 / tol;

    const mwSize m = mxGetM(A);
    const mwSize n = mxGetN(A);

    const mwIndex *jc = mxGetJc(A);
    const mwIndex *ir = mxGetIr(A);

    const double *pr = mxGetPr(A);
    const double *pi = mxGetPi(A); /* NULL for real */
    const bool is_cplx = (pi != NULL);

    /* First pass: count kept entries per column */
    mwIndex *colCounts = (mwIndex*) mxCalloc(n, sizeof(mwIndex));

    for (mwIndex col = 0; col < (mwIndex)n; col++) {
        const mwIndex start = jc[col];
        const mwIndex end   = jc[col + 1];
        mwIndex cnt = 0;
        for (mwIndex k = start; k < end; k++) {
            const double r = q(pr[k], invtol, tol);
            const double i = is_cplx ? q(pi[k], invtol, tol) : 0.0;
            if ((r != 0.0) || (i != 0.0)) cnt++;
        }
        colCounts[col] = cnt;
    }

    /* Prefix sum into jc_out */
    mwIndex *jc_out_tmp = (mwIndex*) mxCalloc(n + 1, sizeof(mwIndex));
    jc_out_tmp[0] = 0;
    for (mwIndex col = 0; col < (mwIndex)n; col++) {
        jc_out_tmp[col + 1] = jc_out_tmp[col] + colCounts[col];
    }
    const mwIndex nnz_out = jc_out_tmp[n];

    /* Allocate output */
    mxArray *Aout = mxCreateSparse(m, n, nnz_out, is_cplx ? mxCOMPLEX : mxREAL);

    mwIndex *jc_out = mxGetJc(Aout);
    mwIndex *ir_out = mxGetIr(Aout);
    double  *pr_out = mxGetPr(Aout);
    double  *pi_out = mxGetPi(Aout);

    /* Copy column pointers */
    for (mwIndex col = 0; col < (mwIndex)(n + 1); col++) jc_out[col] = jc_out_tmp[col];

    /* Second pass: fill */
    for (mwIndex col = 0; col < (mwIndex)n; col++) {
        mwIndex w = jc_out[col];
        const mwIndex start = jc[col];
        const mwIndex end   = jc[col + 1];

        for (mwIndex k = start; k < end; k++) {
            const double r = q(pr[k], invtol, tol);
            const double i = is_cplx ? q(pi[k], invtol, tol) : 0.0;
            if ((r != 0.0) || (i != 0.0)) {
                ir_out[w] = ir[k];
                pr_out[w] = r;
                if (is_cplx) pi_out[w] = i;
                w++;
            }
        }

        /* Sanity check */
        if (w != jc_out[col + 1])
            mexErrMsgIdAndTxt("Spinach:prune_cpu:internal", "Internal error: column write count mismatch.");
    }

    mxFree(jc_out_tmp);
    mxFree(colCounts);

    if (nlhs == 1) plhs[0] = Aout;
}
