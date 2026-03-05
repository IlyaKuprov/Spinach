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
 *   Modern MATLAB releases (including our minimum supported R2024b) do not
 *   provide a supported way to detach and modify prhs[0] sparse storage in-place.
 *   This implementation is therefore "single-extra-copy" (output only), not
 *   true in-place within the caller's sparse buffers.
 */

#include "mex.h"
#include <cmath>

static void grumble(int nlhs, int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("Spinach:prune_cpu:nrhs", "Two inputs required: A, nonzero_tol.");
    if (nlhs != 1) mexErrMsgIdAndTxt("Spinach:prune_cpu:nlhs", "Exactly one output is required: Aout = prune_cpu(A, nonzero_tol).");

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

    const bool is_cplx = mxIsComplex(A);

    /* First pass: count kept entries per column */
    mwIndex *colCounts = (mwIndex*) mxCalloc(n, sizeof(mwIndex));

    if (is_cplx) {
        const mxComplexDouble *px = mxGetComplexDoubles(A);
        for (mwIndex col = 0; col < (mwIndex)n; col++) {
            const mwIndex start = jc[col];
            const mwIndex end   = jc[col + 1];
            mwIndex cnt = 0;
            for (mwIndex k = start; k < end; k++) {
                const double r = q(px[k].real, invtol, tol);
                const double i = q(px[k].imag, invtol, tol);
                if ((r != 0.0) || (i != 0.0)) cnt++;
            }
            colCounts[col] = cnt;
        }
    } else {
        const double *pr = mxGetDoubles(A);
        for (mwIndex col = 0; col < (mwIndex)n; col++) {
            const mwIndex start = jc[col];
            const mwIndex end   = jc[col + 1];
            mwIndex cnt = 0;
            for (mwIndex k = start; k < end; k++) {
                const double r = q(pr[k], invtol, tol);
                if (r != 0.0) cnt++;
            }
            colCounts[col] = cnt;
        }
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

    /* Copy column pointers */
    for (mwIndex col = 0; col < (mwIndex)(n + 1); col++) jc_out[col] = jc_out_tmp[col];

    /* Second pass: fill */
    if (is_cplx) {
        const mxComplexDouble *px_in = mxGetComplexDoubles(A);
        mxComplexDouble *px_out = mxGetComplexDoubles(Aout);

        for (mwIndex col = 0; col < (mwIndex)n; col++) {
            mwIndex w = jc_out[col];
            const mwIndex start = jc[col];
            const mwIndex end   = jc[col + 1];

            for (mwIndex k = start; k < end; k++) {
                const double r = q(px_in[k].real, invtol, tol);
                const double i = q(px_in[k].imag, invtol, tol);
                if ((r != 0.0) || (i != 0.0)) {
                    ir_out[w] = ir[k];
                    px_out[w].real = r;
                    px_out[w].imag = i;
                    w++;
                }
            }

            /* Sanity check */
            if (w != jc_out[col + 1])
                mexErrMsgIdAndTxt("Spinach:prune_cpu:internal", "Internal error: column write count mismatch.");
        }
    } else {
        const double *pr_in = mxGetDoubles(A);
        double *pr_out = mxGetDoubles(Aout);

        for (mwIndex col = 0; col < (mwIndex)n; col++) {
            mwIndex w = jc_out[col];
            const mwIndex start = jc[col];
            const mwIndex end   = jc[col + 1];

            for (mwIndex k = start; k < end; k++) {
                const double r = q(pr_in[k], invtol, tol);
                if (r != 0.0) {
                    ir_out[w] = ir[k];
                    pr_out[w] = r;
                    w++;
                }
            }

            /* Sanity check */
            if (w != jc_out[col + 1])
                mexErrMsgIdAndTxt("Spinach:prune_cpu:internal", "Internal error: column write count mismatch.");
        }
    }

    mxFree(jc_out_tmp);
    mxFree(colCounts);

    if (nlhs == 1) plhs[0] = Aout;
}
