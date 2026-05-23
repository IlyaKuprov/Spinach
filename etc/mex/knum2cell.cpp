/* knum2cell.cpp
 *
 * Row-partitioning MEX utility for Spinach
 *
 * Syntax:
 *
 *                    c=knum2cell(a)
 *
 * Performs:
 *
 *                    c=num2cell(a,2)
 *
 * for real single-precision matrices with at least two rows
 * and at least two columns
 */

#include "mex.h"
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

static void release_rows(mxArray **row_cells,mwSize n_rows);
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
    const mxArray *a=prhs[0];
    const mwSize n_rows=mxGetM(a);
    const mwSize n_cols=mxGetN(a);

    // Get input matrix data pointer
    const float *a_data=mxGetSingles(a);

    // Halt if matrix data access failed
    if (a_data==nullptr)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:dataPtr","Failed to access input matrix data.");

    // Allocate output cell array
    mxArray *c_out=mxCreateCellMatrix(n_rows,1);

    // Halt if output cell allocation failed
    if (c_out==nullptr)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:allocCell","Failed to allocate output cell array.");

    // Allocate temporary row descriptor arrays
    mxArray **row_cells=(mxArray**)mxCalloc(n_rows,sizeof(mxArray*));
    float **row_data=(float**)mxCalloc(n_rows,sizeof(float*));

    // Halt if temporary descriptor allocation failed
    if ((row_cells==nullptr)||(row_data==nullptr))
    {
        mxFree(row_data);
        mxFree(row_cells);
        mxDestroyArray(c_out);
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:allocTemp","Failed to allocate temporary descriptors.");
    }

    // Allocate output row vectors and resolve their data pointers
    for (mwSize row=0;row<n_rows;row++)
    {
        row_cells[row]=mxCreateUninitNumericMatrix(1,n_cols,mxSINGLE_CLASS,mxREAL);

        if (row_cells[row]==nullptr)
        {
            release_rows(row_cells,row);
            mxFree(row_data);
            mxFree(row_cells);
            mxDestroyArray(c_out);
            mexErrMsgIdAndTxt("Spinach:num2cell_2_single:allocRow","Failed to allocate an output row vector.");
        }

        row_data[row]=mxGetSingles(row_cells[row]);

        if (row_data[row]==nullptr)
        {
            release_rows(row_cells,row+1);
            mxFree(row_data);
            mxFree(row_cells);
            mxDestroyArray(c_out);
            mexErrMsgIdAndTxt("Spinach:num2cell_2_single:dataRow","Failed to access output row vector data.");
        }
    }

    // Copy matrix rows into output row vectors
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (mwSignedIndex row=0;row<(mwSignedIndex)n_rows;row++)
    {
        const mwSize row_idx=(mwSize)row;
        float *dst=row_data[row_idx];

        for (mwSize col=0;col<n_cols;col++)
            dst[col]=a_data[row_idx+col*n_rows];
    }

    // Assign row vectors into the output cell array
    for (mwSize row=0;row<n_rows;row++)
        mxSetCell(c_out,row,row_cells[row]);

    // Release temporary row descriptor arrays
    mxFree(row_data);
    mxFree(row_cells);

    // Return output argument
    plhs[0]=c_out;
}

static void release_rows(mxArray **row_cells,mwSize n_rows)
{
    for (mwSize row=0;row<n_rows;row++)
        if (row_cells[row]!=nullptr)
            mxDestroyArray(row_cells[row]);
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=1)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:nrhs","One input is required: A.");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:nlhs","Exactly one output is required: C=num2cell_2_single(A).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:prhsNull","Input argument array pointer is null.");

    const mxArray *a=prhs[0];

    if (a==nullptr)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:argNull","Input argument A is null.");

    if (!mxIsNumeric(a))
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:notNumeric","A must be numeric.");

    if (mxIsSparse(a))
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:sparse","A must be full.");

    if (mxIsComplex(a))
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:complex","A must be real.");

    if (mxIsLogical(a))
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:logical","A must not be logical.");

    if (mxGetClassID(a)!=mxSINGLE_CLASS)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:notSingle","A must be single precision.");

    if (mxGetNumberOfDimensions(a)!=2)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:not2d","A must be a two-dimensional matrix.");

    const mwSize n_rows=mxGetM(a);
    const mwSize n_cols=mxGetN(a);

    if (n_rows<2)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:rows","A must have more than one row.");

    if (n_cols<2)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:cols","A must have more than one column.");

    const mwSize row_lim=(mwSize)std::numeric_limits<mwSignedIndex>::max();

    if (n_rows>row_lim)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:rowLim","A has too many rows for index conversion.");

    const mwSize max_elem=std::numeric_limits<mwSize>::max();

    if ((n_cols!=0)&&(n_rows>(max_elem/n_cols)))
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:elemLim","A dimensions overflow element indexing.");

    if (mxGetSingles(a)==nullptr)
        mexErrMsgIdAndTxt("Spinach:num2cell_2_single:dataAccess","A data pointer is null.");
}
