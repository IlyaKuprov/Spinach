/* fsparse.cpp
 *
 * Fast sparse real double assembly MEX utility for Spinach
 *
 * Syntax:
 *
 *                    A=fsparse(row_idx,col_idx,vals,n_rows,n_cols)
 *
 * The sparse assembly algorithm is derived from Stefan Engblom's
 * fsparse.c implementation in stenglib, commit
 * 4e0183b0b3e6150d4d810c44e7a7151b66410f56.
 *
 * Reference:
 *
 * S. Engblom and D. Lukarski, Fast Matlab compatible sparse assembly
 * on multicore computers, Parallel Computing 56, 1-17 (2016).
 *
 * stenglib licence statement:
 *
 * You may download all of stenglib and use, modify and redistribute
 * it in any way you like. A redistributor must fully attribute the
 * authorship and make a good effort to cite the original location of
 * the software. A researcher making critical use of the software in
 * research is requested to acknowledge this in publications related to
 * the research. A company may use the code in software products
 * provided that the original location and the author is clearly cited.
 */

#include "mex.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);
static int get_dim(const mxArray *array_ptr,const char *arg_name);
static int get_len(const mxArray *array_ptr);
static int *get_indices(const mxArray *array_ptr,int max_idx,const char *arg_name);
static mxArray *call_sparse(int nrhs,const mxArray *prhs[]);
static void squeeze(mxArray *S);
static void sparse_insert(mwIndex *irS,double *prS,const int *irank,
                          const int *rank,const mwIndex *jrS,
                          const int *ii,const double *sr,
                          int smod,int sdiv,int len,int M);
static mxArray *sparse_real(const int *ii,const int *jj,const double *sr,
                            int smod,int sdiv,int len,int M,int N);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // Validate input arguments
    grumble(nlhs,nrhs,prhs);

    // Fall back to Matlab sparse for complex values
    if (mxIsComplex(prhs[2]))
    {
        plhs[0]=call_sparse(nrhs,prhs);
        return;
    }

    try
    {
        // Get output matrix dimensions
        const int M=get_dim(prhs[3],"n_rows");
        const int N=get_dim(prhs[4],"n_cols");

        // Get triplet count and input values
        const int len=get_len(prhs[0]);
        const mwSize val_len=mxGetNumberOfElements(prhs[2]);
        const double *vals=mxGetDoubles(prhs[2]);

        // Convert indices to the internal integer representation
        int *row_idx=get_indices(prhs[0],M,"row_idx");
        int *col_idx=get_indices(prhs[1],N,"col_idx");

        // Select full or scalar value insertion
        const bool scalar_vals=((val_len==1)&&(len!=1));
        const int smod=scalar_vals?1:len;
        const int sdiv=scalar_vals?0:1;

        // Assemble the sparse matrix
        plhs[0]=sparse_real(row_idx,col_idx,vals,smod,sdiv,len,M,N);

        // Remove zeros produced by explicit zero values or cancellation
        squeeze(plhs[0]);

        // Release index work arrays
        mxFree(row_idx);
        mxFree(col_idx);
    }
    catch (const std::bad_alloc&)
    {
        mexErrMsgIdAndTxt("Spinach:fsparse:badAlloc","Memory allocation failed.");
    }
    catch (const std::exception& err)
    {
        mexErrMsgIdAndTxt("Spinach:fsparse:cppException","C++ exception: %s",err.what());
    }
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=5)
        mexErrMsgIdAndTxt("Spinach:fsparse:nrhs","Five inputs are required: A=fsparse(row_idx,col_idx,vals,n_rows,n_cols).");

    if (nlhs!=1)
        mexErrMsgIdAndTxt("Spinach:fsparse:nlhs","Exactly one output is required: A=fsparse(row_idx,col_idx,vals,n_rows,n_cols).");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:prhsNull","Input argument array pointer is null.");

    for (int n=0;n<5;n++)
        if (prhs[n]==nullptr)
            mexErrMsgIdAndTxt("Spinach:fsparse:argNull","Input argument pointer is null.");

    for (int n=0;n<2;n++)
    {
        if ((!mxIsDouble(prhs[n]))&&(mxGetClassID(prhs[n])!=mxINT32_CLASS))
            mexErrMsgIdAndTxt("Spinach:fsparse:indexType","row_idx and col_idx must be double or int32.");

        if (mxIsComplex(prhs[n]))
            mexErrMsgIdAndTxt("Spinach:fsparse:indexComplex","row_idx and col_idx must be real.");

        if (mxIsSparse(prhs[n]))
            mexErrMsgIdAndTxt("Spinach:fsparse:indexSparse","row_idx and col_idx must be full.");

        if (mxGetNumberOfDimensions(prhs[n])!=2)
            mexErrMsgIdAndTxt("Spinach:fsparse:indexShape","row_idx and col_idx must be two-dimensional arrays.");
    }

    if (mxGetNumberOfElements(prhs[0])!=mxGetNumberOfElements(prhs[1]))
        mexErrMsgIdAndTxt("Spinach:fsparse:indexCount","row_idx and col_idx must have the same number of elements.");

    if ((!mxIsDouble(prhs[2]))||mxIsSparse(prhs[2]))
        mexErrMsgIdAndTxt("Spinach:fsparse:valueType","vals must be a full double array.");

    if (mxGetNumberOfDimensions(prhs[2])!=2)
        mexErrMsgIdAndTxt("Spinach:fsparse:valueShape","vals must be a two-dimensional array.");

    const mwSize len=mxGetNumberOfElements(prhs[0]);
    const mwSize val_len=mxGetNumberOfElements(prhs[2]);

    if (!((val_len==len)||(val_len==1)||((val_len==0)&&(len==0))))
        mexErrMsgIdAndTxt("Spinach:fsparse:valueCount","vals must contain one value per index pair, or be scalar.");

    for (int n=3;n<5;n++)
    {
        if ((!mxIsDouble(prhs[n]))||mxIsComplex(prhs[n])||mxIsSparse(prhs[n])||
            (mxGetNumberOfElements(prhs[n])!=1))
            mexErrMsgIdAndTxt("Spinach:fsparse:dimType","n_rows and n_cols must be real double scalars.");
    }

    if (len>(mwSize)std::numeric_limits<int>::max())
        mexErrMsgIdAndTxt("Spinach:fsparse:indexLimit","Number of index pairs exceeds internal integer range.");
}

static int get_dim(const mxArray *array_ptr,const char *arg_name)
{
    const double dim_val=mxGetScalar(array_ptr);

    if ((!std::isfinite(dim_val))||(dim_val<0.0)||
        (dim_val!=std::floor(dim_val))||
        (dim_val>(double)std::numeric_limits<int>::max()))
        mexErrMsgIdAndTxt("Spinach:fsparse:dimValue","%s must be a non-negative integer in internal range.",arg_name);

    return (int)dim_val;
}

static int get_len(const mxArray *array_ptr)
{
    const mwSize len=mxGetNumberOfElements(array_ptr);

    if (len>(mwSize)std::numeric_limits<int>::max())
        mexErrMsgIdAndTxt("Spinach:fsparse:lenLimit","Number of index pairs exceeds internal integer range.");

    return (int)len;
}

static int *get_indices(const mxArray *array_ptr,int max_idx,const char *arg_name)
{
    const int len=get_len(array_ptr);
    int *idx=(int*)mxMalloc((len==0?1:len)*sizeof(int));

    if (idx==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:indexAlloc","Failed to allocate index work array.");

    if (mxIsDouble(array_ptr))
    {
        const double *idx_in=mxGetDoubles(array_ptr);

        if ((len!=0)&&(idx_in==nullptr))
            mexErrMsgIdAndTxt("Spinach:fsparse:indexAccess","Failed to access double index data.");

        for (int n=0;n<len;n++)
        {
            const double idx_val=idx_in[n];

            if ((!std::isfinite(idx_val))||(idx_val<1.0)||
                (idx_val!=std::floor(idx_val))||(idx_val>(double)max_idx))
                mexErrMsgIdAndTxt("Spinach:fsparse:indexValue","%s entries must be positive integers within matrix dimensions.",arg_name);

            idx[n]=(int)idx_val;
        }
    }
    else
    {
        const int32_T *idx_in=(const int32_T*)mxGetData(array_ptr);

        if ((len!=0)&&(idx_in==nullptr))
            mexErrMsgIdAndTxt("Spinach:fsparse:indexAccess","Failed to access int32 index data.");

        for (int n=0;n<len;n++)
        {
            if ((idx_in[n]<1)||(idx_in[n]>max_idx))
                mexErrMsgIdAndTxt("Spinach:fsparse:indexValue","%s entries must be positive integers within matrix dimensions.",arg_name);

            idx[n]=(int)idx_in[n];
        }
    }

    return idx;
}

static mxArray *call_sparse(int nrhs,const mxArray *prhs[])
{
    mxArray *lhs[1]={nullptr};
    mxArray *rhs[5];

    for (int n=0;n<nrhs;n++)
        rhs[n]=const_cast<mxArray*>(prhs[n]);

    if (mexCallMATLAB(1,lhs,nrhs,rhs,"sparse")!=0)
        mexErrMsgIdAndTxt("Spinach:fsparse:sparseFallback","Matlab sparse fallback failed.");

    return lhs[0];
}

static void squeeze(mxArray *S)
{
    const mwSize N=mxGetN(S);
    mwIndex *jcS=mxGetJc(S);
    mwIndex *irS=mxGetIr(S);
    double *prS=mxGetDoubles(S);

    if ((jcS==nullptr)||((mxGetNzmax(S)!=0)&&((irS==nullptr)||(prS==nullptr))))
        mexErrMsgIdAndTxt("Spinach:fsparse:sparseAccess","Failed to access sparse output arrays.");

    mwIndex dest=0;

    for (mwIndex col=0;col<(mwIndex)N;col++)
    {
        const mwIndex col_start=jcS[col];
        const mwIndex col_end=jcS[col+1];
        jcS[col]=dest;

        for (mwIndex ptr=col_start;ptr<col_end;ptr++)
        {
            if (prS[ptr]!=0.0)
            {
                irS[dest]=irS[ptr];
                prS[dest]=prS[ptr];
                dest++;
            }
        }
    }

    jcS[N]=dest;
}

static void sparse_insert(mwIndex *irS,double *prS,const int *irank,
                          const int *rank,const mwIndex *jrS,
                          const int *ii,const double *sr,
                          int smod,int sdiv,int len,int M)
{
    switch (2*(smod==len)+(sdiv==1))
    {
        case 3:
#ifndef _OPENMP
            for (int n=0;n<len;n++)
            {
                irS[irank[n]]=(mwIndex)(ii[n]-1);
                prS[irank[n]]+=sr[n];
            }
#else
            if (rank!=nullptr)
            {
#pragma omp parallel
                {
                    const int n_threads=omp_get_num_threads();
                    const int thread_id=omp_get_thread_num();
                    const int row_start=1+M*thread_id/n_threads;
                    const int row_end=M*(thread_id+1)/n_threads;
                    int insert_start;

                    if (row_start==1)
                        insert_start=0;
                    else
                        insert_start=(int)jrS[row_start-1];

                    if (row_end>=1)
                    {
                        for (int n=insert_start;n<(int)jrS[row_end];n++)
                            irS[irank[n]]=(mwIndex)(ii[rank[n]]-1);

                        for (int n=insert_start;n<(int)jrS[row_end];n++)
                            prS[irank[n]]+=sr[rank[n]];
                    }
                }
            }
            else
            {
#pragma omp parallel
                {
#pragma omp single nowait
                    for (int n=0;n<len;n++)
                        irS[irank[n]]=(mwIndex)(ii[n]-1);

#pragma omp single nowait
                    for (int n=0;n<len;n++)
                        prS[irank[n]]+=sr[n];
                }
            }
#endif
            break;

        case 0:
        {
            const double scalar=sr[0];
#ifndef _OPENMP
            for (int n=0;n<len;n++)
            {
                irS[irank[n]]=(mwIndex)(ii[n]-1);
                prS[irank[n]]+=scalar;
            }
#else
#pragma omp parallel
            {
#pragma omp single nowait
                for (int n=0;n<len;n++)
                    irS[irank[n]]=(mwIndex)(ii[n]-1);

#pragma omp single nowait
                for (int n=0;n<len;n++)
                    prS[irank[n]]+=scalar;
            }
#endif
            break;
        }

        default:
            mexErrMsgIdAndTxt("Spinach:fsparse:valueShape","Unsupported value shape.");
    }
}

#ifndef _OPENMP
static mxArray *sparse_real(const int *ii,const int *jj,const double *sr,
                            int smod,int sdiv,int len,int M,int N)
{
    if (len==0)
        return mxCreateSparse((mwSize)M,(mwSize)N,0,mxREAL);

    // Count and accumulate indices by row
    mwIndex *jrS=(mwIndex*)mxCalloc((mwSize)M+1,sizeof(mwIndex));

    if (jrS==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:rowAlloc","Failed to allocate row counters.");

    for (int n=0;n<len;n++)
        jrS[ii[n]]++;

    for (int row=2;row<=M;row++)
        jrS[row]+=jrS[row-1];

    // Build a row-wise rank array
    int *rank=(int*)mxMalloc((mwSize)len*sizeof(int));

    if (rank==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:rankAlloc","Failed to allocate rank array.");

    jrS--;
    for (int n=0;n<len;n++)
        rank[jrS[ii[n]]++]=n;

    // Count unique row entries within every output column
    mwIndex *jcS=(mwIndex*)mxCalloc((mwSize)N+1,sizeof(mwIndex));
    int *hcol=(int*)mxCalloc((mwSize)N,sizeof(int));
    int *irank=(int*)mxMalloc((mwSize)len*sizeof(int));

    if ((jcS==nullptr)||(hcol==nullptr)||(irank==nullptr))
        mexErrMsgIdAndTxt("Spinach:fsparse:workAlloc","Failed to allocate sparse assembly work arrays.");

    hcol--;
    for (int row=1,pos=0;row<=M;row++)
    {
        for (;pos<(int)jrS[row];pos++)
        {
            const int triplet=rank[pos];
            const int col=jj[triplet];

            if (hcol[col]<row)
            {
                hcol[col]=row;
                jcS[col]++;
            }

            irank[triplet]=(int)jcS[col]-1;
        }
    }

    mxFree(++hcol);
    mxFree(rank);
    mxFree(++jrS);

    // Accumulate column pointers
    for (int col=2;col<=N;col++)
        jcS[col]+=jcS[col-1];

    jcS--;
    for (int n=0;n<len;n++)
        irank[n]+=(int)jcS[jj[n]];
    jcS++;

    // Allocate output and install exact column pointers
    mxArray *S=mxCreateSparse(0,0,jcS[N],mxREAL);
    mxSetM(S,(mwSize)M);
    mxSetN(S,(mwSize)N);
    mxFree(mxGetJc(S));
    mxSetJc(S,jcS);

    // Insert values into the sparse matrix
    sparse_insert(mxGetIr(S),mxGetDoubles(S),irank,nullptr,nullptr,ii,sr,
                  smod,sdiv,len,M);

    mxFree(irank);
    return S;
}
#else
static mxArray *sparse_real(const int *ii,const int *jj,const double *sr,
                            int smod,int sdiv,int len,int M,int N)
{
    if (len==0)
        return mxCreateSparse((mwSize)M,(mwSize)N,0,mxREAL);

    // Count and accumulate indices by row on each thread
    omp_set_dynamic(0);
    const int n_threads=omp_get_max_threads();
    mwIndex **jrS=(mwIndex**)mxMalloc(((mwSize)n_threads+1)*sizeof(mwIndex*));

    if (jrS==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:rowAlloc","Failed to allocate row counter table.");

    for (int k=0;k<=n_threads;k++)
    {
        jrS[k]=(mwIndex*)mxCalloc((mwSize)M+1,sizeof(mwIndex));

        if (jrS[k]==nullptr)
            mexErrMsgIdAndTxt("Spinach:fsparse:rowAlloc","Failed to allocate row counters.");

        jrS[k]--;
    }

#pragma omp parallel num_threads(n_threads)
    {
        const int thread_id=omp_get_thread_num();
        const int idx_start=len*thread_id/n_threads;
        const int idx_end=len*(thread_id+1)/n_threads;

        for (int n=idx_start;n<idx_end;n++)
            jrS[thread_id+1][ii[n]]++;

#pragma omp barrier

#pragma omp for
        for (int row=1;row<=M;row++)
            for (int k=1;k<n_threads;k++)
                jrS[k+1][row]+=jrS[k][row];

#pragma omp single
        for (int row=1;row<=M;row++)
            jrS[0][row+1]+=jrS[0][row]+jrS[n_threads][row];

#pragma omp for
        for (int row=1;row<=M;row++)
            for (int k=1;k<n_threads;k++)
                jrS[k][row]+=jrS[0][row];
    }

    // Build a row-wise rank array
    int *rank=(int*)mxMalloc((mwSize)len*sizeof(int));

    if (rank==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:rankAlloc","Failed to allocate rank array.");

#pragma omp parallel num_threads(n_threads)
    {
        const int thread_id=omp_get_thread_num();
        const int idx_start=len*thread_id/n_threads;
        const int idx_end=len*(thread_id+1)/n_threads;

        for (int n=idx_start;n<idx_end;n++)
            rank[jrS[thread_id][ii[n]]++]=n;
    }

    // Count unique row entries within every output column
    mwIndex **jcS=(mwIndex**)mxMalloc(((mwSize)n_threads+1)*sizeof(mwIndex*));

    if (jcS==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:colAlloc","Failed to allocate column pointer table.");

    for (int k=0;k<=n_threads;k++)
    {
        jcS[k]=(mwIndex*)mxCalloc((mwSize)N+1,sizeof(mwIndex));

        if (jcS[k]==nullptr)
            mexErrMsgIdAndTxt("Spinach:fsparse:colAlloc","Failed to allocate column counters.");
    }

    int *irankp=(int*)mxMalloc((mwSize)len*sizeof(int));
    int *irank=nullptr;

    if (irankp==nullptr)
        mexErrMsgIdAndTxt("Spinach:fsparse:rankAlloc","Failed to allocate inverse rank array.");

    if (2*(smod==len)+(sdiv==1)!=3)
    {
        irank=(int*)mxMalloc((mwSize)len*sizeof(int));

        if (irank==nullptr)
            mexErrMsgIdAndTxt("Spinach:fsparse:rankAlloc","Failed to allocate inverse rank array.");
    }

#pragma omp parallel num_threads(n_threads)
    {
        int *hcol;

#pragma omp critical
        hcol=(int*)mxCalloc((mwSize)N,sizeof(int));

        if (hcol==nullptr)
            mexErrMsgIdAndTxt("Spinach:fsparse:hashAlloc","Failed to allocate column marker array.");

        hcol--;

        const int thread_id=omp_get_thread_num();
        const int row_start=1+M*thread_id/n_threads;
        const int row_end=M*(thread_id+1)/n_threads;
        int insert_start;

        if (row_start==1)
            insert_start=0;
        else
            insert_start=(int)jrS[n_threads-1][row_start-1];

        for (int row=row_start,pos=insert_start;row<=row_end;row++)
        {
            for (;pos<(int)jrS[n_threads-1][row];pos++)
            {
                const int triplet=rank[pos];
                const int col=jj[triplet];

                if (hcol[col]<row)
                {
                    hcol[col]=row;
                    jcS[thread_id+1][col]++;
                }

                irankp[pos]=(int)jcS[thread_id+1][col]-1;
            }
        }

#pragma omp critical
        mxFree(++hcol);

#pragma omp barrier

#pragma omp for
        for (int col=1;col<=N;col++)
            for (int k=1;k<n_threads;k++)
                jcS[k+1][col]+=jcS[k][col];

#pragma omp single
        {
            for (int col=1;col<=N;col++)
                jcS[0][col]+=jcS[0][col-1]+jcS[n_threads][col];

            jcS[0]--;
        }

#pragma omp for
        for (int col=1;col<=N;col++)
            for (int k=1;k<n_threads;k++)
                jcS[k][col]+=jcS[0][col];

        if (row_end>=1)
            for (int pos=insert_start;pos<(int)jrS[n_threads-1][row_end];pos++)
                irankp[pos]+=(int)jcS[thread_id][jj[rank[pos]]];

        if (irank!=nullptr)
            if (row_end>=1)
                for (int pos=insert_start;pos<(int)jrS[n_threads-1][row_end];pos++)
                    irank[rank[pos]]=irankp[pos];
    }

    // Allocate output and install exact column pointers
    jcS[0]++;
    mxArray *S=mxCreateSparse(0,0,jcS[0][N],mxREAL);
    mxSetM(S,(mwSize)M);
    mxSetN(S,(mwSize)N);
    mxFree(mxGetJc(S));
    mxSetJc(S,jcS[0]);

    for (int k=1;k<=n_threads;k++)
        mxFree(jcS[k]);

    mxFree(jcS);

    // Insert values into the sparse matrix
    if (irank==nullptr)
        sparse_insert(mxGetIr(S),mxGetDoubles(S),irankp,rank,
                      jrS[n_threads-1],ii,sr,smod,sdiv,len,M);
    else
        sparse_insert(mxGetIr(S),mxGetDoubles(S),irank,nullptr,nullptr,ii,sr,
                      smod,sdiv,len,M);

    // Release work arrays
    if (irank!=nullptr)
        mxFree(irank);

    mxFree(irankp);
    mxFree(rank);

    for (int k=0;k<=n_threads;k++)
        mxFree(++jrS[k]);

    mxFree(jrS);

    return S;
}
#endif
