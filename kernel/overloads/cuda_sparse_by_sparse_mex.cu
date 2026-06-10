/* cuda_sparse_by_sparse_mex.cu
 *
 * Sparse GPU matrix product through cuSPARSE SpGEMM ALG3
 *
 * Syntax:
 *
 *      [row_c,col_c,val_c]=cuda_sparse_by_sparse_mex(row_a,col_a,val_a,...
 *                                                    row_b,col_b,val_b,...
 *                                                    dims,chunk_fraction)
 *
 * Internal helper for cuda_sparse_by_sparse.m. Inputs row_a and row_b
 * are zero-based int32 COO row-index gpuArrays. Inputs col_a and col_b
 * are zero-based int32 COO column-index gpuArrays. Inputs val_a and
 * val_b are real or complex double gpuArrays in MATLAB find() order. The
 * MEX gateway interprets the column-compressed MATLAB ordering as CSR
 * storage of the transposed matrices, computes B.'*A.' with cuSPARSE
 * SpGEMM ALG3, and returns one-based triplets for C=A*B. dims is
 * uint64([rows_a cols_a cols_b]). Output row_c and col_c are one-based
 * int32 gpuArray indices, and val_c is a real or complex double gpuArray.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cusparse.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

class DeviceBuffer
{
public:
    DeviceBuffer()=default;
    explicit DeviceBuffer(size_t n_bytes) { allocate(n_bytes); }
    DeviceBuffer(const DeviceBuffer&)=delete;
    DeviceBuffer& operator=(const DeviceBuffer&)=delete;
    ~DeviceBuffer() { release(); }

    void allocate(size_t n_bytes)
    {
        release();
        if (n_bytes>0)
        {
            cudaError_t status=cudaMalloc(&ptr,n_bytes);
            if (status!=cudaSuccess)
                throw std::runtime_error(std::string("cudaMalloc failed: ")+cudaGetErrorString(status));
        }
    }

    void release()
    {
        if (ptr!=nullptr)
        {
            cudaFree(ptr);
            ptr=nullptr;
        }
    }

    void* data() { return ptr; }

    template<typename T> T* as()
    {
        return static_cast<T*>(ptr);
    }

private:
    void *ptr=nullptr;
};

class GpuView
{
public:
    explicit GpuView(const mxArray *array)
    {
        ptr=mxGPUCreateFromMxArray(array);
        if (ptr==nullptr)
            throw std::runtime_error("Failed to create an mxGPUArray view.");
    }

    GpuView(const GpuView&)=delete;
    GpuView& operator=(const GpuView&)=delete;
    ~GpuView()
    {
        if (ptr!=nullptr)
            mxGPUDestroyGPUArray(ptr);
    }

    const mxGPUArray* get() const { return ptr; }

private:
    const mxGPUArray *ptr=nullptr;
};

class GpuOutput
{
public:
    GpuOutput(mwSize n_elem,mxClassID class_id,mxComplexity complexity)
    {
        const mwSize dims[2]={n_elem,1};
        ptr=mxGPUCreateGPUArray(2,dims,class_id,complexity,MX_GPU_DO_NOT_INITIALIZE);
        if (ptr==nullptr)
            throw std::runtime_error("Failed to create an output gpuArray.");
    }

    GpuOutput(const GpuOutput&)=delete;
    GpuOutput& operator=(const GpuOutput&)=delete;
    ~GpuOutput()
    {
        if (ptr!=nullptr)
            mxGPUDestroyGPUArray(ptr);
    }

    template<typename T> T* data()
    {
        return static_cast<T*>(mxGPUGetData(ptr));
    }

    mxArray* to_mx() const
    {
        return mxGPUCreateMxArrayOnGPU(ptr);
    }

private:
    mxGPUArray *ptr=nullptr;
};

class CusparseHandle
{
public:
    CusparseHandle()
    {
        cusparseStatus_t status=cusparseCreate(&handle);
        if (status!=CUSPARSE_STATUS_SUCCESS)
            throw std::runtime_error(std::string("cusparseCreate failed: ")+cusparseGetErrorString(status));
    }

    CusparseHandle(const CusparseHandle&)=delete;
    CusparseHandle& operator=(const CusparseHandle&)=delete;
    ~CusparseHandle()
    {
        if (handle!=nullptr)
            cusparseDestroy(handle);
    }

    cusparseHandle_t get() const { return handle; }

private:
    cusparseHandle_t handle=nullptr;
};

class SpMat
{
public:
    ~SpMat()
    {
        if (descr!=nullptr)
            cusparseDestroySpMat(descr);
    }

    cusparseSpMatDescr_t *ptr() { return &descr; }
    operator cusparseSpMatDescr_t() const { return descr; }

private:
    cusparseSpMatDescr_t descr=nullptr;
};

class SpGemmDesc
{
public:
    SpGemmDesc()
    {
        cusparseStatus_t status=cusparseSpGEMM_createDescr(&descr);
        if (status!=CUSPARSE_STATUS_SUCCESS)
            throw std::runtime_error(std::string("cusparseSpGEMM_createDescr failed: ")+cusparseGetErrorString(status));
    }

    SpGemmDesc(const SpGemmDesc&)=delete;
    SpGemmDesc& operator=(const SpGemmDesc&)=delete;
    ~SpGemmDesc()
    {
        if (descr!=nullptr)
            cusparseSpGEMM_destroyDescr(descr);
    }

    operator cusparseSpGEMMDescr_t() const { return descr; }

private:
    cusparseSpGEMMDescr_t descr=nullptr;
};

struct InputCsrMatrix
{
    int rows=0;
    int cols=0;
    int nnz=0;
    const void *row_offsets=nullptr;
    const void *col_indices=nullptr;
    const void *values=nullptr;
    DeviceBuffer row_store;
    DeviceBuffer value_store;
};

struct OutputCsrMatrix
{
    int rows=0;
    int cols=0;
    int nnz=0;
    DeviceBuffer row_offsets;
    DeviceBuffer col_indices;
    DeviceBuffer values;
};

static void grumble(int nlhs,int nrhs,const mxArray *prhs[]);
static void check_cuda(cudaError_t status,const char *call);
static void check_cusparse(cusparseStatus_t status,const char *call);
static void validate_gpu_view(const GpuView &view,const char *name,
                              mxClassID class_id,bool allow_complex);
static int get_nnz(const GpuView &view);
static void prepare_transpose_input(cusparseHandle_t handle,
                                    const GpuView &row_idx,
                                    const GpuView &col_idx,
                                    const GpuView &vals,int n_rows,
                                    int n_cols,bool make_complex,
                                    InputCsrMatrix &mat);
static void run_spgemm(cusparseHandle_t handle,InputCsrMatrix &mat_a,
                       InputCsrMatrix &mat_b,int n_rows,int n_cols,
                       cudaDataType value_type,const void *alpha,
                       const void *beta,float chunk_fraction,
                       OutputCsrMatrix &mat_c);
static void make_outputs(OutputCsrMatrix &mat_c,bool make_complex,
                         mxArray *plhs[]);

__global__ void real_to_complex_kernel(const double *vals_in,
                                       cuDoubleComplex *vals_out,int nnz)
{
    const int idx=blockIdx.x*blockDim.x+threadIdx.x;
    if (idx<nnz)
        vals_out[idx]=make_cuDoubleComplex(vals_in[idx],0.0);
}

__global__ void expand_rows_kernel(const int *row_offsets,int *rows_out,
                                   int n_rows)
{
    const int row=blockIdx.x*blockDim.x+threadIdx.x;
    if (row>=n_rows)
        return;

    const int start=row_offsets[row];
    const int finish=row_offsets[row+1];

    for (int ptr=start;ptr<finish;ptr++)
        rows_out[ptr]=row+1;
}

__global__ void one_based_cols_kernel(const int *cols_in,int *cols_out,int nnz)
{
    const int idx=blockIdx.x*blockDim.x+threadIdx.x;
    if (idx<nnz)
        cols_out[idx]=cols_in[idx]+1;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    try
    {
        if (mxInitGPU()!=MX_GPU_SUCCESS)
            throw std::runtime_error("Failed to initialise the MATLAB GPU API.");

        grumble(nlhs,nrhs,prhs);

        const GpuView row_a(prhs[0]);
        const GpuView col_a(prhs[1]);
        const GpuView val_a(prhs[2]);
        const GpuView row_b(prhs[3]);
        const GpuView col_b(prhs[4]);
        const GpuView val_b(prhs[5]);

        validate_gpu_view(row_a,"row_a",mxINT32_CLASS,false);
        validate_gpu_view(col_a,"col_a",mxINT32_CLASS,false);
        validate_gpu_view(val_a,"val_a",mxDOUBLE_CLASS,true);
        validate_gpu_view(row_b,"row_b",mxINT32_CLASS,false);
        validate_gpu_view(col_b,"col_b",mxINT32_CLASS,false);
        validate_gpu_view(val_b,"val_b",mxDOUBLE_CLASS,true);
        check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");

        const uint64_T *dims=static_cast<const uint64_T*>(mxGetData(prhs[6]));
        const uint64_T max_int=static_cast<uint64_T>(std::numeric_limits<int>::max());

        if ((dims[0]>max_int)||(dims[1]>max_int)||(dims[2]>max_int))
            throw std::runtime_error("Matrix dimensions must fit into int32.");

        const int n_rows=static_cast<int>(dims[0]);
        const int n_inner=static_cast<int>(dims[1]);
        const int n_cols=static_cast<int>(dims[2]);
        const float chunk_fraction=static_cast<float>(mxGetScalar(prhs[7]));
        const bool make_complex=(mxGPUGetComplexity(val_a.get())==mxCOMPLEX)||
                                (mxGPUGetComplexity(val_b.get())==mxCOMPLEX);

        InputCsrMatrix mat_a;
        InputCsrMatrix mat_b;
        OutputCsrMatrix mat_c;
        CusparseHandle handle;

        prepare_transpose_input(handle.get(),row_a,col_a,val_a,
                                n_rows,n_inner,make_complex,mat_a);
        prepare_transpose_input(handle.get(),row_b,col_b,val_b,
                                n_inner,n_cols,make_complex,mat_b);

        if (make_complex)
        {
            const cuDoubleComplex alpha=make_cuDoubleComplex(1.0,0.0);
            const cuDoubleComplex beta=make_cuDoubleComplex(0.0,0.0);

            run_spgemm(handle.get(),mat_b,mat_a,n_cols,n_rows,CUDA_C_64F,
                       &alpha,&beta,chunk_fraction,mat_c);
        }
        else
        {
            const double alpha=1.0;
            const double beta=0.0;

            run_spgemm(handle.get(),mat_b,mat_a,n_cols,n_rows,CUDA_R_64F,
                       &alpha,&beta,chunk_fraction,mat_c);
        }

        make_outputs(mat_c,make_complex,plhs);
        check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");
    }
    catch (const std::exception &err)
    {
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:runtime","%s",err.what());
    }
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=8)
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:nrhs","Eight inputs are required.");

    if (nlhs!=3)
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:nlhs","Three outputs are required.");

    if (prhs==nullptr)
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:prhsNull","Input argument array pointer is null.");

    for (int arg=0;arg<6;arg++)
        if (!mxIsGPUArray(prhs[arg]))
            mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:notGpu","COO arrays must be gpuArrays.");

    if ((mxGetClassID(prhs[6])!=mxUINT64_CLASS)||(mxGetNumberOfElements(prhs[6])!=3))
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:dims","dims must be a uint64 vector with three elements.");

    if (!mxIsDouble(prhs[7])||mxIsComplex(prhs[7])||(mxGetNumberOfElements(prhs[7])!=1))
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:chunk","chunk_fraction must be a real double scalar.");

    const double chunk_fraction=mxGetScalar(prhs[7]);

    if ((!std::isfinite(chunk_fraction))||(chunk_fraction<=0.0)||(chunk_fraction>1.0))
        mexErrMsgIdAndTxt("Spinach:cuda_sparse_by_sparse_mex:chunkRange","chunk_fraction must be in the range (0,1].");
}

static void check_cuda(cudaError_t status,const char *call)
{
    if (status!=cudaSuccess)
        throw std::runtime_error(std::string(call)+" failed: "+cudaGetErrorString(status));
}

static void check_cusparse(cusparseStatus_t status,const char *call)
{
    if (status!=CUSPARSE_STATUS_SUCCESS)
        throw std::runtime_error(std::string(call)+" failed: "+cusparseGetErrorString(status));
}

static void validate_gpu_view(const GpuView &view,const char *name,
                              mxClassID class_id,bool allow_complex)
{
    if (mxGPUIsSparse(view.get()))
        throw std::runtime_error(std::string(name)+" must be full.");

    if (mxGPUGetClassID(view.get())!=class_id)
        throw std::runtime_error(std::string(name)+" has an incompatible underlying type.");

    if ((!allow_complex)&&(mxGPUGetComplexity(view.get())!=mxREAL))
        throw std::runtime_error(std::string(name)+" must be real.");
}

static int get_nnz(const GpuView &view)
{
    const mwSize n_elem=mxGPUGetNumberOfElements(view.get());
    const mwSize max_int=static_cast<mwSize>(std::numeric_limits<int>::max());

    if (n_elem>max_int)
        throw std::runtime_error("COO arrays must have int32-compatible lengths.");

    return static_cast<int>(n_elem);
}

static void prepare_transpose_input(cusparseHandle_t handle,
                                    const GpuView &row_idx,
                                    const GpuView &col_idx,
                                    const GpuView &vals,int n_rows,
                                    int n_cols,bool make_complex,
                                    InputCsrMatrix &mat)
{
    mat.rows=n_cols;
    mat.cols=n_rows;
    mat.nnz=get_nnz(vals);

    if (get_nnz(row_idx)!=mat.nnz)
        throw std::runtime_error("COO row-index and value arrays have inconsistent lengths.");

    if (get_nnz(col_idx)!=mat.nnz)
        throw std::runtime_error("COO column-index and value arrays have inconsistent lengths.");

    mat.row_store.allocate(static_cast<size_t>(n_cols+1)*sizeof(int));

    if (mat.nnz>0)
    {
        check_cusparse(cusparseXcoo2csr(handle,
                                        static_cast<const int*>(mxGPUGetDataReadOnly(col_idx.get())),
                                        mat.nnz,n_cols,mat.row_store.as<int>(),
                                        CUSPARSE_INDEX_BASE_ZERO),
                       "cusparseXcoo2csr");
    }
    else
    {
        check_cuda(cudaMemset(mat.row_store.data(),0,
                              static_cast<size_t>(n_cols+1)*sizeof(int)),
                   "cudaMemset");
    }

    mat.row_offsets=mat.row_store.data();
    mat.col_indices=mxGPUGetDataReadOnly(row_idx.get());

    if (make_complex&&(mxGPUGetComplexity(vals.get())==mxREAL))
    {
        mat.value_store.allocate(static_cast<size_t>(mat.nnz)*sizeof(cuDoubleComplex));

        if (mat.nnz>0)
        {
            const int block_size=256;
            const int grid_size=(mat.nnz+block_size-1)/block_size;

            real_to_complex_kernel<<<grid_size,block_size>>>(
                static_cast<const double*>(mxGPUGetDataReadOnly(vals.get())),
                mat.value_store.as<cuDoubleComplex>(),mat.nnz);
            check_cuda(cudaGetLastError(),"real_to_complex_kernel");
        }

        mat.values=mat.value_store.data();
    }
    else
    {
        mat.values=mxGPUGetDataReadOnly(vals.get());
    }
    check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");
}

static void run_spgemm(cusparseHandle_t handle,InputCsrMatrix &mat_a,
                       InputCsrMatrix &mat_b,int n_rows,int n_cols,
                       cudaDataType value_type,const void *alpha,
                       const void *beta,float chunk_fraction,
                       OutputCsrMatrix &mat_c)
{
    mat_c.rows=n_rows;
    mat_c.cols=n_cols;
    mat_c.row_offsets.allocate(static_cast<size_t>(n_rows+1)*sizeof(int));

    SpMat descr_a;
    SpMat descr_b;
    SpMat descr_c;
    SpGemmDesc spgemm_descr;

    check_cusparse(cusparseCreateCsr(descr_a.ptr(),mat_a.rows,mat_a.cols,mat_a.nnz,
                                     const_cast<void*>(mat_a.row_offsets),
                                     const_cast<void*>(mat_a.col_indices),
                                     const_cast<void*>(mat_a.values),
                                     CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO,value_type),
                   "cusparseCreateCsr");
    check_cusparse(cusparseCreateCsr(descr_b.ptr(),mat_b.rows,mat_b.cols,mat_b.nnz,
                                     const_cast<void*>(mat_b.row_offsets),
                                     const_cast<void*>(mat_b.col_indices),
                                     const_cast<void*>(mat_b.values),
                                     CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO,value_type),
                   "cusparseCreateCsr");
    check_cusparse(cusparseCreateCsr(descr_c.ptr(),n_rows,n_cols,0,
                                     mat_c.row_offsets.data(),nullptr,nullptr,
                                     CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO,value_type),
                   "cusparseCreateCsr");

    size_t buffer_size1=0;
    size_t buffer_size2=0;
    size_t buffer_size3=0;

    check_cusparse(cusparseSpGEMM_workEstimation(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 alpha,descr_a,descr_b,beta,descr_c,
                                                 value_type,CUSPARSE_SPGEMM_ALG3,
                                                 spgemm_descr,&buffer_size1,nullptr),
                   "cusparseSpGEMM_workEstimation");

    DeviceBuffer buffer1(buffer_size1);

    check_cusparse(cusparseSpGEMM_workEstimation(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 alpha,descr_a,descr_b,beta,descr_c,
                                                 value_type,CUSPARSE_SPGEMM_ALG3,
                                                 spgemm_descr,&buffer_size1,buffer1.data()),
                   "cusparseSpGEMM_workEstimation");

    int64_t num_products=0;
    check_cusparse(cusparseSpGEMM_getNumProducts(spgemm_descr,&num_products),
                   "cusparseSpGEMM_getNumProducts");

    if (num_products<0)
        throw std::runtime_error("cuSPARSE reported a negative intermediate product count.");

    check_cusparse(cusparseSpGEMM_estimateMemory(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 alpha,descr_a,descr_b,beta,descr_c,
                                                 value_type,CUSPARSE_SPGEMM_ALG3,
                                                 spgemm_descr,chunk_fraction,
                                                 &buffer_size3,nullptr,nullptr),
                   "cusparseSpGEMM_estimateMemory");

    DeviceBuffer buffer3(buffer_size3);

    check_cusparse(cusparseSpGEMM_estimateMemory(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 alpha,descr_a,descr_b,beta,descr_c,
                                                 value_type,CUSPARSE_SPGEMM_ALG3,
                                                 spgemm_descr,chunk_fraction,
                                                 &buffer_size3,buffer3.data(),&buffer_size2),
                   "cusparseSpGEMM_estimateMemory");

    buffer3.release();

    DeviceBuffer buffer2(buffer_size2);

    check_cusparse(cusparseSpGEMM_compute(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          alpha,descr_a,descr_b,beta,descr_c,
                                          value_type,CUSPARSE_SPGEMM_ALG3,
                                          spgemm_descr,&buffer_size2,buffer2.data()),
                   "cusparseSpGEMM_compute");

    int64_t c_rows=0;
    int64_t c_cols=0;
    int64_t c_nnz=0;

    check_cusparse(cusparseSpMatGetSize(descr_c,&c_rows,&c_cols,&c_nnz),
                   "cusparseSpMatGetSize");

    if ((c_rows!=n_rows)||(c_cols!=n_cols))
        throw std::runtime_error("cuSPARSE returned inconsistent output dimensions.");

    if (c_nnz>static_cast<int64_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("cuSPARSE output nonzero count exceeds int32.");

    mat_c.nnz=static_cast<int>(c_nnz);

    if (mat_c.nnz==0)
        return;

    const size_t value_size=(value_type==CUDA_R_64F)?sizeof(double):sizeof(cuDoubleComplex);

    mat_c.col_indices.allocate(static_cast<size_t>(mat_c.nnz)*sizeof(int));
    mat_c.values.allocate(static_cast<size_t>(mat_c.nnz)*value_size);

    check_cusparse(cusparseCsrSetPointers(descr_c,mat_c.row_offsets.data(),
                                          mat_c.col_indices.data(),mat_c.values.data()),
                   "cusparseCsrSetPointers");
    check_cusparse(cusparseSpGEMM_copy(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       alpha,descr_a,descr_b,beta,descr_c,
                                       value_type,CUSPARSE_SPGEMM_ALG3,
                                       spgemm_descr),
                   "cusparseSpGEMM_copy");
    check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");

    return;
}

static void make_outputs(OutputCsrMatrix &mat_c,bool make_complex,
                         mxArray *plhs[])
{
    GpuOutput row_out(mat_c.nnz,mxINT32_CLASS,mxREAL);
    GpuOutput col_out(mat_c.nnz,mxINT32_CLASS,mxREAL);
    GpuOutput val_out(mat_c.nnz,mxDOUBLE_CLASS,make_complex?mxCOMPLEX:mxREAL);

    if (mat_c.nnz>0)
    {
        const int block_size=256;
        const int row_grid=(mat_c.rows+block_size-1)/block_size;

        const int col_grid=(mat_c.nnz+block_size-1)/block_size;

        one_based_cols_kernel<<<col_grid,block_size>>>(
            mat_c.col_indices.as<int>(),row_out.data<int>(),mat_c.nnz);
        check_cuda(cudaGetLastError(),"one_based_cols_kernel");

        expand_rows_kernel<<<row_grid,block_size>>>(
            mat_c.row_offsets.as<int>(),col_out.data<int>(),mat_c.rows);
        check_cuda(cudaGetLastError(),"expand_rows_kernel");

        const size_t value_size=make_complex?sizeof(cuDoubleComplex):sizeof(double);
        check_cuda(cudaMemcpy(val_out.data<void>(),mat_c.values.data(),
                              static_cast<size_t>(mat_c.nnz)*value_size,
                              cudaMemcpyDeviceToDevice),"cudaMemcpy");
    }

    plhs[0]=row_out.to_mx();
    plhs[1]=col_out.to_mx();
    plhs[2]=val_out.to_mx();
}
