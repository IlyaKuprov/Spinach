/* cuda_sparse_by_sparse_mex.cu
 *
 * Sparse GPU matrix product through cuSPARSE SpGEMM, reading the operands
 * from MATLAB's own sparse gpuArray storage
 *
 * Syntax:
 *
 *      [row_c,col_c,val_c]=cuda_sparse_by_sparse_mex(A,B,alg)
 *
 * Internal helper for cuda_sparse_by_sparse.m. Inputs A and B are sparse
 * double gpuArrays with size(A,2)==size(B,1); alg is 1, 2, or 3 and selects
 * CUSPARSE_SPGEMM_ALG1, ALG2, or ALG3. Outputs row_c and col_c are one-based
 * int32 index column gpuArrays and val_c is a double column gpuArray, all of
 * length nnz(A*B) in row-major order; as in native mtimes, val_c is complex
 * if either input is complex and neither has zero stored entries.
 *
 * MATLAB R2026b keeps a sparse gpuArray as zero-based CSR of the matrix
 * itself with int32 indices. The storage block is two pointer hops from the
 * mxGPUArray handle and holds the value, column index, and row offset device
 * pointers at byte offsets 144, 152, and 160. The layout is undocumented, so
 * before use the MEX checks that all three are device allocations on the
 * current GPU large enough for the matrix, that the first row offset is zero,
 * and that the last equals nzmax(), the stored entry count reported by
 * MATLAB; it exceeds nnz() when GPU arithmetic leaves explicit zeros. Any
 * mismatch raises Spinach:cuda_sparse_by_sparse_mex:layout. Input buffers
 * are never written.
 *
 * ALG1 and ALG2 run on MATLAB's 32-bit indices in place. ALG3 corrupts
 * device memory with 32-bit indices on large products, so its operands get
 * 64-bit index copies; ALG3 chunk fraction is halved from one until its
 * estimation and compute buffers fit into free device memory. When cuSPARSE
 * or CUDA runs out of resources, the rows of A are bisected at half of their
 * nonzeros and the halves are multiplied separately; a row block of A is
 * passed as pointer offsets into MATLAB's index and value arrays with a
 * rebased copy of its row offsets. A real operand of a mixed real-complex
 * product is copied into complex values, because cuSPARSE needs a single
 * value type.
 *
 * ilya.kuprov@weizmann.ac.il
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Resource exhaustion, recoverable by multiplying fewer rows of A
struct OutOfResources : std::runtime_error
{
    using std::runtime_error::runtime_error;
};

// MATLAB sparse gpuArray storage does not have the expected layout
struct LayoutMismatch : std::runtime_error
{
    using std::runtime_error::runtime_error;
};

static void check_cuda(cudaError_t status,const char *call)
{
    if (status==cudaSuccess) return;
    cudaGetLastError();
    const std::string message=std::string(call)+" failed: "+cudaGetErrorString(status);
    if (status==cudaErrorMemoryAllocation) throw OutOfResources(message);
    throw std::runtime_error(message);
}

static void check_cusparse(cusparseStatus_t status,const char *call)
{
    if (status==CUSPARSE_STATUS_SUCCESS) return;
    const std::string message=std::string(call)+" failed: "+cusparseGetErrorString(status);
    if ((status==CUSPARSE_STATUS_INSUFFICIENT_RESOURCES)||
        (status==CUSPARSE_STATUS_ALLOC_FAILED)) throw OutOfResources(message);
    throw std::runtime_error(message);
}

class DeviceBuffer
{
public:
    DeviceBuffer()=default;
    explicit DeviceBuffer(size_t n_bytes)
    {
        if (n_bytes>0) check_cuda(cudaMalloc(&ptr,n_bytes),"cudaMalloc");
    }
    DeviceBuffer(DeviceBuffer &&other) noexcept : ptr(other.ptr) { other.ptr=nullptr; }
    DeviceBuffer& operator=(DeviceBuffer &&other) noexcept
    {
        std::swap(ptr,other.ptr); return *this;
    }
    DeviceBuffer(const DeviceBuffer&)=delete;
    DeviceBuffer& operator=(const DeviceBuffer&)=delete;
    ~DeviceBuffer() { if (ptr!=nullptr) cudaFree(ptr); }
    template<typename T> T* as() const { return static_cast<T*>(ptr); }
private:
    void *ptr=nullptr;
};

class CusparseHandle
{
public:
    CusparseHandle() { check_cusparse(cusparseCreate(&handle),"cusparseCreate"); }
    CusparseHandle(const CusparseHandle&)=delete;
    CusparseHandle& operator=(const CusparseHandle&)=delete;
    ~CusparseHandle() { cusparseDestroy(handle); }
    operator cusparseHandle_t() const { return handle; }
private:
    cusparseHandle_t handle=nullptr;
};

class SpMat
{
public:
    SpMat(int64_t rows,int64_t cols,int64_t nnz,const void *offsets,const void *indices,
          const void *values,cusparseIndexType_t index_type,cudaDataType value_type)
    {
        check_cusparse(cusparseCreateCsr(&descr,rows,cols,nnz,const_cast<void*>(offsets),
                                         const_cast<void*>(indices),const_cast<void*>(values),
                                         index_type,index_type,CUSPARSE_INDEX_BASE_ZERO,
                                         value_type),"cusparseCreateCsr");
    }
    SpMat(const SpMat&)=delete;
    SpMat& operator=(const SpMat&)=delete;
    ~SpMat() { cusparseDestroySpMat(descr); }
    operator cusparseSpMatDescr_t() const { return descr; }
private:
    cusparseSpMatDescr_t descr=nullptr;
};

class SpGemmDesc
{
public:
    SpGemmDesc() { check_cusparse(cusparseSpGEMM_createDescr(&descr),"cusparseSpGEMM_createDescr"); }
    SpGemmDesc(const SpGemmDesc&)=delete;
    SpGemmDesc& operator=(const SpGemmDesc&)=delete;
    ~SpGemmDesc() { cusparseSpGEMM_destroyDescr(descr); }
    operator cusparseSpGEMMDescr_t() const { return descr; }
private:
    cusparseSpGEMMDescr_t descr=nullptr;
};

class GpuInput
{
public:
    explicit GpuInput(const mxArray *array) : ptr(mxGPUCreateFromMxArray(array)) {}
    GpuInput(const GpuInput&)=delete;
    GpuInput& operator=(const GpuInput&)=delete;
    ~GpuInput() { mxGPUDestroyGPUArray(ptr); }
    const mxGPUArray* get() const { return ptr; }
private:
    const mxGPUArray *ptr;
};

class GpuOutput
{
public:
    GpuOutput(int64_t n,mxClassID class_id,mxComplexity complexity)
    {
        const mwSize dims[2]={static_cast<mwSize>(n),1};
        ptr=mxGPUCreateGPUArray(2,dims,class_id,complexity,MX_GPU_DO_NOT_INITIALIZE);
    }
    GpuOutput(const GpuOutput&)=delete;
    GpuOutput& operator=(const GpuOutput&)=delete;
    ~GpuOutput() { mxGPUDestroyGPUArray(ptr); }
    template<typename T> T* data() const { return static_cast<T*>(mxGPUGetData(ptr)); }
    mxArray* to_matlab() const { return mxGPUCreateMxArrayOnGPU(ptr); }
private:
    mxGPUArray *ptr=nullptr;
};

// Sparse gpuArray storage as found in MATLAB memory; GPU arithmetic
// can leave explicit zeros, so stored entries may outnumber nonzeros
struct MatlabCsr
{
    int64_t rows,cols,nnz;
    const int32_t *offsets,*indices;
    const void *values;
    bool is_complex;
};

// CSR operand in cuSPARSE index type I
template<typename I> struct Operand
{
    int64_t rows,cols,nnz;
    const I *offsets,*indices;
    const char *values;
};

// CSR product of a row block of A with B
template<typename I> struct Block
{
    int64_t first_row,rows,nnz;
    DeviceBuffer offsets,indices,values;
};

static unsigned int grid_size(int64_t n)
{
    const int64_t blocks=(n+255)/256;
    return static_cast<unsigned int>(std::max<int64_t>(1,std::min<int64_t>(blocks,std::numeric_limits<int>::max())));
}

// Shifted copy of an index array, optionally into a wider type
template<typename In,typename Out>
__global__ void rebase_kernel(const In *in,Out *out,int64_t n,int64_t origin)
{
    for (int64_t k=blockIdx.x*(int64_t)blockDim.x+threadIdx.x;k<n;k+=(int64_t)gridDim.x*blockDim.x)
        out[k]=static_cast<Out>(in[k]-origin);
}

// Real values as complex values with zero imaginary parts
__global__ void promote_kernel(const double *in,double2 *out,int64_t n)
{
    for (int64_t k=blockIdx.x*(int64_t)blockDim.x+threadIdx.x;k<n;k+=(int64_t)gridDim.x*blockDim.x)
        out[k]=make_double2(in[k],0.0);
}

// One-based row and column indices of each CSR nonzero of a row block
template<typename I>
__global__ void triplet_kernel(const I *offsets,const I *indices,int64_t rows,int64_t nnz,
                               int64_t first_row,int32_t *row_c,int32_t *col_c)
{
    for (int64_t k=blockIdx.x*(int64_t)blockDim.x+threadIdx.x;k<nnz;k+=(int64_t)gridDim.x*blockDim.x)
    {
        int64_t lo=0,hi=rows;
        while (hi-lo>1)
        {
            const int64_t mid=(lo+hi)/2;
            if (offsets[mid]<=k) lo=mid; else hi=mid;
        }
        row_c[k]=static_cast<int32_t>(first_row+lo+1);
        col_c[k]=static_cast<int32_t>(indices[k]+1);
    }
}

template<typename In,typename Out>
static void rebase(const In *in,Out *out,int64_t n,int64_t origin)
{
    rebase_kernel<In,Out><<<grid_size(n),256>>>(in,out,n,origin);
    check_cuda(cudaGetLastError(),"rebase_kernel");
}

static void promote(const void *in,DeviceBuffer &out,int64_t n)
{
    out=DeviceBuffer(n*sizeof(double2));
    promote_kernel<<<grid_size(n),256>>>(static_cast<const double*>(in),out.as<double2>(),n);
    check_cuda(cudaGetLastError(),"promote_kernel");
}

// Host memory test that cannot fault: the operating system copies the bytes into a pipe
static bool host_readable(const void *ptr,size_t n_bytes)
{
    if ((ptr==nullptr)||(reinterpret_cast<uintptr_t>(ptr)%alignof(void*)!=0)) return false;
    int fd[2];
    if (pipe(fd)!=0) return false;
    const bool readable=(write(fd[1],ptr,n_bytes)==static_cast<ssize_t>(n_bytes));
    close(fd[0]); close(fd[1]);
    return readable;
}

// Pointer must be a device allocation on this GPU spanning n_bytes
static void check_extent(const void *ptr,size_t n_bytes,int device,const char *name)
{
    cudaPointerAttributes attributes;
    if ((cudaPointerGetAttributes(&attributes,ptr)!=cudaSuccess)||
        (attributes.type!=cudaMemoryTypeDevice)||(attributes.device!=device))
    {
        cudaGetLastError();
        throw LayoutMismatch(std::string(name)+" pointer is not device memory on the current GPU.");
    }
    CUdeviceptr base=0; size_t size=0;
    const CUdeviceptr start=reinterpret_cast<CUdeviceptr>(ptr);
    if ((cuMemGetAddressRange(&base,&size,start)!=CUDA_SUCCESS)||(start+n_bytes>base+size))
        throw LayoutMismatch(std::string(name)+" allocation is smaller than the matrix requires.");
}

// Stored entry count from nzmax() in MATLAB
static int64_t matlab_nzmax(const mxArray *array)
{
    mxArray *input=const_cast<mxArray*>(array),*count=nullptr;
    mexCallMATLAB(1,&count,1,&input,"nzmax");
    if (mxIsGPUArray(count))
    {
        mxArray *host=nullptr;
        mexCallMATLAB(1,&host,1,&count,"gather");
        mxDestroyArray(count); count=host;
    }
    const int64_t n=static_cast<int64_t>(mxGetScalar(count));
    mxDestroyArray(count);
    return n;
}

static MatlabCsr read_storage(const mxArray *array,const GpuInput &gpu,int device,const char *name)
{
    const mwSize *dims=mxGPUGetDimensions(gpu.get());
    MatlabCsr csr={static_cast<int64_t>(dims[0]),static_cast<int64_t>(dims[1]),matlab_nzmax(array),
                   nullptr,nullptr,nullptr,mxGPUGetComplexity(gpu.get())==mxCOMPLEX};
    mxFree(const_cast<mwSize*>(dims));
    if (csr.nnz==0) return csr;

    // Follow the two pointer hops to the storage block
    const void *handle=gpu.get();
    if (!host_readable(handle,sizeof(void*))) throw LayoutMismatch(std::string(name)+" handle is unreadable.");
    const void *object=*static_cast<void* const*>(handle);
    if (!host_readable(object,sizeof(void*))) throw LayoutMismatch(std::string(name)+" object is unreadable.");
    const char *storage=*static_cast<char* const*>(object);
    if (!host_readable(storage,168)) throw LayoutMismatch(std::string(name)+" storage block is unreadable.");
    std::memcpy(&csr.values,storage+144,sizeof(void*));
    std::memcpy(&csr.indices,storage+152,sizeof(void*));
    std::memcpy(&csr.offsets,storage+160,sizeof(void*));

    // Validate the buffers against the matrix MATLAB reports
    const size_t value_size=csr.is_complex?sizeof(double2):sizeof(double);
    check_extent(csr.values,csr.nnz*value_size,device,name);
    check_extent(csr.indices,csr.nnz*sizeof(int32_t),device,name);
    check_extent(csr.offsets,(csr.rows+1)*sizeof(int32_t),device,name);
    int32_t first=-1,last=-1;
    check_cuda(cudaMemcpy(&first,csr.offsets,sizeof(int32_t),cudaMemcpyDeviceToHost),"cudaMemcpy");
    check_cuda(cudaMemcpy(&last,csr.offsets+csr.rows,sizeof(int32_t),cudaMemcpyDeviceToHost),"cudaMemcpy");
    if ((first!=0)||(last!=csr.nnz))
        throw LayoutMismatch(std::string(name)+" row offsets disagree with MATLAB nzmax.");
    return csr;
}

// MATLAB storage in place for 32-bit indices, widened copies for 64-bit
template<typename I>
static Operand<I> make_operand(const MatlabCsr &csr,const void *values,
                               DeviceBuffer &offsets,DeviceBuffer &indices)
{
    Operand<I> op={csr.rows,csr.cols,csr.nnz,nullptr,nullptr,static_cast<const char*>(values)};
    if constexpr (sizeof(I)==sizeof(int32_t))
    {
        op.offsets=csr.offsets; op.indices=csr.indices;
    }
    else
    {
        offsets=DeviceBuffer((csr.rows+1)*sizeof(I)); indices=DeviceBuffer(csr.nnz*sizeof(I));
        rebase(csr.offsets,offsets.as<I>(),csr.rows+1,0);
        rebase(csr.indices,indices.as<I>(),csr.nnz,0);
        op.offsets=offsets.as<I>(); op.indices=indices.as<I>();
    }
    return op;
}

template<typename I>
static Block<I> spgemm_block(cusparseHandle_t handle,const Operand<I> &a,const std::vector<int64_t> &a_offsets,
                             const Operand<I> &b,int64_t first_row,int64_t end_row,
                             cusparseSpGEMMAlg_t alg,cudaDataType value_type,size_t value_size)
{
    const cusparseIndexType_t index_type=(sizeof(I)==sizeof(int32_t))?CUSPARSE_INDEX_32I:CUSPARSE_INDEX_64I;
    const cusparseOperation_t op=CUSPARSE_OPERATION_NON_TRANSPOSE;
    const double real_one=1.0,real_zero=0.0;
    const double2 complex_one=make_double2(1.0,0.0),complex_zero=make_double2(0.0,0.0);
    const void *alpha=(value_type==CUDA_R_64F)?static_cast<const void*>(&real_one):&complex_one;
    const void *beta=(value_type==CUDA_R_64F)?static_cast<const void*>(&real_zero):&complex_zero;

    // Row block of A: offsets into MATLAB arrays, rebased row offsets
    const int64_t rows=end_row-first_row,origin=a_offsets[first_row];
    DeviceBuffer block_offsets; const I *offsets=a.offsets;
    if (rows<a.rows)
    {
        block_offsets=DeviceBuffer((rows+1)*sizeof(I));
        rebase(a.offsets+first_row,block_offsets.as<I>(),rows+1,origin);
        offsets=block_offsets.as<I>();
    }
    const SpMat mat_a(rows,a.cols,a_offsets[end_row]-origin,offsets,a.indices+origin,
                      a.values+origin*value_size,index_type,value_type);
    const SpMat mat_b(b.rows,b.cols,b.nnz,b.offsets,b.indices,b.values,index_type,value_type);
    Block<I> c={first_row,rows,0,DeviceBuffer((rows+1)*sizeof(I)),DeviceBuffer(),DeviceBuffer()};
    const SpMat mat_c(rows,b.cols,0,c.offsets.template as<I>(),nullptr,nullptr,index_type,value_type);
    const SpGemmDesc desc;

    // Work estimation
    size_t size_1=0;
    check_cusparse(cusparseSpGEMM_workEstimation(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                                 alg,desc,&size_1,nullptr),"cusparseSpGEMM_workEstimation");
    const DeviceBuffer buffer_1(size_1);
    check_cusparse(cusparseSpGEMM_workEstimation(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                                 alg,desc,&size_1,buffer_1.as<void>()),
                   "cusparseSpGEMM_workEstimation");

    // Compute buffer size; ALG3 chunk fraction is halved until both
    // its estimation and its compute buffers fit into free memory
    size_t size_2=0;
    if (alg==CUSPARSE_SPGEMM_ALG1)
    {
        check_cusparse(cusparseSpGEMM_compute(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                              alg,desc,&size_2,nullptr),"cusparseSpGEMM_compute");
    }
    else
    {
        const bool chunked=(alg==CUSPARSE_SPGEMM_ALG3);
        int64_t n_products=0;
        if (chunked) check_cusparse(cusparseSpGEMM_getNumProducts(desc,&n_products),"cusparseSpGEMM_getNumProducts");
        for (float chunk_fraction=1.0f;;chunk_fraction*=0.5f)
        {
            if (chunked&&(chunk_fraction<1)&&(chunk_fraction*n_products<1))
                throw OutOfResources("ALG3 chunks do not fit into free device memory.");
            size_t size_3=0,free_bytes=0,total_bytes=0;
            check_cusparse(cusparseSpGEMM_estimateMemory(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                                         alg,desc,chunk_fraction,&size_3,nullptr,nullptr),
                           "cusparseSpGEMM_estimateMemory");
            check_cuda(cudaMemGetInfo(&free_bytes,&total_bytes),"cudaMemGetInfo");
            if (chunked&&(size_3>free_bytes)) continue;
            const DeviceBuffer buffer_3(size_3);
            check_cusparse(cusparseSpGEMM_estimateMemory(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                                         alg,desc,chunk_fraction,&size_3,buffer_3.as<void>(),
                                                         &size_2),"cusparseSpGEMM_estimateMemory");
            if ((!chunked)||(size_2<=free_bytes)) break;
        }
    }
    const DeviceBuffer buffer_2(size_2);
    check_cusparse(cusparseSpGEMM_compute(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,
                                          alg,desc,&size_2,buffer_2.as<void>()),"cusparseSpGEMM_compute");

    // Copy the product into its own CSR arrays
    int64_t c_rows=0,c_cols=0;
    check_cusparse(cusparseSpMatGetSize(mat_c,&c_rows,&c_cols,&c.nnz),"cusparseSpMatGetSize");
    c.indices=DeviceBuffer(c.nnz*sizeof(I)); c.values=DeviceBuffer(c.nnz*value_size);
    check_cusparse(cusparseCsrSetPointers(mat_c,c.offsets.template as<void>(),c.indices.template as<void>(),
                                          c.values.template as<void>()),"cusparseCsrSetPointers");
    check_cusparse(cusparseSpGEMM_copy(handle,op,op,alpha,mat_a,mat_b,beta,mat_c,value_type,alg,desc),
                   "cusparseSpGEMM_copy");
    return c;
}

// Multiply rows [first_row,end_row) of A by B, bisecting on resource exhaustion
template<typename I>
static void multiply_rows(cusparseHandle_t handle,const Operand<I> &a,const std::vector<int64_t> &a_offsets,
                          const Operand<I> &b,int64_t first_row,int64_t end_row,cusparseSpGEMMAlg_t alg,
                          cudaDataType value_type,size_t value_size,std::vector<Block<I>> &blocks)
{
    try
    {
        blocks.push_back(spgemm_block(handle,a,a_offsets,b,first_row,end_row,alg,value_type,value_size));
        return;
    }
    catch (const OutOfResources&)
    {
        if (end_row-first_row<2) throw;
    }

    // Refuse to continue after an asynchronous fault
    check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");

    // Split where the block has half of its nonzeros
    const int64_t half=(a_offsets[first_row]+a_offsets[end_row])/2;
    int64_t mid=std::upper_bound(a_offsets.begin()+first_row,a_offsets.begin()+end_row,half)-a_offsets.begin()-1;
    mid=std::min(std::max(mid,first_row+1),end_row-1);
    multiply_rows(handle,a,a_offsets,b,first_row,mid,alg,value_type,value_size,blocks);
    multiply_rows(handle,a,a_offsets,b,mid,end_row,alg,value_type,value_size,blocks);
}

template<typename I>
static void run(cusparseHandle_t handle,const MatlabCsr &csr_a,const void *values_a,
                const MatlabCsr &csr_b,const void *values_b,cusparseSpGEMMAlg_t alg,
                bool is_complex,mxArray *plhs[])
{
    const cudaDataType value_type=is_complex?CUDA_C_64F:CUDA_R_64F;
    const size_t value_size=is_complex?sizeof(double2):sizeof(double);

    // Operands, sharing index arrays when B is A
    DeviceBuffer wide[4];
    const Operand<I> a=make_operand<I>(csr_a,values_a,wide[0],wide[1]);
    const Operand<I> b=((csr_b.offsets==csr_a.offsets)&&(values_b==values_a))?a:
                       make_operand<I>(csr_b,values_b,wide[2],wide[3]);

    // Host copy of the row offsets of A for block boundaries
    std::vector<int32_t> offsets_32(csr_a.rows+1);
    check_cuda(cudaMemcpy(offsets_32.data(),csr_a.offsets,offsets_32.size()*sizeof(int32_t),
                          cudaMemcpyDeviceToHost),"cudaMemcpy");
    const std::vector<int64_t> a_offsets(offsets_32.begin(),offsets_32.end());

    // Row blocks of the product
    std::vector<Block<I>> blocks;
    multiply_rows(handle,a,a_offsets,b,0,a.rows,alg,value_type,value_size,blocks);
    int64_t nnz=0;
    for (const Block<I> &block : blocks) nnz+=block.nnz;
    if (nnz>std::numeric_limits<int32_t>::max())
        throw std::runtime_error("The product has more nonzeros than a MATLAB sparse gpuArray can hold.");

    // One-based triplets in row-major order
    const GpuOutput row_c(nnz,mxINT32_CLASS,mxREAL),col_c(nnz,mxINT32_CLASS,mxREAL);
    const GpuOutput val_c(nnz,mxDOUBLE_CLASS,is_complex?mxCOMPLEX:mxREAL);
    int64_t base=0;
    for (const Block<I> &block : blocks)
    {
        if (block.nnz==0) continue;
        triplet_kernel<I><<<grid_size(block.nnz),256>>>(block.offsets.template as<I>(),
                                                         block.indices.template as<I>(),
                                                         block.rows,block.nnz,block.first_row,
                                                         row_c.data<int32_t>()+base,col_c.data<int32_t>()+base);
        check_cuda(cudaGetLastError(),"triplet_kernel");
        check_cuda(cudaMemcpy(val_c.data<char>()+base*value_size,block.values.template as<void>(),
                              block.nnz*value_size,cudaMemcpyDeviceToDevice),"cudaMemcpy");
        base+=block.nnz;
    }
    check_cuda(cudaDeviceSynchronize(),"cudaDeviceSynchronize");
    plhs[0]=row_c.to_matlab(); plhs[1]=col_c.to_matlab(); plhs[2]=val_c.to_matlab();
}

static void grumble(int nlhs,int nrhs,const mxArray *prhs[])
{
    if (nrhs!=3) throw std::invalid_argument("Three inputs are required.");
    if (nlhs!=3) throw std::invalid_argument("Three outputs are required.");
    if ((!mxIsGPUArray(prhs[0]))||(!mxIsGPUArray(prhs[1])))
        throw std::invalid_argument("A and B must be gpuArrays.");
    if ((!mxIsDouble(prhs[2]))||mxIsComplex(prhs[2])||(mxGetNumberOfElements(prhs[2])!=1))
        throw std::invalid_argument("alg must be a real double scalar.");
    const double alg=mxGetScalar(prhs[2]);
    if ((alg!=1)&&(alg!=2)&&(alg!=3)) throw std::invalid_argument("alg must be 1, 2, or 3.");
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    static char message[1024];
    const char *identifier=nullptr;
    try
    {
        if (mxInitGPU()!=MX_GPU_SUCCESS) throw std::runtime_error("Failed to initialise the MATLAB GPU API.");
        grumble(nlhs,nrhs,prhs);
        const GpuInput gpu_a(prhs[0]),gpu_b(prhs[1]);
        if ((!mxGPUIsSparse(gpu_a.get()))||(!mxGPUIsSparse(gpu_b.get()))||
            (mxGPUGetClassID(gpu_a.get())!=mxDOUBLE_CLASS)||(mxGPUGetClassID(gpu_b.get())!=mxDOUBLE_CLASS))
            throw std::invalid_argument("A and B must be sparse double gpuArrays.");
        int device=0;
        check_cuda(cudaGetDevice(&device),"cudaGetDevice");
        const MatlabCsr a=read_storage(prhs[0],gpu_a,device,"A");
        const MatlabCsr b=read_storage(prhs[1],gpu_b,device,"B");
        if (a.cols!=b.rows) throw std::invalid_argument("A and B dimensions are inconsistent.");
        const bool is_complex=a.is_complex||b.is_complex;

        // Operands without stored entries give a real empty product, as in native mtimes
        if ((a.nnz==0)||(b.nnz==0))
        {
            const GpuOutput row_c(0,mxINT32_CLASS,mxREAL),col_c(0,mxINT32_CLASS,mxREAL);
            const GpuOutput val_c(0,mxDOUBLE_CLASS,mxREAL);
            plhs[0]=row_c.to_matlab(); plhs[1]=col_c.to_matlab(); plhs[2]=val_c.to_matlab();
            return;
        }

        // Complex copies of the values of a real operand in a mixed product
        DeviceBuffer complex_a,complex_b;
        const void *values_a=a.values,*values_b=b.values;
        if (is_complex&&!a.is_complex) { promote(a.values,complex_a,a.nnz); values_a=complex_a.as<void>(); }
        if (is_complex&&!b.is_complex) { promote(b.values,complex_b,b.nnz); values_b=complex_b.as<void>(); }

        // ALG3 needs 64-bit indices, ALG1 and ALG2 use MATLAB's own
        const CusparseHandle handle;
        const int alg=static_cast<int>(mxGetScalar(prhs[2]));
        if (alg==3)
            run<int64_t>(handle,a,values_a,b,values_b,CUSPARSE_SPGEMM_ALG3,is_complex,plhs);
        else
            run<int32_t>(handle,a,values_a,b,values_b,(alg==1)?CUSPARSE_SPGEMM_ALG1:CUSPARSE_SPGEMM_ALG2,
                         is_complex,plhs);
    }
    catch (const LayoutMismatch &err)
    {
        identifier="Spinach:cuda_sparse_by_sparse_mex:layout";
        std::snprintf(message,sizeof(message),"%s",err.what());
    }
    catch (const std::invalid_argument &err)
    {
        identifier="Spinach:cuda_sparse_by_sparse_mex:input";
        std::snprintf(message,sizeof(message),"%s",err.what());
    }
    catch (const std::exception &err)
    {
        identifier="Spinach:cuda_sparse_by_sparse_mex:runtime";
        std::snprintf(message,sizeof(message),"%s",err.what());
    }
    if (identifier!=nullptr) mexErrMsgIdAndTxt(identifier,"%s",message);
}
