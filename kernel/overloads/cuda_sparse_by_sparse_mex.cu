/* Low-level CUDA CSR-by-CSR product for MATLAB R2026b.
 * C=cuda_sparse_by_sparse_mex(A,B), sparse double gpuArrays, real or complex.
 * Gustavson row-wise symbolic/numeric phases; no cuSPARSE or BLAS calls.
 * Shared-memory dense accumulators for narrow column ranges, bounded hash
 * accumulators otherwise. Overflowed symbolic tasks split their column range.
 * Inputs are read-only MATLAB CSR buffers. A fresh MATLAB-owned CSR structure
 * is allocated once from the symbolic pattern; numeric values are written
 * directly to it, and the ready object is returned without a numeric copy.
 * The allocation currently uses MATLAB sparse() with placeholder values.
 * Layout offsets are undocumented and checked before use. Input/output storage
 * uses MATLAB's int32 ABI; all offsets, task sizes, and host sums use int64.
 * Sources: Gustavson (1978), DOI 10.1145/355791.355796; Davis et al., Sparse
 * Direct Methods, Algorithm 2.2; GraphBLAS GB_AxB_saxpy3.c; spECK HashSpGEMM.
 */
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <climits>
#include <stdexcept>
#include <string>
#include <vector>

struct LayoutMismatch:std::runtime_error { using std::runtime_error::runtime_error; };
static void check_cuda(cudaError_t code,const char *stage)
{
    if (code!=cudaSuccess) throw std::runtime_error(std::string(stage)+": "+cudaGetErrorString(code));
}
struct Buffer
{
    void *ptr=nullptr;
    explicit Buffer(size_t bytes) { if (bytes) check_cuda(cudaMalloc(&ptr,bytes),"cudaMalloc"); }
    ~Buffer() { if (ptr) cudaFree(ptr); }
    Buffer(const Buffer&)=delete;
    Buffer& operator=(const Buffer&)=delete;
    template<class T> T* as() const { return static_cast<T*>(ptr); }
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


struct MatlabCsr
{
    int64_t rows,cols,nnz;
    const int32_t *offsets,*indices;
    const void *values;
    bool is_complex;
};

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

// Allocated entry capacity from nzmax() in MATLAB
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
    if ((first!=0)||(last<0)||(last>csr.nnz))
        throw LayoutMismatch(std::string(name)+" row offsets exceed MATLAB nzmax capacity.");
    csr.nnz=last;
    return csr;
}


// Symbolic tasks partition output rows into disjoint column ranges
struct Task { int64_t row,lo,hi,base; int32_t count,overflow; };
__host__ __device__ int table_capacity(int64_t bound,int limit)
{
    int slots=1;
    while ((slots<limit)&&(slots<2*bound)) slots*=2;
    return slots;
}
__device__ int64_t lower_bound_col(const int32_t *cols,int64_t first,int64_t last,int64_t key)
{
    while (first<last) {
        const int64_t mid=first+(last-first)/2;
        if (cols[mid]<key) first=mid+1; else last=mid;
    }
    return first;
}
__device__ int hash_slot(int col,int capacity)
{
    return (static_cast<uint32_t>(col)*2654435761u)&(capacity-1);
}
__device__ int insert_key(int *keys,int capacity,int col)
{
    int slot=hash_slot(col,capacity);
    for (int k=0;k<capacity;k++) {
        const int old=atomicCAS(keys+slot,-1,col);
        if ((old==-1)||(old==col)) return slot;
        slot=(slot+1)&(capacity-1);
    }
    return -1;
}
__device__ void sort_keys(int *keys,int capacity)
{
    for (int k=threadIdx.x;k<capacity;k+=blockDim.x)
        if (keys[k]<0) keys[k]=INT_MAX;
    __syncthreads();
    for (int width=2;width<=capacity;width*=2)
        for (int gap=width/2;gap;gap/=2) {
            for (int k=threadIdx.x;k<capacity;k+=blockDim.x) {
                const int other=k^gap;
                if (other>k) {
                    const int a=keys[k],b=keys[other];
                    if (((k&width)==0)?(a>b):(a<b)) { keys[k]=b; keys[other]=a; }
                }
            }
            __syncthreads();
        }
}
// Validate every CSR index and row boundary before launching arithmetic
__global__ void validate_csr(MatlabCsr a,int *invalid)
{
    for (int64_t row=blockIdx.x*(int64_t)blockDim.x+threadIdx.x;row<a.rows;row+=(int64_t)blockDim.x*gridDim.x) {
        const int64_t first=a.offsets[row],last=a.offsets[row+1];
        if ((first<0)||(last<first)||(last>a.nnz)) { atomicExch(invalid,1); continue; }
        for (int64_t p=first;p<last;p++)
            if ((a.indices[p]<0)||(a.indices[p]>=a.cols)||((p>first)&&(a.indices[p]<=a.indices[p-1])))
                atomicExch(invalid,1);
    }
}
// Count or emit each task's structural columns using a bitset or a hash table
__global__ void symbolic(MatlabCsr a,MatlabCsr b,Task *tasks,int n_tasks,int dense_width,int capacity,
                         int32_t *rows_c,int32_t *cols_c,bool emit)
{
    extern __shared__ int keys[];
    __shared__ int count,overflow,local_capacity,prefix[256];
    const int tid=threadIdx.x,lane=tid&31,warp=tid/32;
    for (int task_id=blockIdx.x;task_id<n_tasks;task_id+=gridDim.x) {
        Task &task=tasks[task_id];
        const bool dense=task.hi-task.lo<=dense_width;
        const int width=static_cast<int>(task.hi-task.lo);
        if (tid==0) {
            int64_t bound=emit?task.count:0;
            if (!emit&&!dense)
                for (int64_t p=a.offsets[task.row];(p<a.offsets[task.row+1])&&(bound<capacity/2);p++)
                    bound+=static_cast<int64_t>(b.offsets[a.indices[p]+1])-b.offsets[a.indices[p]];
            local_capacity=table_capacity(bound,capacity);
        }
        __syncthreads();
        const int slots=dense?(width+31)/32:local_capacity;
        for (int k=tid;k<slots;k+=blockDim.x) keys[k]=dense?0:-1;
        if (tid==0) { count=0; overflow=0; }
        __syncthreads();
        const int64_t first_a=a.offsets[task.row],last_a=a.offsets[task.row+1];
        for (int64_t p=first_a+warp;p<last_a;p+=blockDim.x/32) {
            const int64_t row_b=a.indices[p];
            const int64_t first_b=b.offsets[row_b],last_b=b.offsets[row_b+1];
            const int64_t begin=(task.lo==0)?first_b:lower_bound_col(b.indices,first_b,last_b,task.lo);
            const int64_t end=(task.hi==b.cols)?last_b:lower_bound_col(b.indices,begin,last_b,task.hi);
            for (int64_t q=begin+lane;q<end;q+=32) {
                const int col=b.indices[q];
                if (dense) atomicOr(reinterpret_cast<unsigned int*>(keys)+(col-task.lo)/32,1u<<((col-task.lo)%32));
                else if ((atomicAdd(&overflow,0)==0)&&(insert_key(keys,local_capacity,col)<0)) atomicExch(&overflow,1);
            }
        }
        __syncthreads();
        int local=0;
        for (int k=tid;k<slots;k+=blockDim.x) local+=dense?__popc(keys[k]):(keys[k]!=-1);
        atomicAdd(&count,local);
        __syncthreads();
        if (!emit) {
            if (tid==0) { task.count=count; task.overflow=overflow||(!dense&&(count>local_capacity/2)); }
        } else if (dense) {
            int output=0;
            for (int start=0;start<width;start+=blockDim.x) {
                const int k=start+tid;
                const int present=(k<width)?((static_cast<uint32_t>(keys[k/32])>>(k%32))&1u):0;
                prefix[tid]=present;
                __syncthreads();
                for (int gap=1;gap<blockDim.x;gap*=2) {
                    const int prev=(tid>=gap)?prefix[tid-gap]:0;
                    __syncthreads();
                    prefix[tid]+=prev;
                    __syncthreads();
                }
                if (present) {
                    const int64_t pos=task.base+output+prefix[tid]-1;
                    rows_c[pos]=static_cast<int32_t>(task.row+1);
                    cols_c[pos]=static_cast<int32_t>(task.lo+k+1);
                }
                output+=prefix[blockDim.x-1];
                __syncthreads();
            }
        } else {
            sort_keys(keys,local_capacity);
            for (int k=tid;k<task.count;k+=blockDim.x) {
                rows_c[task.base+k]=static_cast<int32_t>(task.row+1);
                cols_c[task.base+k]=keys[k]+1;
            }
        }
        __syncthreads();
    }
}
__device__ double load_value(const double *values,int64_t k) { return values[k]; }
__device__ double2 load_value(const double2 *values,int64_t k) { return values[k]; }
__device__ double multiply(double a,double b) { return a*b; }
__device__ double2 multiply(double a,double2 b) { return make_double2(a*b.x,a*b.y); }
__device__ double2 multiply(double2 a,double b) { return make_double2(a.x*b,a.y*b); }
__device__ double2 multiply(double2 a,double2 b)
{
    return make_double2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}
__device__ double zero_value(double) { return 0.0; }
__device__ double2 zero_value(double2) { return make_double2(0.0,0.0); }
__device__ double sum_value(double a,double b) { return a+b; }
__device__ double2 sum_value(double2 a,double2 b) { return make_double2(a.x+b.x,a.y+b.y); }
__device__ void add_value(double *dest,double x) { atomicAdd(dest,x); }
__device__ void add_value(double2 *dest,double2 x) { atomicAdd(&dest->x,x.x); atomicAdd(&dest->y,x.y); }
// Each block accumulates one row/range; all matrix values stay on the GPU
// Dense ranges index shared accumulators directly; wider ranges use bounded hashing
// Shared-memory storage is reused across the task fleet, never per scalar product
// Template arguments avoid promoting real inputs to complex buffers
// Row-major sorted symbolic offsets address the newly allocated MATLAB object
// No input storage is written
// (These explanatory details are part of the algorithm, not runtime options.)
template<class AValue,class BValue,class CValue>
__global__ void numeric(MatlabCsr a,MatlabCsr b,const Task *tasks,int n_tasks,int dense_width,int max_capacity,int acc_stride,
                        const int32_t *cols_c,CValue *values_c)
{
    extern __shared__ double2 shared[];
    CValue *acc=reinterpret_cast<CValue*>(shared);
    int *keys=reinterpret_cast<int*>(acc+acc_stride);
    const int tid=threadIdx.x,lane=tid&31,warp=tid/32;
    for (int task_id=blockIdx.x;task_id<n_tasks;task_id+=gridDim.x) {
        const Task task=tasks[task_id];
        const bool dense=task.hi-task.lo<=dense_width;
        const int capacity=table_capacity(task.count,max_capacity);
        const int width=dense?static_cast<int>(task.hi-task.lo):capacity;
        for (int k=tid;k<width;k+=blockDim.x) { acc[k]=zero_value(CValue{}); if (!dense) keys[k]=-1; }
        __syncthreads();
        if (!dense)
            for (int k=tid;k<task.count;k+=blockDim.x) insert_key(keys,capacity,cols_c[task.base+k]);
        __syncthreads();
        const int64_t first_a=a.offsets[task.row],last_a=a.offsets[task.row+1];
        const int n_warps=blockDim.x/32;
        const int64_t col_lo=dense?task.lo+(task.hi-task.lo)*warp/n_warps:task.lo;
        const int64_t col_hi=dense?task.lo+(task.hi-task.lo)*(warp+1)/n_warps:task.hi;
        for (int64_t p=first_a+(dense?0:warp);p<last_a;p+=dense?1:n_warps) {
            const AValue av=load_value(static_cast<const AValue*>(a.values),p);
            const int64_t row_b=a.indices[p];
            const int64_t first_b=b.offsets[row_b],last_b=b.offsets[row_b+1];
            const int64_t begin=(col_lo==0)?first_b:lower_bound_col(b.indices,first_b,last_b,col_lo);
            const int64_t end=(col_hi==b.cols)?last_b:lower_bound_col(b.indices,begin,last_b,col_hi);
            for (int64_t q=begin+lane;q<end;q+=32) {
                const int col=b.indices[q];
                int slot=static_cast<int>(col-task.lo);
                const CValue value=multiply(av,load_value(static_cast<const BValue*>(b.values),q));
                if (dense) acc[slot]=sum_value(acc[slot],value);
                else {
                    slot=hash_slot(col,capacity);
                    while (keys[slot]!=col) slot=(slot+1)&(capacity-1);
                    add_value(acc+slot,value);
                }
            }
            __syncwarp();
        }
        __syncthreads();
        for (int k=tid;k<task.count;k+=blockDim.x) {
            const int col=cols_c[task.base+k];
            int slot=static_cast<int>(col-task.lo);
            if (!dense) { slot=hash_slot(col,capacity); while (keys[slot]!=col) slot=(slot+1)&(capacity-1); }
            values_c[task.base+k]=acc[slot];
        }
        __syncthreads();
    }
}
template<class T> __global__ void fill_ones(T *values,int64_t n)
{
    for (int64_t k=blockIdx.x*(int64_t)blockDim.x+threadIdx.x;k<n;k+=(int64_t)gridDim.x*blockDim.x) {
        if constexpr (sizeof(T)==sizeof(double)) values[k]=1.0;
        else values[k]=make_double2(1.0,1.0);
    }
}
static mxArray *allocate_pattern(const GpuOutput &rows,const GpuOutput &cols,int64_t nnz,
                                int64_t m,int64_t n,bool complex)
{
    const GpuOutput ones(nnz,mxDOUBLE_CLASS,complex?mxCOMPLEX:mxREAL);
    if (complex) fill_ones<<<256,256>>>(ones.data<double2>(),nnz);
    else fill_ones<<<256,256>>>(ones.data<double>(),nnz);
    check_cuda(cudaDeviceSynchronize(),"fill_ones");
    mxArray *args[5]={rows.to_matlab(),cols.to_matlab(),ones.to_matlab(),
                     mxCreateDoubleScalar(static_cast<double>(m)),mxCreateDoubleScalar(static_cast<double>(n))};
    mxArray *result=nullptr;
    mxArray *error=mexCallMATLABWithTrap(1,&result,5,args,"sparse");
    for (auto *arg:args) mxDestroyArray(arg);
    if (error) { mxDestroyArray(error); throw std::runtime_error("MATLAB sparse pattern allocation failed."); }
    return result;
}
static void product(const mxArray *array_a,const mxArray *array_b,mxArray **result)
{
    const GpuInput gpu_a(array_a),gpu_b(array_b);
    if ((!mxGPUIsSparse(gpu_a.get()))||(!mxGPUIsSparse(gpu_b.get()))||
        (mxGPUGetClassID(gpu_a.get())!=mxDOUBLE_CLASS)||(mxGPUGetClassID(gpu_b.get())!=mxDOUBLE_CLASS))
        throw std::invalid_argument("Inputs must be sparse double gpuArrays.");
    int device; check_cuda(cudaGetDevice(&device),"cudaGetDevice");
    const MatlabCsr a=read_storage(array_a,gpu_a,device,"A"),b=read_storage(array_b,gpu_b,device,"B");
    if (a.cols!=b.rows) throw std::invalid_argument("Inner matrix dimensions must agree.");
    if ((a.rows>INT_MAX)||(b.cols>INT_MAX)) throw std::invalid_argument("Dimensions exceed MATLAB's sparse GPU index ABI.");
    const bool complex=a.is_complex||b.is_complex;
    if ((!a.nnz)||(!b.nnz)||(!a.rows)||(!b.cols)) {
        const GpuOutput rows(0,mxINT32_CLASS,mxREAL),cols(0,mxINT32_CLASS,mxREAL);
        *result=allocate_pattern(rows,cols,0,a.rows,b.cols,false); return;
    }
    const Buffer invalid(sizeof(int));
    check_cuda(cudaMemset(invalid.ptr,0,sizeof(int)),"cudaMemset");
    validate_csr<<<256,256>>>(a,invalid.as<int>());
    validate_csr<<<256,256>>>(b,invalid.as<int>());
    int bad=0; check_cuda(cudaMemcpy(&bad,invalid.ptr,sizeof(int),cudaMemcpyDeviceToHost),"validate_csr");
    if (bad) throw LayoutMismatch("CSR row offsets or sorted column indices are invalid.");
    cudaDeviceProp prop; check_cuda(cudaGetDeviceProperties(&prop,device),"cudaGetDeviceProperties");
    const int value_size=complex?sizeof(double2):sizeof(double);
    const int budget=std::min(prop.sharedMemPerBlockOptin,prop.sharedMemPerMultiprocessor/2)-2048;
    int capacity=1;
    while (2*capacity*(value_size+sizeof(int))<=budget) capacity*=2;
    // Reserve room for the hash keys when deriving the actual dense width
    const int dense_cols=((budget-capacity*sizeof(int))/value_size/256)*256;
    std::vector<Task> tasks;
    tasks.reserve(a.rows);
    for (int64_t row=0;row<a.rows;row++) tasks.push_back({row,0,b.cols,0,0,0});
    const size_t sym_shared=std::max(capacity,(dense_cols+31)/32)*sizeof(int);
    // Symbolic failures split only column ranges, without materialising products
    for (;;) {
        const Buffer device_tasks(tasks.size()*sizeof(Task));
        check_cuda(cudaMemcpy(device_tasks.ptr,tasks.data(),tasks.size()*sizeof(Task),cudaMemcpyHostToDevice),"tasks upload");
        symbolic<<<std::min<size_t>(tasks.size(),prop.multiProcessorCount*8),256,sym_shared>>>
            (a,b,device_tasks.as<Task>(),tasks.size(),dense_cols,capacity,nullptr,nullptr,false);
        check_cuda(cudaMemcpy(tasks.data(),device_tasks.ptr,tasks.size()*sizeof(Task),cudaMemcpyDeviceToHost),"symbolic counts");
        std::vector<Task> next;
        bool split=false;
        for (const Task &task:tasks) {
            if (task.overflow) {
                const int64_t mid=task.lo+(task.hi-task.lo)/2;
                if ((mid==task.lo)||(mid==task.hi)) throw std::runtime_error("Symbolic task cannot be split.");
                next.push_back({task.row,task.lo,mid,0,0,0}); next.push_back({task.row,mid,task.hi,0,0,0}); split=true;
            } else if (task.count) next.push_back(task);
        }
        tasks.swap(next);
        if (!split) break;
    }
    int64_t nnz=0;
    int max_hash=0,acc_stride=0;
    for (Task &task:tasks) {
        task.base=nnz; nnz+=task.count;
        if (task.hi-task.lo<=dense_cols) acc_stride=std::max(acc_stride,static_cast<int>(task.hi-task.lo));
        else max_hash=std::max(max_hash,table_capacity(task.count,capacity));
    }
    acc_stride=std::max(acc_stride,max_hash);
    const size_t shared_bytes=acc_stride*value_size+max_hash*sizeof(int);
    if (nnz>INT_MAX) throw std::runtime_error("Output exceeds MATLAB's sparse GPU index ABI.");
    const GpuOutput rows(nnz,mxINT32_CLASS,mxREAL),cols(nnz,mxINT32_CLASS,mxREAL);
    if (!nnz) { *result=allocate_pattern(rows,cols,0,a.rows,b.cols,complex); return; }
    const Buffer device_tasks(tasks.size()*sizeof(Task));
    check_cuda(cudaMemcpy(device_tasks.ptr,tasks.data(),tasks.size()*sizeof(Task),cudaMemcpyHostToDevice),"tasks upload");
    const int grid=std::min<size_t>(tasks.size(),prop.multiProcessorCount*8);
    symbolic<<<grid,256,sym_shared>>>(a,b,device_tasks.as<Task>(),tasks.size(),dense_cols,capacity,rows.data<int32_t>(),cols.data<int32_t>(),true);
    check_cuda(cudaDeviceSynchronize(),"symbolic pattern");
    *result=allocate_pattern(rows,cols,nnz,a.rows,b.cols,complex);
    const GpuInput gpu_c(*result);
    const MatlabCsr c=read_storage(*result,gpu_c,device,"C");
    if (c.nnz!=nnz) throw LayoutMismatch("Allocated output structure has an unexpected entry count.");
    #define RUN_NUMERIC(AV,BV,CV) \
        check_cuda(cudaFuncSetAttribute(numeric<AV,BV,CV>,cudaFuncAttributeMaxDynamicSharedMemorySize,shared_bytes),"shared memory opt-in"); \
        numeric<AV,BV,CV><<<grid,512,shared_bytes>>>(a,b,device_tasks.as<Task>(),tasks.size(),dense_cols,capacity,acc_stride,c.indices,const_cast<CV*>(static_cast<const CV*>(c.values)))
    if (a.is_complex&&b.is_complex) { RUN_NUMERIC(double2,double2,double2); }
    else if (a.is_complex) { RUN_NUMERIC(double2,double,double2); }
    else if (b.is_complex) { RUN_NUMERIC(double,double2,double2); }
    else { RUN_NUMERIC(double,double,double); }
    #undef RUN_NUMERIC
    check_cuda(cudaDeviceSynchronize(),"numeric product");
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    static char message[1024]; const char *identifier=nullptr;
    try {
        if (nrhs!=2||nlhs!=1) throw std::invalid_argument("Two inputs and one output are required.");
        if (!mxIsGPUArray(prhs[0])||!mxIsGPUArray(prhs[1])) throw std::invalid_argument("Inputs must be gpuArrays.");
        if (mxInitGPU()!=MX_GPU_SUCCESS) throw std::runtime_error("GPU API initialisation failed.");
        product(prhs[0],prhs[1],plhs);
    } catch (const LayoutMismatch &err) {
        identifier="Spinach:cuda_sparse_by_sparse_mex:layout"; std::snprintf(message,sizeof(message),"%s",err.what());
    } catch (const std::exception &err) {
        identifier="Spinach:cuda_sparse_by_sparse_mex:runtime"; std::snprintf(message,sizeof(message),"%s",err.what());
    }
    if (identifier) { if ((nlhs>0)&&plhs[0]) { mxDestroyArray(plhs[0]); plhs[0]=nullptr; } mexErrMsgIdAndTxt(identifier,"%s",message); }
}
