#ifndef CONVENIENCE_CUH_
#define CUDA_CONVENIENCE_CUH_

#include <cuda_runtime_api.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#define STR2(x) #x
#define STRINGIFY(x) STR2(x)

#define FILE_LINE __FILE__ ":" STRINGIFY(__LINE__)

static inline int getGridDim(int x, int y)
{
    return (x + y - 1) / y;
}

/*static inline void cudaCheckError(std::string fileline){
    cudaSafeCall(cudaGetLastError(), fileline);
}*/

#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cudaCheckError() __cudaSafeCall(cudaGetLastError(), __FILE__, __LINE__)
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
                printf("%s(%i) : cudaSafeCall() Runtime API error : %s.\n",
                file, line, cudaGetErrorString(err) );
        exit(-1);
    }
}

static inline void cudaPrintMemory(){
    size_t free, total;
    cudaMemGetInfo(&free,&total);
    std::cout << "Cuda memory: " << float(free)/(1024.0*1024.0) << " / " << float(total)/(1024.0*1024.0) << "MB" << std::endl;
}

#endif /* CONVENIENCE_CUH_ */
