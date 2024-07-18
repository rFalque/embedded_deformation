#ifndef __SAFE_CALL_
#define __SAFE_CALL_

#include "cuda_runtime.h"
#include <iostream>
#include <stdlib.h>
#include <cublas_v2.h>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <cuda.h>

#if defined(__GNUC__)
    #define cudaSafeCall(expr)  ___cudaSafeCall(expr, __FILE__, __LINE__, __func__)
#else /* defined(__CUDACC__) || defined(__MSVC__) */
    #define cudaSafeCall(expr)  ___cudaSafeCall(expr, __FILE__, __LINE__)    
#endif


static inline void ___cudaSafeCall(cudaError_t err, const char *file, const int line, const char *func = "")
{
    if (cudaSuccess != err){
        std::cout << "Error: " << err <<"\t" << file << ":" << line << std::endl;
        printf("%s",cudaGetErrorString(err));
        std::cout<<std::endl;
        exit(0);
    }
}

static inline int divUp(int total, int grain) 
{ 
    return (total + grain - 1) / grain; 
}

#endif 

#ifndef cublasSafeCall
#define cublasSafeCall(err)  __cublasSafeCall(err, __FILE__, __LINE__)
static inline void __cublasSafeCall(cublasStatus_t err, const char *file, const int line)
{
    if(CUBLAS_STATUS_SUCCESS != err) {
        fprintf(stderr, "CUBLAS error in file '%s', line %d\n \nerror %d \nterminating!\n",__FILE__, __LINE__, err);
        cudaDeviceReset();
        std::exit(1);
    }
}
#endif

#ifndef cuSafeCall
#define cuSafeCall(err) __cuSafeCall(err, __FILE__, __LINE__)
static inline void __cuSafeCall(CUresult err, const char* file, const int line) {
    if(err != CUDA_SUCCESS) {
        //Query the name and string of the error
        const char* error_name;
        cuGetErrorName(err, &error_name);
        const char* error_string;
        cuGetErrorString(err, &error_string);
        fprintf(stderr, "CUDA driver error %s: %s in the line %d of file %s \n", error_name, error_string, line, file);
        cudaDeviceReset();
        std::exit(1);
    }
}
	
#endif





