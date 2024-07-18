#ifndef KERNEL_CONTAINERS_HPP_
#define KERNEL_CONTAINERS_HPP_

#include <cstddef>
#include "imageProcessor/safe_call.hpp"

#if defined(__CUDACC__)
#define GPU_HOST_DEVICE__ __host__ __device__ __forceinline__
#else
#define GPU_HOST_DEVICE__
#endif

namespace CUDA{
  
template <typename T>
struct DevPtr {
  typedef T elem_type;
  const static size_t elem_size = sizeof(elem_type);

  T* data;

  GPU_HOST_DEVICE__ DevPtr() : data(0) {}                                        
  GPU_HOST_DEVICE__ DevPtr(T* data_arg) : data(data_arg) {}                       

  GPU_HOST_DEVICE__ size_t elemSize() const { return elem_size; }                 
  GPU_HOST_DEVICE__ operator T*() { return data; }
  GPU_HOST_DEVICE__ operator const T*() const { return data; }
};

// 继承父类 加入元素数量
template <typename T> 
struct PtrSz : public DevPtr<T> {
  GPU_HOST_DEVICE__ PtrSz() : size(0) {}
  GPU_HOST_DEVICE__ PtrSz(T* data_arg, size_t size_arg) : DevPtr<T>(data_arg), size(size_arg) {}

  size_t size;
};

// 继承父类 加入表示两个连续行之间的步长 适合用于表示二维数组
template <typename T>
struct PtrStep : public DevPtr<T> {
  GPU_HOST_DEVICE__ PtrStep() : step(0) {}
  GPU_HOST_DEVICE__ PtrStep(T* data_arg, size_t step_arg) : DevPtr<T>(data_arg), step(step_arg) {}              // 调用基类的构造函数进行初始化

  /** \brief stride between two consecutive rows in bytes. Step is stored always and everywhere in bytes!!! */
  size_t step;

  GPU_HOST_DEVICE__ T* ptr(int y = 0) { return (T*)((char*)DevPtr<T>::data + y * step); }                          // 用于获取指向指定行的指针
  GPU_HOST_DEVICE__ const T* ptr(int y = 0) const { return (const T*)((const char*)DevPtr<T>::data + y * step); }
};

// 用于表示和操作大小和布局都已知的二维数组
template <typename T>
struct PtrStepSz : public PtrStep<T> {
  GPU_HOST_DEVICE__ PtrStepSz() : cols(0), rows(0) {}
  GPU_HOST_DEVICE__ PtrStepSz(int rows_arg, int cols_arg, T* data_arg, size_t step_arg)                           // 指定行，列，数据类型，
      : PtrStep<T>(data_arg, step_arg), cols(cols_arg), rows(rows_arg) {}

  int cols;
  int rows;
};
}




#endif /* KERNEL_CONTAINERS_HPP_ */
