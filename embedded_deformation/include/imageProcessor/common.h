#pragma once
#include "common_types.hpp"
#include "cuda_runtime_api.h"
#include "safe_call.hpp"
#include "global_configs.h"
#include <vector>
#include <glm/vec3.hpp>
#include <math/vector_ops.hpp>
#include "containers/device_array.hpp"

using namespace CUDA;

// 创建深度图纹理表面
void createDepthTextureSurface(const unsigned width, 
                            const unsigned height, 
                            cudaTextureObject_t& texture, 
                            cudaSurfaceObject_t& surface,
                            cudaArray_t& array);
// 创建纹理描述资源
void createDefault2DTextureDesc(cudaTextureDesc &desc);

// 释放纹理内存
void releaseTexture(CudaTextureSurface &texture);


// 深度图处理 滤波
void BilateralFilterDepth(
    cudaTextureObject_t raw_depth,
    cudaSurfaceObject_t filter_depth,
    const unsigned raw_width,
    const unsigned raw_height,
    const unsigned clip_width,
    const unsigned clip_height,
    cudaStream_t stream = 0
);

void createFloat4TextureSurface(const unsigned width, 
                                const unsigned height,
                                cudaTextureObject_t &texture,
                                cudaSurfaceObject_t &surface,
                                cudaArray_t &array
                                );

// 计算顶点坐标
void createVertexMap(
    cudaTextureObject_t depth_texture,
    Intrinsic intrinsic,
    const unsigned width,
    const unsigned height,
    cudaSurfaceObject_t vertex_surface,
    cudaStream_t stream = 0
);

// 计算法向量
void createNormalMap(
    cudaTextureObject_t vertex_map,
    const unsigned width,
    const unsigned height,
    cudaSurfaceObject_t normal_map,
    cudaStream_t stream = 0
);

// 收集有效数据
void CollectValidData(
    cudaTextureObject_t vertex_map,
    cudaTextureObject_t normal_map,
    const unsigned width,
    const unsigned height,
    DeviceArray<Surfel>& valid_surfel,
    int& surfel_count,
    cudaStream_t stream = 0
);
