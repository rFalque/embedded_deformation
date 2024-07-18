#pragma once
#include <opencv2/opencv.hpp>
#include "common_types.hpp"
#include "common.h"
#include "fetchInterface.h"
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "containers/DeviceBufferArray.h"

using namespace cv;
using namespace CUDA;

class imageProcessor
{
public:
    imageProcessor( Intrinsic intrinsic, 
                    const unsigned width, 
                    const unsigned height, 
                    int start_frame,
                    std::shared_ptr<FetchInterface> image_fetcher
                    );
    ~imageProcessor();

private:
    // 从主机读取的图像
    cv::Mat m_depth_image, m_depth_image_prev;
    cv::Mat m_rgb_image, m_rgb_image_prev;
    
    // 相机内参 
    Intrinsic m_raw_intrinsic;
    // 修剪后的相机内参
    Intrinsic m_clip_intrinsic;

    // 批量图像读取对象
    std::shared_ptr<FetchInterface> m_image_fetcher;

    // 原始图像的宽和高
    const unsigned m_raw_width;
    const unsigned m_raw_height;
    // 修剪后图像的宽和高
    const unsigned m_clip_width;
    const unsigned m_clip_height;

    // 锁页内存 从主机到设备数据传输的中转站 
    // 它提供了更快的主机到设备的数据传输
    void* m_depth_buffer_pagelock;
    void* m_depth_prev_buffer_pagelock;
    void* m_rgb_buffer_pagelock;
    void* m_rgb_prev_buffer_pagelock;

    // 纹理和表面
    // 用于存储原始深度图 
    CudaTextureSurface depth_texture_surface_raw;
    CudaTextureSurface depth_texture_surface_raw_prev;
    // 用于存储滤波后的深度图
    CudaTextureSurface depth_texture_surface_filter;
    CudaTextureSurface depth_texture_surface_filter_prev;
    // 用于存储顶点图
    CudaTextureSurface vertex_texture_surface;
    CudaTextureSurface vertex_texture_surface_prev;
    // 用于存储法向量图
    CudaTextureSurface normal_texture_surface;
    CudaTextureSurface normal_texture_surface_prev;
    // surfel
    DeviceBufferArray<Surfel> m_surfels;
    DeviceBufferArray<Surfel> m_surfels_prev;

    // 用于下载设备端的数据
    unsigned short* h_depthData;
    float4* h_vertexData;
    float4* h_vertexData_prev;
    float4* h_normalData;
    float4* h_normalData_prev;

    // 当前图像的索引
    int m_frame_idx;

    bool no_prev_frame = false;

public:

    // 分配用于读取图像的内存
    void allocateFetchBuffer();
    // 释放用于读取图像的内存
    void releaseFetchBuffer();
    // 分配存储深度图的纹理和表面的内存
    void allocateDepthTexture();
    // 释放存储深度图的纹理和表面的内存
    void releaseDepthTexture();   
    // 分配存储顶点图的纹理和表面的内存
    void allocateVertexTexture();
    // 释放存储顶点图的纹理和表面的内存
    void releaseVertexTexture();
    // 分配存储法向量图的纹理和表面的内存
    void allocateNormalTexture();
    // 释放存储法向量图的纹理和表面的内存
    void releaseNormalTexture();
    // 分配存储Surfel的内存
    void allocateSurfelBuffer();

    // 读取图像
    void FetchFrame(size_t frame_idx);
    // 读取深度图
    void FetchDepthImage(size_t frame_idx);
    void FetchDepthPrevImage(size_t frame_idx);
    // 读取彩色图
    void FetchRGBImage(size_t frame_idx);
    void FetchRGBPrevImage(size_t frame_idx);

    // 将深度图上传到设备
    void UploadDepthImage(cudaStream_t stream = 0);
    // 深度图滤波 这里滤掉了>1000的深度值
    void FilterDepthImage(cudaStream_t stream = 0);
    // 构建顶点坐标，需要根据数据集设置scale_factor
    void BuildVertexMap(cudaStream_t stream = 0);
    // 根据顶点计算法向量
    void BuildNormalMap(cudaStream_t stream = 0);

    // 将设备上的数据下载到主机
    void DownloadDepthData(cudaStream_t stream = 0);
    void DownloadVertexNormalData(cudaStream_t stream = 0);

    // 收集有效的surfel数据 包括Vertx和Normal
    void CollectValidSurfelData(cudaStream_t stream = 0);

    // 接口函数
    vertex_normal_maps getVertexNormalMaps() {return {vertex_texture_surface.texture, normal_texture_surface.texture};}
    vertex_normal_maps getVertexNormalMapsPrev() {return {vertex_texture_surface_prev.texture, normal_texture_surface_prev.texture};}
    unsigned clip_width() const { return m_clip_width; }
    unsigned clip_height() const { return m_clip_height; }
    Intrinsic clip_intrinsic() const { return m_clip_intrinsic; }
    float4* getVertexData() const { return h_vertexData; }
    std::vector<Surfel> getSurfelData() const { DeviceArrayView<Surfel> surfel_view = m_surfels.ArrayView();
                                                std::vector<Surfel> h_surfels;
                                                surfel_view.Download(h_surfels);
                                                return h_surfels;
                                                }
    std::vector<Surfel> getSurfelDataPrev() const { DeviceArrayView<Surfel> surfel_view = m_surfels_prev.ArrayView();
                                                std::vector<Surfel> h_surfels;
                                                surfel_view.Download(h_surfels);
                                                return h_surfels;
                                                }

    // 测试函数
    void test();


};