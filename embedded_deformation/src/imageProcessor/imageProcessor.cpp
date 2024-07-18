#include "imageProcessor/imageProcessor.h"

imageProcessor::imageProcessor( Intrinsic raw_Intrinsic,
                                const unsigned width, 
                                const unsigned height, 
                                int start_frame,
                                std::shared_ptr<FetchInterface> image_fetcher
                                ):
m_raw_intrinsic(raw_Intrinsic),
m_raw_width(width),
m_raw_height(height),
m_clip_width(width-2*boundary_clip),
m_clip_height(height-2*boundary_clip),
m_frame_idx(start_frame),
m_image_fetcher(image_fetcher)
{   
    m_clip_intrinsic = Intrinsic(m_raw_intrinsic.focal_x, 
                                m_raw_intrinsic.focal_y,
                                m_raw_intrinsic.principal_x-boundary_clip, 
                                m_raw_intrinsic.principal_y-boundary_clip);

    // m_depth_image = cv::imread(m_depth_path, cv::IMREAD_UNCHANGED);
    // m_mask_image = cv::imread(m_mask_path, cv::IMREAD_GRAYSCALE);
    // if(m_depth_image.empty()) {
    //     std::cerr << "Error: Unable to load depth image from " << m_depth_path << std::endl;
    // }
    // if (m_mask_image.empty()) {
    //     std::cerr << "Error: Unable to load mask image from " << m_mask_path<< std::endl;
    // }
}


imageProcessor::~imageProcessor()
{
    releaseDepthTexture();
    releaseFetchBuffer();
    releaseVertexTexture();
    releaseNormalTexture();
}

// 分配与释放内存
// 分配用于主机与设备之间数据传输的图像和锁页
void imageProcessor::allocateFetchBuffer()
{
    const auto raw_img_size = m_raw_width * m_raw_height;
    cudaSafeCall(cudaMallocHost(&m_depth_buffer_pagelock, sizeof(unsigned short)*raw_img_size));
    cudaSafeCall(cudaMallocHost(&m_depth_prev_buffer_pagelock,sizeof(unsigned short)*raw_img_size));
    cudaSafeCall(cudaDeviceSynchronize());
    cudaSafeCall(cudaGetLastError());
    m_depth_image = cv::Mat(cv::Size(m_raw_width, m_raw_height), CV_16UC1);
    m_depth_image_prev = cv::Mat(cv::Size(m_raw_width, m_raw_height), CV_16UC1);
}

void imageProcessor::releaseFetchBuffer()
{
    cudaSafeCall(cudaFreeHost(m_depth_buffer_pagelock));
    // cudaSafeCall(cudaFreeHost(m_rgb_buffer_packlock));
    cudaSafeCall(cudaFreeHost(m_depth_prev_buffer_pagelock));
    cudaSafeCall(cudaDeviceSynchronize());
    cudaSafeCall(cudaGetLastError());
}

void imageProcessor::allocateDepthTexture()
{
    createDepthTextureSurface(m_raw_width, m_raw_height,
                                depth_texture_surface_raw.texture,
                                depth_texture_surface_raw.surface,
                                depth_texture_surface_raw.array);
    
    createDepthTextureSurface(m_raw_width, m_raw_height,
                                depth_texture_surface_raw_prev.texture,
                                depth_texture_surface_raw_prev.surface,
                                depth_texture_surface_raw_prev.array);

    createDepthTextureSurface(m_clip_width, m_clip_height,
                                depth_texture_surface_filter.texture,
                                depth_texture_surface_filter.surface,
                                depth_texture_surface_filter.array);
    
    createDepthTextureSurface(m_clip_width, m_clip_height,
                                depth_texture_surface_filter_prev.texture,
                                depth_texture_surface_filter_prev.surface,
                                depth_texture_surface_filter_prev.array);
}

void imageProcessor::releaseDepthTexture()
{
    releaseTexture(depth_texture_surface_raw);
    releaseTexture(depth_texture_surface_filter);
    releaseTexture(depth_texture_surface_raw_prev);
    releaseTexture(depth_texture_surface_filter_prev);

    if(h_depthData != nullptr)
    {
        delete[] h_depthData;
        h_depthData = nullptr;
    }

}

void imageProcessor::allocateVertexTexture()
{
    createFloat4TextureSurface(m_clip_width, m_clip_height,
                                vertex_texture_surface.texture,
                                vertex_texture_surface.surface,
                                vertex_texture_surface.array);
    
    createFloat4TextureSurface(m_clip_width, m_clip_height,
                                vertex_texture_surface_prev.texture,
                                vertex_texture_surface_prev.surface,
                                vertex_texture_surface_prev.array);
}

void imageProcessor::releaseVertexTexture()
{
    releaseTexture(vertex_texture_surface);
    releaseTexture(vertex_texture_surface_prev);

    if(h_vertexData != nullptr){
        delete[] h_vertexData;
        h_vertexData = nullptr;
    }

    if(h_vertexData_prev != nullptr){
        delete[] h_vertexData_prev;
        h_vertexData_prev = nullptr;
    }
}

void imageProcessor::allocateNormalTexture()
{
    createFloat4TextureSurface(m_clip_width, m_clip_height,
                                normal_texture_surface.texture,
                                normal_texture_surface.surface,
                                normal_texture_surface.array);
    
    createFloat4TextureSurface(m_clip_width, m_clip_height,
                                normal_texture_surface_prev.texture,
                                normal_texture_surface_prev.surface,
                                normal_texture_surface_prev.array);
}

void imageProcessor::releaseNormalTexture()
{
    releaseTexture(normal_texture_surface);
    releaseTexture(normal_texture_surface_prev);

    if(h_normalData != nullptr){
        delete[] h_normalData;
        h_normalData = nullptr;
    }

    if(h_normalData_prev != nullptr){
        delete[] h_normalData_prev;
        h_normalData_prev = nullptr;
    }
}

void imageProcessor::allocateSurfelBuffer()
{
    const auto num_pixels = m_clip_height*m_clip_width;
    m_surfels.AllocateBuffer(num_pixels);
    m_surfels_prev.AllocateBuffer(num_pixels);
}

// 图像抓取
void imageProcessor::FetchDepthImage(size_t frame_idx)
{
    m_image_fetcher->FetchDepthImage(frame_idx, m_depth_image);
    memcpy(m_depth_buffer_pagelock, m_depth_image.data, 
            sizeof(unsigned short)*m_raw_width*m_raw_height);
    cudaSafeCall(cudaGetLastError());
}

void imageProcessor::FetchDepthPrevImage(size_t frame_idx)
{
    size_t prev_frame_idx = frame_idx - 10;
    m_image_fetcher->FetchDepthImage(prev_frame_idx, m_depth_image_prev);
    memcpy(m_depth_prev_buffer_pagelock, m_depth_image_prev.data,
            sizeof(unsigned short)*m_raw_width*m_raw_height);
}

void imageProcessor::FetchRGBImage(size_t frame_idx)
{
    m_image_fetcher->FetchRGBImage(frame_idx, m_rgb_image);
    memcpy(m_rgb_buffer_pagelock, m_rgb_image.data,
                sizeof(uchar3) * m_raw_width * m_raw_height);
}

void imageProcessor::FetchRGBPrevImage(size_t frame_idx)
{
    size_t prev_frame_idx = frame_idx - 3;
    m_image_fetcher->FetchRGBImage(prev_frame_idx, m_rgb_image_prev);
    memcpy(m_rgb_prev_buffer_pagelock, m_rgb_image_prev.data,
            sizeof(uchar3)*m_raw_width*m_raw_height);
}

void imageProcessor::FetchFrame(size_t frame_idx)
{   
    if(frame_idx == 0){
        no_prev_frame = true;
        return;
    } 

    FetchDepthImage(frame_idx);

    if(!no_prev_frame){
        FetchDepthPrevImage(frame_idx);
    }
    
}

void imageProcessor::UploadDepthImage(cudaStream_t stream)
{
    cudaSafeCall(cudaMemcpyToArrayAsync(
        depth_texture_surface_raw.array,
        0,0,
        m_depth_buffer_pagelock,
        sizeof(unsigned short)*m_raw_width*m_raw_height,
        cudaMemcpyHostToDevice,
        stream
    ));
    
    if(!no_prev_frame){
        cudaSafeCall(cudaMemcpyToArrayAsync(
            depth_texture_surface_raw_prev.array,
            0,0,
            m_depth_prev_buffer_pagelock,
            sizeof(unsigned short)*m_raw_width*m_raw_height,
            cudaMemcpyHostToDevice,
            stream
        )); 
    }
    
}

void imageProcessor::FilterDepthImage(cudaStream_t stream)
{
    BilateralFilterDepth(
        depth_texture_surface_raw.texture,
        depth_texture_surface_filter.surface,
        m_raw_width, m_raw_height,
        m_clip_width, m_clip_height,
        stream
    );

    if(!no_prev_frame){
        BilateralFilterDepth(
            depth_texture_surface_raw_prev.texture,
            depth_texture_surface_filter_prev.surface,
            m_raw_width, m_raw_height,
            m_clip_width, m_clip_height,
            stream
        );
    }
}

void imageProcessor::BuildVertexMap(cudaStream_t stream)
{
    createVertexMap(
        depth_texture_surface_filter.texture,
        m_clip_intrinsic,
        m_clip_width, m_clip_height,
        vertex_texture_surface.surface,
        stream
    );

    if(!no_prev_frame){
        createVertexMap(
            depth_texture_surface_filter_prev.texture,
            m_clip_intrinsic,
            m_clip_width, m_clip_height,
            vertex_texture_surface_prev.surface,
            stream
        );
    }
}

void imageProcessor::BuildNormalMap(cudaStream_t stream)
{
    createNormalMap(
        vertex_texture_surface.texture,
        m_clip_width, m_clip_height,
        normal_texture_surface.surface,
        stream
    ); 

    if(!no_prev_frame){
        createNormalMap(
            vertex_texture_surface_prev.texture,
            m_clip_width, m_clip_height,
            normal_texture_surface_prev.surface,
            stream
        );
    }
}

void imageProcessor::DownloadDepthData(cudaStream_t stream)
{
    h_depthData = (unsigned short*)malloc(sizeof(unsigned short)*m_clip_width*m_clip_height);
    cudaSafeCall(cudaMemcpyFromArray(h_depthData, depth_texture_surface_filter_prev.array, 
                        0, 0, sizeof(unsigned short)*m_clip_width*m_clip_height, 
                        cudaMemcpyDeviceToHost));
    
    cv::Mat filtered_depth(m_clip_height, m_clip_width, CV_16UC1, h_depthData);
    cv::Mat normalizedFilteredDepth;
    cv::normalize(filtered_depth, normalizedFilteredDepth, 0, 255, cv::NORM_MINMAX, CV_8UC1);
    cv::namedWindow("Normalized Filtered Depth Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Normalized Filtered Depth Image", normalizedFilteredDepth);
    
    cv::Mat normalizedDepth;
    cv::normalize(m_depth_image_prev, normalizedDepth, 0, 255, cv::NORM_MINMAX, CV_8UC1);
    cv::namedWindow("Normalized Depth Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Normalized Depth Image", normalizedDepth);
    cv::waitKey(0); // 等待按键继续
    cv::destroyAllWindows(); // 关闭所有窗口
}

void imageProcessor::DownloadVertexNormalData(cudaStream_t stream)
{
    
    // h_vertexData = (float4*)malloc(sizeof(float4)*m_clip_width*m_clip_height);
    // cudaSafeCall(cudaMemcpyFromArray(h_vertexData, vertex_texture_surface.array, 
    //                     0, 0, sizeof(float4)*m_clip_width*m_clip_height, 
    //                     cudaMemcpyDeviceToHost));
    
    // h_vertexData_prev = (float4*)malloc(sizeof(float4)*m_clip_width*m_clip_height);
    // cudaSafeCall(cudaMemcpyFromArray(h_vertexData_prev, vertex_texture_surface_prev.array, 
    //                     0, 0, sizeof(float4)*m_clip_width*m_clip_height, 
    //                     cudaMemcpyDeviceToHost));
    
    // h_normalData = (float4*) malloc(sizeof(float4)* m_clip_width* m_clip_height);
    // cudaSafeCall(cudaMemcpyFromArray(h_normalData, normal_texture_surface.array,
    //                     0, 0, sizeof(float4)* m_clip_width* m_clip_height,
    //                     cudaMemcpyDeviceToHost));
    // h_normalData_prev = (float4*) malloc(sizeof(float4)* m_clip_width* m_clip_height);
    // cudaSafeCall(cudaMemcpyFromArray(h_normalData_prev, normal_texture_surface_prev.array,
    //                     0, 0, sizeof(float4)* m_clip_width* m_clip_height,
    //                     cudaMemcpyDeviceToHost));

    
    // std::vector<glm::vec3> points;
    // std::vector<glm::vec3> points_prev;
    // std::vector<glm::vec3> normals;
    // std::vector<glm::vec3> normals_prev;

    // for(int i=0; i<m_clip_height; i++)
    // {
    //     for(int j=0; j<m_clip_width; j++)
    //     {   
    //         float x = h_vertexData[i*m_clip_width+j].x;
    //         float y = h_vertexData[i*m_clip_width+j].y;
    //         float z = h_vertexData[i*m_clip_width+j].z;
            // float x_prev = h_vertexData_prev[i*m_clip_width+j].x;
            // float y_prev = h_vertexData_prev[i*m_clip_width+j].y;
            // float z_prev = h_vertexData_prev[i*m_clip_width+j].z;
            // float nx = h_normalData[i*m_clip_width+j].x;
            // float ny = h_normalData[i*m_clip_width+j].y;
            // float nz = h_normalData[i*m_clip_width+j].z;
            // float nx_prev = h_normalData_prev[i*m_clip_width+j].x;
            // float ny_prev = h_normalData_prev[i*m_clip_width+j].y;
            // float nz_prev = h_normalData_prev[i*m_clip_width+j].z;

            // points.push_back(glm::vec3(x, y, z));
            // normals.push_back(glm::vec3(nx, ny, nz));
            // points_prev.push_back(glm::vec3(x_prev, y_prev, z_prev));
            // normals_prev.push_back(glm::vec3(nx_prev, ny_prev, nz_prev));
    //     }
    // }
    // polyscope::init();
    // 注册第一组点云并设置为红色
    // auto* psCloud1 = polyscope::registerPointCloud("Point Cloud 1", points);
    // polyscope::getPointCloud("Point Cloud 1")->addVectorQuantity("normals 1", normals);
    // psCloud1->setPointColor(glm::vec3(1.0, 0.0, 0.0)); // RGB颜色，这里是红色

    // // 注册第二组点云并设置为蓝色
    // auto* psCloud2 = polyscope::registerPointCloud("Point Cloud 2", points_prev);
    // psCloud2->setPointColor(glm::vec3(0.0, 0.0, 1.0)); // RGB颜色，这里是蓝色
    // polyscope::getPointCloud("Point Cloud 2")->addVectorQuantity("normals 2", normals_prev);
    // polyscope::show();

    
}

void imageProcessor::CollectValidSurfelData(cudaStream_t stream)
{
    int valid_surfel_count = 0;
    DeviceArray<Surfel> valid_surfels = m_surfels.Array();
    CollectValidData(
        vertex_texture_surface.texture,
        normal_texture_surface.texture,
        m_clip_width, m_clip_height,
        valid_surfels,
        valid_surfel_count,
        stream
    );
    m_surfels.ResizeArrayOrException(valid_surfel_count);

    int valid_surfel_count_prev = 0;
    DeviceArray<Surfel> valid_surfels_prev = m_surfels_prev.Array();
    if(!no_prev_frame) {
    CollectValidData(
        vertex_texture_surface_prev.texture,
        normal_texture_surface_prev.texture,
        m_clip_width, m_clip_height,
        valid_surfels_prev,
        valid_surfel_count_prev,
        stream
    );
    m_surfels_prev.ResizeArrayOrException(valid_surfel_count_prev);
    }
}

void imageProcessor::test()
{
    allocateFetchBuffer();
    allocateDepthTexture();
    allocateVertexTexture();
    allocateNormalTexture();
    allocateSurfelBuffer();
    FetchFrame(m_frame_idx);
    UploadDepthImage();
    FilterDepthImage();
    // DownloadDepthData()
    BuildVertexMap();
    BuildNormalMap();
    CollectValidSurfelData();
    // DownloadVertexNormalData();
    
    
}
