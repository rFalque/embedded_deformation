#include "imageProcessor/common.h"
#include "math/eigen33.h"

namespace device{


    __global__ void CollectValidDataKernel(
        cudaTextureObject_t vertex_map,
        cudaTextureObject_t normal_map,
        const unsigned width,
        const unsigned height,
        PtrSz<Surfel> valid_surfel,
        int* valid_surfel_counter,
        cudaStream_t stream = 0
    )
    {
        const auto x = threadIdx.x + blockIdx.x * blockDim.x;
        const auto y = threadIdx.y + blockIdx.y * blockDim.y;
        if(x >= width || y >= height) return;

        float4 vertex = tex2D<float4>(vertex_map, x, y);
        // 只判断了顶点是否为0 其中许多normal也是0但没有被判断出来
        if (is_zero_vertex(vertex)) return;
        Surfel surfel;
        surfel.pixel_coord.x() = x;
        surfel.pixel_coord.y() = y;
        surfel.vertex = vertex;
        surfel.normal = tex2D<float4>(normal_map, x, y);
        const int index = atomicAdd(valid_surfel_counter, 1);
        valid_surfel[index] = surfel;
    }

    enum {
		window_dim = 7,
		halfsize = 3,
		window_size = window_dim * window_dim
	};

    __global__ void createVertexMapKernel(
        cudaTextureObject_t depth_texture,
        const Intrinsic intrinsic,
        const unsigned width,
        const unsigned height,
        cudaSurfaceObject_t vertex_surface
    )
    {
        const auto x = threadIdx.x + blockIdx.x * blockDim.x;
        const auto y = threadIdx.y + blockIdx.y * blockDim.y;

        if(x >= width || y >= height) return;

        const unsigned short depth =  tex2D<unsigned short>(depth_texture,x,y);
        float4 vertex;

        // scale the depth to [m]
        // bonn 数据集对应的应该是5000
        vertex.z = float(depth) * 0.0010000000474974513;
        vertex.x = vertex.z *(x- intrinsic.principal_x) / intrinsic.focal_x;
        vertex.y = vertex.z *(y- intrinsic.principal_y) / intrinsic.focal_y;
        vertex.w = 0.0f;
        surf2Dwrite(vertex, vertex_surface, x * sizeof(float4), y);
    }

    __global__ void createNormalMapKernel(
        cudaTextureObject_t vertex_map,
        const unsigned width,
        const unsigned height,
        cudaSurfaceObject_t normal_map
    ) {
        
        const auto x = threadIdx.x + blockIdx.x * blockDim.x;
        const auto y = threadIdx.y + blockIdx.y * blockDim.y;
        if(x >= width || y >= height) return;
        
        const float4 vertex_center = tex2D<float4>(vertex_map, x, y);

        float4 normal = make_float4(0, 0, 0, 0);
        if(!is_zero_vertex(vertex_center)) {
            float4 centeroid = make_float4(0, 0, 0, 0);
            int counter = 0;
            for (int cy = y - halfsize; cy <= y + halfsize; cy += 1) {                   // 7*7窗口求和
				for (int cx = x - halfsize; cx <= x + halfsize; cx += 1) {
					const float4 p = tex2D<float4>(vertex_map, cx, cy);
					if (!is_zero_vertex(p)) {
						centeroid.x += p.x;
						centeroid.y += p.y;
						centeroid.z += p.z;
						counter++;
					}
				}
			}

            //At least half of the window is valid
			if(counter > (window_size / 2)) {                                           
				centeroid *= (1.0f / counter);
				float covariance[6] = { 0 };

				//Second window search to compute the normal
				for (int cy = y - halfsize; cy < y + halfsize; cy += 1) {
					for (int cx = x - halfsize; cx < x + halfsize; cx += 1) {
						const float4 p = tex2D<float4>(vertex_map, cx, cy);
						if (!is_zero_vertex(p)) {
							const float4 diff = p - centeroid;
							//Compute the covariance
							covariance[0] += diff.x * diff.x; //(0, 0)
							covariance[1] += diff.x * diff.y; //(0, 1)
							covariance[2] += diff.x * diff.z; //(0, 2)
							covariance[3] += diff.y * diff.y; //(1, 1)
							covariance[4] += diff.y * diff.z; //(1, 2)
							covariance[5] += diff.z * diff.z; //(2, 2)
						}
					}
				}
                eigen33 eigen(covariance);
				float3 normal_value;
				eigen.compute(normal_value);
				if (dotxyz(normal_value, vertex_center) >= 0.0f) normal *= -1;

				//The radius
				const float radius = 0.0;

				//Write to local variable
				normal.x = normal_value.x;
				normal.y = normal_value.y;
				normal.z = normal_value.z;
				normal.w = radius;
            }   
        }
        surf2Dwrite(normal, normal_map, x * sizeof(float4), y);
    }
    
    __global__ void bilateralFilterKernel(
        cudaTextureObject_t raw_depth,
        const unsigned width,
        const unsigned height,
        const unsigned clip_width,
        const unsigned clip_height,
        const float sigma_s_inv_square,
        const float sigma_r_inv_square,
        cudaSurfaceObject_t filter_depth
    )
    {
        const auto x = threadIdx.x + blockIdx.x * blockDim.x;           // 0~640
        const auto y = threadIdx.y + blockIdx.y * blockDim.y;           // 0~480
        if(y >= clip_height || x >= clip_width) return;                 // 0~620 0~460
        
        const auto half_width = 5;
        const auto raw_x = x + boundary_clip;                           // 10~630
        const auto raw_y = y + boundary_clip;                           // 10~470
        const unsigned short center_depth = tex2D<unsigned short> (raw_depth, raw_x, raw_y);

        float sum_all = 0.0f;
        float sum_weight = 0.0f;
        for(auto y_idx = raw_y-half_width; y_idx <= raw_y + half_width; y_idx++)
        {
            for(auto x_idx = raw_x-half_width; x_idx <= raw_x + half_width; x_idx++)
            {
                const unsigned short depth = tex2D<unsigned short> (raw_depth, x_idx, y_idx);
                const float depth_diff2 = (depth -center_depth) * (depth - center_depth);
                const float pixel_diff2 = (x_idx - raw_x) * (x_idx - raw_x) + (y_idx - raw_y) * (y_idx - raw_y);
                const float this_weight = (depth > 0) * expf(-sigma_s_inv_square * pixel_diff2) * expf(-sigma_r_inv_square * depth_diff2);
                sum_weight += this_weight;
                sum_all += this_weight * depth;
            }
        }

        unsigned short filtered_depth_value = __float2uint_rn(sum_all / sum_weight);
        if(filtered_depth_value > 1000) filtered_depth_value = 0;
        surf2Dwrite(filtered_depth_value,filter_depth,x*sizeof(unsigned short),y);
    }
}; // namespace device

void createDefault2DTextureDesc(cudaTextureDesc &desc)
{
    memset(&desc,0,sizeof(desc));
    desc.addressMode[0] = cudaAddressModeBorder;     // 定义了纹理在各个维度（X、Y、Z）上的寻址模式
    desc.addressMode[1] = cudaAddressModeBorder;
    desc.addressMode[2] = cudaAddressModeBorder;
    desc.filterMode = cudaFilterModePoint;
    desc.readMode = cudaReadModeElementType;
    desc.normalizedCoords = 0;
}

// 创建纹理和表面对象
// 三步：1. 创建数组 1.1 创建资源描述符 2. 创建纹理对象 3. 创建表面对象
void createDepthTextureSurface(const unsigned width, const unsigned height, 
                                cudaTextureObject_t& texture, 
                                cudaSurfaceObject_t& surface, 
                                cudaArray_t& array)
{
    // 为array分配内存
    cudaChannelFormatDesc depth_channel_desc = cudaCreateChannelDesc(16,0,0,0,cudaChannelFormatKindUnsigned);
    cudaSafeCall(cudaMallocArray(&array, &depth_channel_desc, width, height));

    // 创建资源描述符
    cudaResourceDesc res_desc;
    memset(&res_desc,0,sizeof(cudaResourceDesc));
    res_desc.resType = cudaResourceTypeArray;         // 指定资源类型为数组
    res_desc.res.array.array = array;                 // 指定数组句柄

    // 创建纹理对象 需要纹理描述符
    cudaTextureDesc depth_texture_desc;
    createDefault2DTextureDesc(depth_texture_desc);
    cudaSafeCall(cudaCreateTextureObject(&texture, &res_desc, &depth_texture_desc, NULL));

    // 创建表面对象
    cudaSafeCall(cudaCreateSurfaceObject(&surface, &res_desc));
}

void createFloat4TextureSurface(const unsigned width, 
                                const unsigned height,
                                cudaTextureObject_t &texture,
                                cudaSurfaceObject_t &surface,
                                cudaArray_t &array)
{
    cudaChannelFormatDesc float4_channel_desc = cudaCreateChannelDesc(32,32,32,32,cudaChannelFormatKindFloat);
    cudaSafeCall(cudaMallocArray(&array, &float4_channel_desc, width, height));

    cudaResourceDesc res_desc;
    memset(&res_desc,0,sizeof(cudaResourceDesc));
    res_desc.resType = cudaResourceTypeArray;
    res_desc.res.array.array = array;  

    cudaTextureDesc float4_texture_desc;
    createDefault2DTextureDesc(float4_texture_desc);
    cudaSafeCall(cudaCreateTextureObject(&texture, &res_desc, &float4_texture_desc, NULL));

    cudaSafeCall(cudaCreateSurfaceObject(&surface, &res_desc));
}

// 释放纹理和表面
void releaseTexture(CudaTextureSurface &texture)
{
    cudaSafeCall(cudaDestroyTextureObject(texture.texture));
    cudaSafeCall(cudaDestroySurfaceObject(texture.surface));
    cudaSafeCall(cudaFreeArray(texture.array));
}

// 深度图双边滤波滤波处理
void BilateralFilterDepth(
    cudaTextureObject_t raw_depth,
    cudaSurfaceObject_t filter_depth,
    const unsigned raw_width,
    const unsigned raw_height,
    const unsigned clip_wdith,
    const unsigned clip_height,
    cudaStream_t stream
)
{
    const float sigma_s = 100.f;
    const float sigma_r = 100.f;
    const float sigma_s_inv_square = 1.0f / (sigma_s * sigma_s);
    const float sigma_r_inv_square = 1.0f / (sigma_r * sigma_r);

    dim3 blk(16,16);
    dim3 grid(divUp(raw_width,blk.x),divUp(raw_height,blk.y));
    device::bilateralFilterKernel<<<grid,blk,0,stream>>>(
        raw_depth,
        raw_width, raw_height,
        clip_wdith, clip_height,
        sigma_s_inv_square, sigma_r_inv_square,
        filter_depth
    );
}

// 计算顶点坐标
void createVertexMap(
    cudaTextureObject_t depth_texture,
    Intrinsic intrinsic,
    const unsigned width,
    const unsigned height,
    cudaSurfaceObject_t vertex_surface,
    cudaStream_t stream
)
{
    dim3 blk(16,16);
    dim3 grid(divUp(width,blk.x),divUp(height,blk.y));
    device::createVertexMapKernel<<<grid,blk,0,stream>>>(
        depth_texture,
        intrinsic,
        width,height,
        vertex_surface
    );
}

void createNormalMap(
    cudaTextureObject_t vertex_map,
    const unsigned width,
    const unsigned height,
    cudaSurfaceObject_t normal_map,
    cudaStream_t stream
) {
    dim3 blk(16,16);
    dim3 grid(divUp(width,blk.x),divUp(height,blk.y));
    device::createNormalMapKernel<<<grid, blk, 0, stream>>>(
        vertex_map, width, height, normal_map
    );
}

void CollectValidData(
    cudaTextureObject_t vertex_map,
    cudaTextureObject_t normal_map,
    const unsigned width,
    const unsigned height,
    DeviceArray<Surfel>& valid_surfel,
    int& surfel_counter,
    cudaStream_t stream
)
{
    int* valid_surfel_counter;
    cudaMalloc((void**)& valid_surfel_counter, sizeof(int));
    cudaMemset(valid_surfel_counter, 0, sizeof(int));
    dim3 blk(16,16);
    dim3 grid(divUp(width,blk.x),divUp(height,blk.y));
    device::CollectValidDataKernel<<<grid,blk,0,stream>>>(
        vertex_map, 
        normal_map, 
        width, height, 
        valid_surfel, 
        valid_surfel_counter
    );
    int h_surfel_counter;
    cudaSafeCall(cudaMemcpy(&h_surfel_counter, valid_surfel_counter, sizeof(int), cudaMemcpyDeviceToHost));
    surfel_counter = h_surfel_counter;
    cudaFree(valid_surfel_counter);
}
