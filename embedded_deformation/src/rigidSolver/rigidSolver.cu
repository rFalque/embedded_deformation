#include "imageProcessor/common_types.hpp"
#include "math/device_mat.h"
#include "rigidSolver/device_intrinsics.h"
#include "rigidSolver/rigidSolver.h"

namespace device {
    struct RigidSolverDevice {

		//The constants for matrix size and blk size
		enum
		{
			//The stoge layout
            // 27个数据项
			lhs_matrix_size = 21,
			rhs_vector_size = 6,
			total_shared_size = lhs_matrix_size + rhs_vector_size,

			//The block size
			block_size = 256,
			num_warps = block_size / 32,
		};
		
		//The map from the renderer
		struct {
			cudaTextureObject_t vertex_map;
			cudaTextureObject_t normal_map;
		} reference_maps;

		//The map from the depth image
		struct {
			cudaTextureObject_t vertex_map;
			cudaTextureObject_t normal_map;
		} live_maps;

		//The camera information
		mat34 init_refer2live;
		Intrinsic intrinsic;

		//The image information
		unsigned image_rows;
		unsigned image_cols;

		//The processing interface
		__device__ __forceinline__ void solverIteration(
			PtrStep<float> reduce_buffer
		) const {
			const auto flatten_pixel_idx = threadIdx.x + blockDim.x * blockIdx.x;
			const auto x = flatten_pixel_idx % image_cols;
			const auto y = flatten_pixel_idx / image_cols;

			//Prepare the jacobian and err
			float jacobian[6] = {0};
			float err = 0.0f;

			// 范围检查
			if(x < image_cols && y < image_rows)
			{
				// 参考帧 应该是从全局模型渲染得到
				const float4 reference_v4 = tex2D<float4>(reference_maps.vertex_map, x, y);                           // 参考帧顶点的索引
				const float4 reference_n4 = tex2D<float4>(reference_maps.normal_map, x, y);

				// 转换到相机坐标系下
				const auto reference_v = init_refer2live.rot * reference_v4 + init_refer2live.trans;
				const auto reference_n = init_refer2live.rot * reference_n4;

				// 反投影 当前帧顶点的索引
				const ushort2 img_coord = {
					__float2uint_rn(((reference_v.x / (reference_v.z + 1e-10)) * intrinsic.focal_x) + intrinsic.principal_x),
					__float2uint_rn(((reference_v.y / (reference_v.z + 1e-10)) * intrinsic.focal_y) + intrinsic.principal_y)
				};

				// 范围检查
				if(img_coord.x < image_cols && img_coord.y < image_rows)
				{
					// 当前帧定点图
					const float4 live_v4 = tex2D<float4>(live_maps.vertex_map, img_coord.x, img_coord.y);
					const float4 live_n4 = tex2D<float4>(live_maps.normal_map, img_coord.x, img_coord.y);

					//Check correspondence
					if(dotxyz(reference_n, live_n4) < 0.8f || squared_distance(reference_v, live_v4) > (0.01f * 0.01f) || is_zero_vertex(live_v4)) {
						//Pass
					}
					else {
		
						err = dotxyz(live_n4, make_float4(live_v4.x - reference_v.x, live_v4.y - reference_v.y,  live_v4.z - reference_v.z, 0.0f));
						*(float3*)jacobian = cross_xyz(reference_v, live_n4);
						*(float3*)(jacobian + 3) = make_float3(live_n4.x, live_n4.y, live_n4.z);
		
					}
				}
			}

			//Time to do reduction
			__shared__ float reduce_mem[total_shared_size][num_warps];
			unsigned shift = 0;
			const auto warp_id = threadIdx.x >> 5;                              // 当前线程所在的warp ID
			const auto lane_id = threadIdx.x & 31;                              // 当前线程在其所在warp中的位置

			//Reduce on matrix
			for (int i = 0; i < 6; i++) { //Row index
				for (int j = i; j < 6; j++) { //Column index, the matrix is symmetry
					float data = (jacobian[i] * jacobian[j]);
					data = warp_scan(data);
					// 检查当前线程是否是warp中的最后一个线程
					if (lane_id == 31) {
						reduce_mem[shift++][warp_id] = data;
					}
					//Another sync here for reduced mem
					__syncthreads();
				}
			}

			//Reduce on vector
			for (int i = 0; i < 6; i++) {
				float data = (err * jacobian[i]);
				data = warp_scan(data);
				if (lane_id == 31) {
					reduce_mem[shift++][warp_id] = data;
				}
				//Another sync here for reduced mem
				__syncthreads();
			}

			//Store the result to global memory
			const auto flatten_blk = blockIdx.x;
			for (int i = threadIdx.x; i < total_shared_size; i += 32) {
				if (warp_id == 0) {
					const auto warp_sum = reduce_mem[i][0] + reduce_mem[i][1] + reduce_mem[i][2] + reduce_mem[i][3] 
							+ reduce_mem[i][4] + reduce_mem[i][5] + reduce_mem[i][6] + reduce_mem[i][7];
					reduce_buffer.ptr(i)[flatten_blk] = warp_sum;
				}
			}
		}
	};

    __global__ void rigidSolveIterationKernel(
		const RigidSolverDevice solver,
		PtrStep<float> reduce_buffer
	) {
		solver.solverIteration(reduce_buffer);
	}

    __global__ void columnReduceKernel(
        const PtrStepSz<const float> global_buffer,
        float* target
    ) {
        const auto idx = threadIdx.x;
        const auto y = threadIdx.y + blockIdx.y * blockDim.y;
        float sum = 0.0f;
        for (auto i = threadIdx.x; i < global_buffer.cols; i += 32) {
            sum += global_buffer.ptr(y)[i];
        }

        sum = warp_scan(sum);
        if(idx==31){
            target[y] = sum;
        }
    }

	__global__ void FindCorresponceKernel(
		cudaTextureObject_t reference_vertex_map,
		cudaTextureObject_t reference_normal_map,
		cudaTextureObject_t live_vertex_map,
		cudaTextureObject_t live_normal_map,
		PtrSz<correspondence> correspondence,
		const unsigned width,
		const unsigned height,
		const Intrinsic intrinsic,
		int* correspondence_counter,
		const mat34 refer2live
	) {
		const auto x = threadIdx.x + blockIdx.x * blockDim.x;
		const auto y = threadIdx.y + blockIdx.y * blockDim.y;
		if( x >= width || y >= height) return;
		const float4 reference_v4 = tex2D<float4>(reference_vertex_map, x, y);
		const float4 reference_n4 = tex2D<float4>(reference_normal_map, x, y);
		const float3 reference_v = refer2live.rot * reference_v4 + refer2live.trans;
		const float3 reference_n = refer2live.rot * reference_n4;

		const ushort2 img_coord = {
					__float2uint_rn(((reference_v.x / (reference_v.z + 1e-10)) * intrinsic.focal_x) + intrinsic.principal_x),
					__float2uint_rn(((reference_v.y / (reference_v.z + 1e-10)) * intrinsic.focal_y) + intrinsic.principal_y)
				};
		
		
		if(img_coord.x < width && img_coord.y < height)
		{
			// 当前帧定点图
			const float4 live_v4 = tex2D<float4>(live_vertex_map, img_coord.x, img_coord.y);
			const float4 live_n4 = tex2D<float4>(live_normal_map, img_coord.x, img_coord.y);
			//Check correspondence
			if(dotxyz(reference_n, live_n4) < 0.8f || squared_distance(reference_v, live_v4) > (0.01f * 0.01f) || is_zero_vertex(live_v4)) {
				//Pass
			}
			else {
				const int index = atomicAdd(correspondence_counter, 1);
				correspondence[index].src_vertex = make_float4(reference_v.x, reference_v.y, reference_v.z, 0.0f);
				correspondence[index].src_normal = make_float4(reference_n.x, reference_n.y, reference_n.z, 0.0f);
				correspondence[index].tag_vertex = live_v4;
				correspondence[index].tag_normal = live_n4;
			}
		}
	}

	__global__ void FindCorresponceKernel(
		double* src_vertex,
		double* src_normal,
		cudaTextureObject_t live_vertex_map,
		cudaTextureObject_t live_normal_map,
		const unsigned width,
		const unsigned height,
		int boundary,
		PtrSz<correspondence> correspondence,
		const Intrinsic intrinsic,
		int* correspondence_counter
	) {
		const auto x = threadIdx.x + blockIdx.x * blockDim.x;
		if (x >= boundary) return;
		
		const float3 reference_v = make_float3(src_vertex[x], src_vertex[boundary+x], src_vertex[boundary*2+x]);
		const float3 reference_n = make_float3(src_normal[x], src_normal[boundary+x], src_normal[boundary*2+x]);
		const ushort2 img_coord = {
					__float2uint_rn(((reference_v.x / (reference_v.z + 1e-10)) * intrinsic.focal_x) + intrinsic.principal_x),
					__float2uint_rn(((reference_v.y / (reference_v.z + 1e-10)) * intrinsic.focal_y) + intrinsic.principal_y)
				};
		if(img_coord.x < width && img_coord.y < height)
		{	
			// 当前帧定点图
			const float4 live_v4 = tex2D<float4>(live_vertex_map, img_coord.x, img_coord.y);
			const float4 live_n4 = tex2D<float4>(live_normal_map, img_coord.x, img_coord.y);
			//Check correspondence
			if(dotxyz(reference_n, live_n4) < 0.8f || squared_distance(reference_v, live_v4) > (0.01f * 0.01f) || is_zero_vertex(live_v4)) {
				//Pass
			}
			else {
				const int index = atomicAdd(correspondence_counter, 1);
				correspondence[index].src_vertex = make_float4(reference_v.x, reference_v.y, reference_v.z, 0.0f);
				correspondence[index].src_normal = make_float4(reference_n.x, reference_n.y, reference_n.z, 0.0f);
				correspondence[index].tag_vertex = live_v4;
				correspondence[index].tag_normal = live_n4;

			}
		}
	}
}

void RigidSolver::allocateReduceBuffer() {
    m_reduced_matrix_vector.AllocateBuffer(device::RigidSolverDevice::total_shared_size);
    m_reduced_matrix_vector.ResizeArrayOrException(device::RigidSolverDevice::total_shared_size);

    const auto& pixel_size = m_image_height*m_image_width;
    // 27项数据的global_buffer 共有27*num of blocks个数据 每一行代表27项中的一项数据
    m_reduce_buffer.create(device::RigidSolverDevice::total_shared_size, 
                            divUp(pixel_size,device::RigidSolverDevice::block_size));
	m_correspondences.AllocateBuffer(pixel_size);
}

void RigidSolver::rigidSolverDeviceIteration(cudaStream_t stream) {
    device::RigidSolverDevice solver;

    solver.intrinsic = m_project_intrinsic;
    solver.init_refer2live = m_curr_refer2live;
    solver.image_rows = m_image_height;
    solver.image_cols = m_image_width;
    solver.reference_maps.vertex_map = m_reference_maps.vertex_map;
    solver.reference_maps.normal_map = m_reference_maps.normal_map;
    solver.live_maps.vertex_map = m_live_maps.vertex_map;
    solver.live_maps.normal_map = m_live_maps.normal_map;

    dim3 blk(device::RigidSolverDevice::block_size);
    dim3 grid(divUp(m_image_width*m_image_height, device::RigidSolverDevice::block_size));

    device::rigidSolveIterationKernel<<<grid, blk, 0, stream>>>(solver, m_reduce_buffer);

    device::columnReduceKernel<<<dim3(1,1,1),dim3(32,device::RigidSolverDevice::total_shared_size,1)>>>(
        m_reduce_buffer,
        m_reduced_matrix_vector.DevicePtr()
    );

    m_reduced_matrix_vector.SynchronizeToHost(stream, false);

}

void RigidSolver::FindCorrespondences(DeviceArray<correspondence>& correspondences,
										int& correspondence_count,
										cudaStream_t stream
) {
	int* d_correspondence_count;
	cudaMalloc((void**)&d_correspondence_count, sizeof(int));
	cudaMemset(d_correspondence_count, 0, sizeof(int));
	dim3 block(16, 16);
	dim3 grid(divUp(m_image_width, block.x), divUp(m_image_height, block.y));
	device::FindCorresponceKernel<<<grid, block, 0, stream>>>(
		m_reference_maps.vertex_map,
		m_reference_maps.normal_map,
		m_live_maps.vertex_map,
		m_live_maps.normal_map,
		correspondences,
		m_image_width,
		m_image_height,
		m_project_intrinsic,
		d_correspondence_count,
		m_curr_refer2live
	);

	int h_correspondence_count;
	cudaSafeCall(cudaMemcpy(&h_correspondence_count, d_correspondence_count, sizeof(int), cudaMemcpyDeviceToHost));
	correspondence_count = h_correspondence_count;
	cudaFree(d_correspondence_count);
}

void RigidSolver::FindCorrespondences(DeviceArray<correspondence>& correspondences,
									Eigen::MatrixXd& src_vertex,
									Eigen::MatrixXd& src_normal,
                            		int& correspondence_count,
                            		cudaStream_t stream
) {
	int* d_correspondence_count;
	cudaMalloc((void**)&d_correspondence_count, sizeof(int));
	cudaMemset(d_correspondence_count, 0, sizeof(int));

	double* d_src_vertex;
	double* d_src_normal;
	size_t src_size = sizeof(double) * src_vertex.rows() * src_vertex.cols();
	// for(int i = 0; i < src_vertex.rows(); i++)
	// {
	// 	printf("src_vertex[%d]: %f %f %f\n", i, src_vertex(i, 0), src_vertex(i, 1), src_vertex(i, 2));
	// }
	cudaSafeCall(cudaMalloc((void**)&d_src_vertex, src_size));
	cudaSafeCall(cudaMemcpy(d_src_vertex, src_vertex.data(),src_size, cudaMemcpyHostToDevice));
	cudaSafeCall(cudaMalloc((void**)&d_src_normal, src_size));
	cudaSafeCall(cudaMemcpy(d_src_normal, src_normal.data(),src_size, cudaMemcpyHostToDevice));

	dim3 block(256,1);
	dim3 grid(divUp(src_vertex.rows(), block.x));
	device::FindCorresponceKernel<<<grid, block, 0, stream>>>(
		d_src_vertex,
		d_src_normal,
		m_live_maps.vertex_map,
		m_live_maps.normal_map,
		m_image_width,
		m_image_height,
		src_vertex.rows(),
		correspondences,
		m_project_intrinsic,
		d_correspondence_count
	);

	int h_correspondence_count;
	cudaSafeCall(cudaMemcpy(&h_correspondence_count, d_correspondence_count, sizeof(int), cudaMemcpyDeviceToHost));
	correspondence_count = h_correspondence_count;
	cudaFree(d_correspondence_count);
}
