#pragma once
#include <cuda_runtime.h>
#include <Eigen/Eigen>

using Matrix3f = Eigen::Matrix3f;
using Vector3f = Eigen::Vector3f;
using Matrix4f = Eigen::Matrix4f;
using Vector4f = Eigen::Vector4f;
using Matrix6f = Eigen::Matrix<float, 6, 6>;
using Vector6f = Eigen::Matrix<float, 6, 1>;
using MatrixXf = Eigen::MatrixXf;
using VectorXf = Eigen::VectorXf;
using Isometry3f = Eigen::Isometry3f;

struct vertex_normal_maps {
    cudaTextureObject_t vertex_map;
    cudaTextureObject_t normal_map;
};

// 相机内参
struct Intrinsic
{   
    __host__ __device__ Intrinsic() 
		: principal_x(0), principal_y(0), focal_x(0), focal_y(0) 
        {}
    
    __host__ __device__ Intrinsic(
				const float focal_x_, const float focal_y_,
				const float principal_x_, const float principal_y_
		) : principal_x(principal_x_), principal_y(principal_y_),
			focal_x(focal_x_), focal_y(focal_y_) {}

    //Cast to float4
    __host__ operator float4() {
        return make_float4(principal_x, principal_y, focal_x, focal_y);
    }

    float principal_x, principal_y;
    float focal_x, focal_y;
};

// 纹理和表面结构
struct CudaTextureSurface
{
    cudaArray_t array;
    cudaTextureObject_t texture;
    cudaSurfaceObject_t surface;
};

struct PixelCoordinate {
		unsigned row;
		unsigned col;
		__host__ __device__ PixelCoordinate(): row(0), col(0) {}
		__host__ __device__ PixelCoordinate(const unsigned row_, const unsigned col_) 
		: row(row_), col(col_) {}

		__host__ __device__ const unsigned& x() const { return col; }
		__host__ __device__ const unsigned& y() const { return row; }
		__host__ __device__ unsigned& x() { return col; }
		__host__ __device__ unsigned& y() { return row; }
	};

// 面元数据
struct Surfel {
    PixelCoordinate pixel_coord;
    float4 vertex;
    float4 normal;
    float4 color;
};


struct correspondence {
    // PixelCoordinate source_index;
    // PixelCoordinate target_index;
    float4 src_vertex;
    float4 src_normal;
    float4 tag_vertex;
    float4 tag_normal;
};

