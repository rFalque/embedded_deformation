#pragma once

#include "imageProcessor/common_types.hpp"
#include "math/device_mat.h"
#include "containers/SynchronizeArray.h"
#include "containers/device_array.hpp"

using namespace CUDA;

class RigidSolver {
private:
    Intrinsic m_project_intrinsic;
    unsigned m_image_width, m_image_height;

    // 当前帧
    struct {
        cudaTextureObject_t vertex_map;
        cudaTextureObject_t normal_map;
    } m_live_maps;

    // 参考帧 这里是上一帧 按道理应该是从model中获取
    struct {
        cudaTextureObject_t vertex_map;
        cudaTextureObject_t normal_map;
    } m_reference_maps;

    mat34 m_curr_refer2live;


public:
    RigidSolver(const unsigned width,
                const unsigned height,
                const Intrinsic& intrinsic);
    ~RigidSolver();

    void SetInputMaps(
        const vertex_normal_maps& live_maps,
        const vertex_normal_maps& reference_maps,
        const mat34& init_world2camera
    );
    Eigen::Matrix<float, 3, 4> test(
        const vertex_normal_maps& live_maps,
        const vertex_normal_maps& reference_maps);

    mat34 Solve(int max_iters = 3, cudaStream_t stream = 0);

    

private:
    DeviceArray2D<float> m_reduce_buffer;
    SynchronizeArray<float> m_reduced_matrix_vector;
    DeviceBufferArray<correspondence> m_correspondences;

    void allocateReduceBuffer();
    void rigidSolverDeviceIteration(cudaStream_t stream = 0);
    // 找匹配点
    // 图像处理后找匹配点
    void findCorrespondences(cudaStream_t stream = 0);
    void FindCorrespondences(DeviceArray<correspondence>& correspondences,
                            int& correspondence_count,
                            cudaStream_t stream = 0);
    

    

    Eigen::Matrix<float, 6, 6> JtJ_;  // Ci       Cx = d
    Eigen::Matrix<float, 6, 1> JtErr_;  // di
    void rigidSolverHostIterationSync(cudaStream_t stream = 0);

public:
    // 函数接口
    std::vector<correspondence> getCorrespondence() const {
        DeviceArrayView<correspondence> correspondence_view = m_correspondences.ArrayView();
        std::vector<correspondence> correspondences;
        correspondence_view.Download(correspondences);
        return correspondences;
        }

    // 一次非刚性求解后找匹配点
    void findCorrespondences(Eigen::MatrixXd& source_vertex,
                            Eigen::MatrixXd& source_normal,
                            DeviceBufferArray<correspondence>& correspondences,
                            cudaStream_t stream = 0);
    void FindCorrespondences(DeviceArray<correspondence>& correspondences,
                            Eigen::MatrixXd& source_vertex,
                            Eigen::MatrixXd& source_normal,
                            int& correspondence_count,
                            cudaStream_t stream = 0); 
};