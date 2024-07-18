#include "rigidSolver/rigidSolver.h"

RigidSolver::RigidSolver(const unsigned width,
                        const unsigned height,
                        const Intrinsic& intrinsic
) {
    m_project_intrinsic = intrinsic;
    m_image_width = width;
    m_image_height = height;
    m_curr_refer2live = mat34::identity();

    allocateReduceBuffer();

}

RigidSolver::~RigidSolver(){

}

void RigidSolver::SetInputMaps(
    const vertex_normal_maps& live_maps,
    const vertex_normal_maps& reference_maps,
    const mat34& init_refer2live
) {
    m_live_maps.vertex_map = live_maps.vertex_map;
    m_live_maps.normal_map = live_maps.normal_map;
    m_reference_maps.vertex_map = reference_maps.vertex_map;
    m_reference_maps.normal_map = reference_maps.normal_map;
    m_curr_refer2live = init_refer2live;
}

mat34 RigidSolver::Solve(int max_iters, cudaStream_t stream) {
    for(int i = 0; i< max_iters; i++){
        rigidSolverDeviceIteration(stream);
        rigidSolverHostIterationSync(stream);
    }
    return m_curr_refer2live;
}

void RigidSolver::rigidSolverHostIterationSync(cudaStream_t stream) {
    cudaSafeCall(cudaStreamSynchronize(stream));

    const auto& host_array = m_reduced_matrix_vector.HostArray();

    auto shift = 0;
    for(int i = 0; i < 6; i++){
        for(int j = i; j < 6; j++){
            const float value = host_array[shift++];
            JtJ_(i,j) = value;
            JtJ_(j,i) = value;
        }
    }
    for (int i = 0; i < 6; i++) {
		const float value = host_array[shift++];
		JtErr_[i] = value;
	}

    //Solve it
	Eigen::Matrix<float, 6, 1> x = JtJ_.llt().solve(JtErr_).cast<float>();

	//Update the se3
	const float3 twist_rot = make_float3(x(0), x(1), x(2));
	const float3 twist_trans = make_float3(x(3), x(4), x(5));
	const mat34 se3_update(twist_rot, twist_trans);
	m_curr_refer2live = se3_update * m_curr_refer2live;

}

void RigidSolver::findCorrespondences(cudaStream_t stream) {
    int correspondence_count = 0;
    DeviceArray<correspondence> correspondences = m_correspondences.Array();
    FindCorrespondences(correspondences,
                        correspondence_count);
    m_correspondences.ResizeArrayOrException(correspondence_count);
    std::cout<<"find " << correspondence_count << " correspondences"<<std::endl;

}

void RigidSolver::findCorrespondences(Eigen::MatrixXd& source_vertex,
                                    Eigen::MatrixXd& source_normal,
                                    DeviceBufferArray<correspondence>& correspondences,
                                    cudaStream_t stream) {
    int correspondence_count = 0;
    DeviceArray<correspondence> correspondences_array = correspondences.Array();
    FindCorrespondences(correspondences_array,
                        source_vertex,
                        source_normal,
                        correspondence_count,
                        stream);
    correspondences.ResizeArrayOrException(correspondence_count);
    std::cout<<"find new " << correspondence_count << " correspondences"<<std::endl;
}

Eigen::Matrix<float, 3, 4> RigidSolver::test(
    const vertex_normal_maps& live_maps,
    const vertex_normal_maps& reference_maps
){
    SetInputMaps(live_maps, reference_maps, mat34::identity());
    const mat34 refer2live = Solve();
    Eigen::Matrix<float, 3, 4> init_guess;
    init_guess << refer2live.rot.cols[0].x, refer2live.rot.cols[1].x, refer2live.rot.cols[2].x, refer2live.trans.x,
                refer2live.rot.cols[0].y, refer2live.rot.cols[1].y, refer2live.rot.cols[2].y, refer2live.trans.y,
                refer2live.rot.cols[0].z, refer2live.rot.cols[1].z, refer2live.rot.cols[2].z, refer2live.trans.z;
    // std::cout << "refer2live: " <<std::endl;
    // std::cout << refer2live.rot.cols[0].x << " " 
    //         << refer2live.rot.cols[0].y << " " 
    //         << refer2live.rot.cols[0].z << " " 
    //         << refer2live.trans.x << std::endl;
    // std::cout<< refer2live.rot.cols[1].x << " " 
    //         << refer2live.rot.cols[1].y << " " 
    //         << refer2live.rot.cols[1].z << " "
    //         << refer2live.trans.y << std::endl;
    // std::cout<< refer2live.rot.cols[2].x << " " 
    //         << refer2live.rot.cols[2].y << " " 
    //         << refer2live.rot.cols[2].z << " "
    //         << refer2live.trans.z << std::endl;
    findCorrespondences();
    return init_guess;

}