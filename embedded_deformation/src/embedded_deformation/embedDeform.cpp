/*
*   embedded deformation implementation
*   by R. Falque
*   14/11/2018
*/



/* TODO list:
 *  - add the distance in the creation of the graph
 *  - unidirectional versus bidirectional
 *  - what about the point N+1 (last one used for the weights)
 *  - 
 */


// dependencies
#include "embedDeform.hpp"

#include "nanoflannWrapper.hpp"
#include "greedySearch.hpp"
#include "downsampling.hpp"

//#include "costFunction.hpp"

#include <vector>
#include <chrono>
#include <cmath>

// Check gradient flag
bool const kCheckGradient = false;
bool const kCheckGradientReport = false;


EmbeddedDeformation::EmbeddedDeformation(Eigen::MatrixXd V_in, 
										   Eigen::MatrixXi F_in,
										   Eigen::MatrixXd N_in, 
										   Eigen::MatrixXi E_in,
										   options opts)
{
	std::cout << "embedded deformation constructor: graph provided\n";

	// check the data
	if (V_in.rows() == 3 && F_in.rows() == 3 && N_in.rows() == 3 && E_in.rows() == 2 ) {
		V_ = V_in.transpose();
		F_ = F_in.transpose();
		deformation_graph_ptr_ = new libgraphcpp::Graph(N_in.transpose(), E_in.transpose());
		transpose_input_and_output_ = true;
	} else if (V_in.cols() == 3 && F_in.cols() == 3 && N_in.cols() == 3 && E_in.cols() == 2 ) {		
		V_ = V_in;
		F_ = F_in;
		deformation_graph_ptr_ = new libgraphcpp::Graph(N_in, E_in);
		transpose_input_and_output_ = false;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}
	
	use_knn_ = true;
	use_dijkstra_ = false;
	
	// extract point cloud
	nodes_connectivity_ = opts.graph_connectivity;
	verbose_ = opts.verbose;

}


EmbeddedDeformation::EmbeddedDeformation(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in,
											double grid_resolution,
											int graph_connectivity)
{
	std::cout << "use geodesic distance to look for closest point\n";

	// save the data
	if (V_in.rows() == 3 && F_in.rows() == 3) {
		V_ = V_in.transpose();
		F_ = F_in.transpose();
		transpose_input_and_output_ = true;
	} else if (V_in.cols() == 3 && F_in.cols() == 3) {		
		V_ = V_in;
		F_ = F_in;
		transpose_input_and_output_ = false;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}

	// set up class variables
	verbose_ = true;
	use_knn_ = false;
	use_dijkstra_ = true;
	nodes_connectivity_ = graph_connectivity;

	// define nodes as subset
	Eigen::MatrixXd N;
	//downsampling(V_, N, indexes_of_deformation_graph_in_V_, grid_resolution, 0, false, true);
	downsampling(V_, N, indexes_of_deformation_graph_in_V_, grid_resolution, 0, use_farthest_sampling_, true);

	// define edges
	Eigen::MatrixXi E(N.rows()*(nodes_connectivity_+1), 2);
	greedy_search search_object(V_, F_, indexes_of_deformation_graph_in_V_);
    std::vector<int> closest_points;
	int counter = 0;
    for (int i = 0; i < N.rows(); ++i) {
		closest_points = search_object.return_k_closest_points( indexes_of_deformation_graph_in_V_[i], nodes_connectivity_+2 );
		closest_points.erase( closest_points.begin() );

		for (int j=0; j< closest_points.size(); j++) {
			E(counter, 0) = i;
			E(counter, 1) = closest_points[j];
			counter ++;
		}
    }
	
	// define deformation graph
	deformation_graph_ptr_ = new libgraphcpp::Graph(N, E);
}

EmbeddedDeformation::EmbeddedDeformation(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in,
										   options opts)
{
	std::cout << "use geodesic distance to look for closest point\n";

	// save the data
	if (V_in.rows() == 3 && F_in.rows() == 3) {
		V_ = V_in.transpose();
		F_ = F_in.transpose();
		transpose_input_and_output_ = true;
	} else if (V_in.cols() == 3 && F_in.cols() == 3) {		
		V_ = V_in;
		F_ = F_in;
		transpose_input_and_output_ = false;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}

	// set up class variables
	verbose_ = true;
	use_knn_ = false;
	use_dijkstra_ = true;
	nodes_connectivity_ = opts.graph_connectivity;
	use_farthest_sampling_ = opts.use_farthest_sampling;

	w_rot_ = opts.w_rot;
	w_reg_ = opts.w_reg;
	w_rig_ = opts.w_rig;
	w_con_ = opts.w_con;

	// define nodes as subset
	Eigen::MatrixXd N;
	//downsampling(V_, N, indexes_of_deformation_graph_in_V_, opts.grid_resolution, use_farthest_sampling_);
	downsampling(V_, N, indexes_of_deformation_graph_in_V_, 
	             opts.grid_resolution, opts.grid_size, 
				 opts.use_farthest_sampling, opts.use_relative);

	// define edges
	Eigen::MatrixXi E(N.rows()*(nodes_connectivity_+1), 2);
	greedy_search search_object(V_, F_, indexes_of_deformation_graph_in_V_);
    std::vector<int> closest_points;
	int counter = 0;
    for (int i = 0; i < N.rows(); ++i) {
		closest_points = search_object.return_k_closest_points( indexes_of_deformation_graph_in_V_[i], nodes_connectivity_+2 );
		closest_points.erase( closest_points.begin() );

		for (int j=0; j< closest_points.size(); j++) {
			E(counter, 0) = i;
			E(counter, 1) = closest_points[j];
			counter ++;
		}
    }
	
	// define deformation graph
	deformation_graph_ptr_ = new libgraphcpp::Graph(N, E);
}


EmbeddedDeformation::EmbeddedDeformation(Eigen::MatrixXd V_in,
											double grid_resolution,
											int nodes_connectivity)
{
	std::cout << "use knn to look for closest point\n";

	// save the data
	if (V_in.rows() == 3) {
		V_ = V_in.transpose();
		transpose_input_and_output_ = true;
	} else if (V_in.cols() == 3) {		
		V_ = V_in;
		transpose_input_and_output_ = false;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}

	// set up class variables
	verbose_ = true;
	use_knn_ = true;
	use_dijkstra_ = false;
	nodes_connectivity_ = nodes_connectivity;

	// extract point cloud
	Eigen::MatrixXd N;
	//downsampling(V_, N, indexes_of_deformation_graph_in_V_, grid_resolution, use_farthest_sampling_);
	downsampling(V_, N, indexes_of_deformation_graph_in_V_, grid_resolution, 0, use_farthest_sampling_, true);

	// build edges:
	Eigen::MatrixXi E(N.rows()*(nodes_connectivity_+1), 2);
	nanoflann_wrapper tree_2(N);
	int counter = 0;
	for (int i = 0; i < N.rows(); ++i) {
		std::vector< int > closest_points;
		closest_points = tree_2.return_k_closest_points(N.row(i), nodes_connectivity_+ 2 );

		for (int j = 0; j < closest_points.size(); ++j)
			if (i != closest_points[j]) {
				E(counter, 0) = i;
				E(counter, 1) = closest_points[j];
				counter ++;
			}
	}

	deformation_graph_ptr_ = new libgraphcpp::Graph(N, E);
}


EmbeddedDeformation::EmbeddedDeformation(Eigen::MatrixXd V_in,
										   options opts)
{
	std::cout << "use knn to look for closest point\n";

	// save the data
	if (V_in.rows() == 3) {
		V_ = V_in.transpose();
		transpose_input_and_output_ = true;
	} else if (V_in.cols() == 3) {		
		V_ = V_in;
		transpose_input_and_output_ = false;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}

	// set up class variables
	verbose_ = true;
	use_knn_ = true;
	use_dijkstra_ = false;
	nodes_connectivity_ = opts.graph_connectivity;
	use_farthest_sampling_ = opts.use_farthest_sampling;

	w_rot_ = opts.w_rot;
	w_reg_ = opts.w_reg;
	w_rig_ = opts.w_rig;
	w_con_ = opts.w_con;

	// extract point cloud
	Eigen::MatrixXd N;
	downsampling(V_, N, indexes_of_deformation_graph_in_V_, 
	             opts.grid_resolution, opts.grid_size, 
				 opts.use_farthest_sampling, opts.use_relative);

	// build edges:
	Eigen::MatrixXi E(N.rows()*(nodes_connectivity_+1), 2);
	nanoflann_wrapper tree_2(N);
	int counter = 0;
	for (int i = 0; i < N.rows(); ++i) {
		std::vector< int > closest_points;
		closest_points = tree_2.return_k_closest_points(N.row(i), nodes_connectivity_+ 2 );

		for (int j = 0; j < closest_points.size(); ++j)
			if (i != closest_points[j]) {
				E(counter, 0) = i;
				E(counter, 1) = closest_points[j];
				counter ++;
			}
	}

	deformation_graph_ptr_ = new libgraphcpp::Graph(N, E);
}





void EmbeddedDeformation::deform(Eigen::MatrixXd sources_in, Eigen::MatrixXd targets_in, Eigen::MatrixXd & V_deformed)
{
	Eigen::MatrixXd sources, targets;

	if (sources_in.cols() == 3 && targets_in.cols() == 3) {		
		sources = sources_in;
		targets = targets_in;
	} else if (sources_in.rows() == 3 && targets_in.rows() == 3) {
		sources = sources_in.transpose();
		targets = targets_in.transpose();
	} else {
		throw std::invalid_argument( "wrong input size" );
	}

	if (sources_in.cols() == 3 && targets_in.cols() == 3 && sources_in.rows() == 3 && targets_in.rows() == 3)
		throw std::invalid_argument( "the system does not work for 3 input points" );

	std::vector< std::vector<int> > sources_nodes_neighbours;

	// declare the greedy_search object
	greedy_search *deformation_graph_greedy_search = nullptr;
	nanoflann_wrapper *deformation_graph_kdtree = nullptr;

	if (use_dijkstra_)
		deformation_graph_greedy_search = new greedy_search(V_, F_, indexes_of_deformation_graph_in_V_);
	if (use_knn_)
		deformation_graph_kdtree = new nanoflann_wrapper(deformation_graph_ptr_->get_nodes());

	// set the sources neighbours
    sources_nodes_neighbours.resize(sources.rows());
	nanoflann_wrapper tree(V_);
	for (int i = 0; i < sources.rows(); ++i)
	{
		// find closest point on the mesh
		std::vector< int > closest_point;
		closest_point = tree.return_k_closest_points(sources.row(i), 1);

		if (use_dijkstra_)
			sources_nodes_neighbours.at(i) = deformation_graph_greedy_search->return_k_closest_points(closest_point.at(0), nodes_connectivity_+1);
		if (use_knn_)
			sources_nodes_neighbours.at(i) = deformation_graph_kdtree->return_k_closest_points(V_.row(closest_point.at(0)), nodes_connectivity_+1);
	}

	// initialize the parameters that needs to be optimized (nodes rotations and translations)
	std::vector< std::vector <double> > params;

	params.resize(deformation_graph_ptr_->num_nodes());
	for (int i = 0; i < params.size(); ++i)
	{
		params[i].resize(12);

		params[i][0 ] = 1;
		params[i][1 ] = 0;
		params[i][2 ] = 0;
		params[i][3 ] = 0;
		params[i][4 ] = 1;
		params[i][5 ] = 0;
		params[i][6 ] = 0;
		params[i][7 ] = 0;
		params[i][8 ] = 1;
		params[i][9 ] = 0;
		params[i][10] = 0;
		params[i][11] = 0;
	}

	// run levenberg marquardt optimization
	ceres::Problem problem;

	if (verbose_)
		std::cout << "progress : add the residuals for RotCostFunction ..." << std::endl;
		
	for (int i = 0; i < deformation_graph_ptr_->num_nodes(); ++i)
	{
		//RotCostFunction* cost_function = new RotCostFunction(w_rot_);
		RotCostFunction_2* cost_function = new RotCostFunction_2(w_rot_);
		problem.AddResidualBlock(cost_function, NULL, &params[i][0]);
	}
		
	if (verbose_)
		std::cout << "progress : add the residuals for RegCostFunction ..." << std::endl;

	for (int i = 0; i < deformation_graph_ptr_->num_nodes(); ++i)
	{
		for (int j = 0; j < deformation_graph_ptr_->get_adjacency_list(i).size(); ++j)
		{
			int k = deformation_graph_ptr_->get_adjacency_list(i, j);
			Eigen::Vector3d g_j = deformation_graph_ptr_->get_node(i);
			Eigen::Vector3d g_k = deformation_graph_ptr_->get_node(k);
			RegCostFunction* cost_function = new RegCostFunction(w_reg_, g_j, g_k);

			// add residual block
			problem.AddResidualBlock(cost_function, NULL, &params[i][0], &params[k][0]);
		}
	}
	
	if (verbose_)
		std::cout << "progress : add the residuals for RigCostFunction ..." << std::endl;

	std::vector<int> bridges;
	deformation_graph_ptr_->has_bridges(bridges);
	for (int i=0; i<bridges.size(); i++)
	{
		Eigen::Vector2i nodes_list = deformation_graph_ptr_->get_edge(bridges[i]);
		RigCostFunction* cost_function = new RigCostFunction(w_rig_);

		// add residual block
		problem.AddResidualBlock(cost_function, NULL, &params[nodes_list(0)][0], &params[nodes_list(1)][0]);
	}

	if (verbose_)
	std::cout << "progress : add the residuals for ConCostFunction ..." << std::endl;
	for (int i = 0; i < sources.rows(); ++i)
	{
		// stack the parameters
		std::vector<double *> parameter_blocks;
		std::vector<Eigen::Vector3d> vector_g;

		// define cost function and associated parameters
		for (int j = 0; j < nodes_connectivity_; ++j)
			parameter_blocks.push_back(&params[ sources_nodes_neighbours[i][j] ][0]);

		for (int j = 0; j < sources_nodes_neighbours[i].size(); ++j)
			vector_g.push_back( deformation_graph_ptr_->get_node(sources_nodes_neighbours[i][j]) );

		// create the cost function
		ConCostFunction* cost_function = new ConCostFunction(w_con_, vector_g, sources.row(i), targets.row(i));

		// add residual block
		problem.AddResidualBlock(cost_function, NULL, parameter_blocks);
	}

	// Run the solver
	ceres::Solver::Options options;
	options.check_gradients = false;
	options.gradient_check_relative_precision = 0.01;
	options.minimizer_progress_to_stdout = verbose_;
	options.num_threads = 1;
	options.max_num_iterations = 100;
	options.function_tolerance = 1e-10;
	options.parameter_tolerance = 1e-10;
	ceres::Solver::Summary summary;
	Solve(options, &problem, &summary);

	if (verbose_)
		std::cout << summary.FullReport() << "\n";

	// Extract the optimized parameters
	rotation_matrices_.resize(deformation_graph_ptr_->num_nodes());
	translations_.resize(deformation_graph_ptr_->num_nodes());
	for (int i = 0; i < deformation_graph_ptr_->num_nodes(); ++i)
	{
		rotation_matrices_[i] = Eigen::Map<Eigen::Matrix3d>( &(params[i][0]) );
		translations_[i] = Eigen::Map<Eigen::Vector3d> ( &(params[i][0]) + 9);
	}

	// redefine the deformed mesh
	Eigen::MatrixXd V_temp = V_;

	std::vector<double> w_j;
	Eigen::Vector3d new_point;
	for (int i = 0; i < V_.rows(); ++i)
	{
		std::vector< int > neighbours_nodes;

		// get closest nodes
		if (use_dijkstra_)
			neighbours_nodes = deformation_graph_greedy_search->return_k_closest_points(i, nodes_connectivity_+1);
		if (use_knn_)
			neighbours_nodes = deformation_graph_kdtree->return_k_closest_points(V_.row(i), nodes_connectivity_+1);

		// equation (3)
		w_j.clear();
		for (int j = 0; j < nodes_connectivity_; ++j)
			w_j.push_back( pow(1 - (V_.row(i) - deformation_graph_ptr_->get_node(   neighbours_nodes[j]   ).transpose() ).squaredNorm() 
				                 / (V_.row(i) - deformation_graph_ptr_->get_node( neighbours_nodes.back() ).transpose() ).squaredNorm(), 2) );

		double normalization_factor = 0;
		for (int j = 0; j < nodes_connectivity_; ++j)
			normalization_factor += w_j[j];

		for (int j = 0; j < nodes_connectivity_; ++j)
			w_j[j] /= normalization_factor;

		new_point << 0, 0, 0;
		for (int j = 0; j < nodes_connectivity_; ++j)
			new_point += w_j[j] * (rotation_matrices_[ neighbours_nodes[j] ] * (V_.row(i).transpose() - deformation_graph_ptr_->get_node( neighbours_nodes[j] ) ) 
				                                                           + deformation_graph_ptr_->get_node( neighbours_nodes[j] ) + translations_[ neighbours_nodes[j] ]);

		V_temp.row(i) = new_point;

	}

	if (transpose_input_and_output_)
		V_deformed = V_temp.transpose();
	else
		V_deformed = V_temp;

	if (deformation_graph_greedy_search != nullptr) delete deformation_graph_greedy_search;
    if (deformation_graph_kdtree != nullptr) delete deformation_graph_kdtree;
}

// Equation (3) in Sumner et al
void EmbeddedDeformation::update_normals(Eigen::MatrixXd & normals) {
	
	// check inputs
	Eigen::MatrixXd N_input;
	transpose_input_and_output_ = false;

	// test the input size
	if (normals.rows() == 3) {
		N_input = normals.transpose();
		transpose_input_and_output_ = true;
	} else if (normals.cols() == 3) {		
		N_input = normals;
	} else {
		throw std::invalid_argument( "wrong input size" );
	}
	if ( N_input.rows() != V_.rows() || N_input.cols() != V_.cols() ) {
		throw std::invalid_argument( "N is supposed to contain the prior normals of the pointcloud (defined per vertex)" );
	}

	// declare the greedy_search object
	greedy_search *search_object;
	nanoflann_wrapper *tree2;

	if (use_dijkstra_)
		search_object = new greedy_search(V_, F_, indexes_of_deformation_graph_in_V_);
	if (use_knn_)
		tree2 = new nanoflann_wrapper(deformation_graph_ptr_->get_nodes());

	std::vector<double> w_j;
	for (int i = 0; i < V_.rows(); ++i)
	{
		std::vector< int > neighbours_nodes;

		// get closest nodes
		if (use_dijkstra_)
			neighbours_nodes = search_object->return_k_closest_points(i, nodes_connectivity_+1);
		if (use_knn_)
			neighbours_nodes = tree2->return_k_closest_points(V_.row(i), nodes_connectivity_+1);

		// equation (3)
		w_j.clear();
		for (int j = 0; j < nodes_connectivity_; ++j)
			w_j.push_back( pow(1 - (V_.row(i) - deformation_graph_ptr_->get_node(   neighbours_nodes[j]   ).transpose() ).squaredNorm() 
				                 / (V_.row(i) - deformation_graph_ptr_->get_node( neighbours_nodes.back() ).transpose() ).squaredNorm(), 2) );

		double normalization_factor = 0;
		for (int j = 0; j < nodes_connectivity_; ++j)
			normalization_factor += w_j[j];

		for (int j = 0; j < nodes_connectivity_; ++j)
			w_j[j] /= normalization_factor;

		Eigen::Vector3d normal_temp = {0, 0, 0};
		for (int j = 0; j < nodes_connectivity_; ++j)
			normal_temp += w_j[j] * rotation_matrices_[ neighbours_nodes[j] ] * N_input.row(i).transpose();

		N_input.row(i) = normal_temp;

	}

	if (transpose_input_and_output_) {
		normals = N_input.transpose();
	} else {
		normals = N_input;
	}
}

void EmbeddedDeformation::show_deformation_graph()
{
	visualization::plot(*deformation_graph_ptr_, "deformation graph");
}
