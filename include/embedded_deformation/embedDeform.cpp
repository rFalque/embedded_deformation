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


embedded_deformation::embedded_deformation(Eigen::MatrixXd V_in, 
										   Eigen::MatrixXi F_in,
										   Eigen::MatrixXd N_in, 
										   Eigen::MatrixXi E_in,
										   options opts)
{
	std::cout << "embedded deformation constructor: graph provided\n";
	use_knn_ = true;
	use_dijkstra_ = false;
	
	// extract point cloud
	nodes_connectivity_ = opts.graph_connectivity;
	verbose_ = opts.verbose;

	// save the data
	V_ = V_in;
	F_ = F_in;
	deformation_graph = new libgraphcpp::Graph(N_in, E_in);
}


embedded_deformation::embedded_deformation(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in,
											double grid_resolution,
											int graph_connectivity)
{
	std::cout << "use geodesic distance to look for closest point\n";

	// set up class variables
	verbose_ = true;
	use_knn_ = false;
	use_dijkstra_ = true;
	nodes_connectivity_ = graph_connectivity;
	V_ = V_in;
	F_ = F_in;

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
	deformation_graph = new libgraphcpp::Graph(N, E);
}

embedded_deformation::embedded_deformation(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in,
										   options opts)
{
	std::cout << "use geodesic distance to look for closest point\n";

	// set up class variables
	verbose_ = true;
	use_knn_ = false;
	use_dijkstra_ = true;
	nodes_connectivity_ = opts.graph_connectivity;
	use_farthest_sampling_ = opts.use_farthest_sampling;
	V_ = V_in;
	F_ = F_in;

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
	deformation_graph = new libgraphcpp::Graph(N, E);
}


embedded_deformation::embedded_deformation(Eigen::MatrixXd V_in,
											double grid_resolution,
											int nodes_connectivity)
{
	std::cout << "use knn to look for closest point\n";

	// set up class variables
	verbose_ = true;
	use_knn_ = true;
	use_dijkstra_ = false;
	nodes_connectivity_ = nodes_connectivity;
	V_ = V_in;

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

	deformation_graph = new libgraphcpp::Graph(N, E);
}


embedded_deformation::embedded_deformation(Eigen::MatrixXd V_in,
										   options opts)
{
	std::cout << "use knn to look for closest point\n";

	// set up class variables
	verbose_ = true;
	use_knn_ = true;
	use_dijkstra_ = false;
	nodes_connectivity_ = opts.graph_connectivity;
	use_farthest_sampling_ = opts.use_farthest_sampling;
	V_ = V_in;

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

	deformation_graph = new libgraphcpp::Graph(N, E);
}


void embedded_deformation::deform(Eigen::MatrixXd sources, Eigen::MatrixXd targets, Eigen::MatrixXd & V_deformed)
{
	std::vector< std::vector<int> > sources_nodes_neighbours;

	// declare the greedy_search object
	greedy_search *search_object;
	nanoflann_wrapper *tree2;

	if (use_dijkstra_)
		search_object = new greedy_search(V_, F_, indexes_of_deformation_graph_in_V_);
	if (use_knn_)
		tree2 = new nanoflann_wrapper(deformation_graph->get_nodes());

	// set the sources neighbours
    sources_nodes_neighbours.resize(sources.rows());
	for (int i = 0; i < sources.rows(); ++i)
	{
		// find closest point on the mesh
		nanoflann_wrapper tree(V_);
		std::vector< int > closest_point;
		closest_point = tree.return_k_closest_points(sources.row(i), 1);

		if (use_dijkstra_)
			sources_nodes_neighbours.at(i) = search_object->return_k_closest_points(closest_point.at(0), nodes_connectivity_+1);
		if (use_knn_)
			sources_nodes_neighbours.at(i) = tree2->return_k_closest_points(V_.row(closest_point.at(0)), nodes_connectivity_+1);
	}

	// initialize the parameters that needs to be optimized (nodes rotations and translations)
	std::vector< std::vector <double> > params;

	params.resize(deformation_graph->num_nodes());
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
		
	for (int i = 0; i < deformation_graph->num_nodes(); ++i)
	{
		//RotCostFunction* cost_function = new RotCostFunction(w_rot_);
		RotCostFunction_2* cost_function = new RotCostFunction_2(w_rot_);
		problem.AddResidualBlock(cost_function, NULL, &params[i][0]);
	}
		
	if (verbose_)
		std::cout << "progress : add the residuals for RegCostFunction ..." << std::endl;

	for (int i = 0; i < deformation_graph->num_nodes(); ++i)
	{
		for (int j = 0; j < deformation_graph->get_adjacency_list(i).size(); ++j)
		{
			int k = deformation_graph->get_adjacency_list(i, j);
			Eigen::Vector3d g_j = deformation_graph->get_node(i);
			Eigen::Vector3d g_k = deformation_graph->get_node(k);
			RegCostFunction* cost_function = new RegCostFunction(w_reg_, g_j, g_k);

			// add residual block
			problem.AddResidualBlock(cost_function, NULL, &params[i][0], &params[k][0]);
		}
	}
	
	if (verbose_)
		std::cout << "progress : add the residuals for RigCostFunction ..." << std::endl;

	std::vector<int> bridges;
	deformation_graph->has_bridges(bridges);
	for (int i=0; i<bridges.size(); i++)
	{
		Eigen::Vector2i nodes_list = deformation_graph->get_edge(bridges[i]);
		Eigen::Vector3d g_j = deformation_graph->get_node(nodes_list(0));
		Eigen::Vector3d g_k = deformation_graph->get_node(nodes_list(1));
		RigCostFunction* cost_function = new RigCostFunction(w_rig_, g_j, g_k);

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
			vector_g.push_back( deformation_graph->get_node(sources_nodes_neighbours[i][j]) );

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
	rotation_matrices_.resize(deformation_graph->num_nodes());
	translations_.resize(deformation_graph->num_nodes());
	for (int i = 0; i < deformation_graph->num_nodes(); ++i)
	{
		rotation_matrices_[i] = Eigen::Map<Eigen::Matrix3d>( &(params[i][0]) );
		translations_[i] = Eigen::Map<Eigen::Vector3d> ( &(params[i][0]) + 9);
	}

	// redefine the deformed mesh
	V_deformed = V_;

	std::vector<double> w_j;
	Eigen::Vector3d new_point;
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
			w_j.push_back( pow(1 - (V_.row(i) - deformation_graph->get_node(   neighbours_nodes[j]   ).transpose() ).squaredNorm() 
				                 / (V_.row(i) - deformation_graph->get_node( neighbours_nodes.back() ).transpose() ).squaredNorm(), 2) );

		double normalization_factor = 0;
		for (int j = 0; j < nodes_connectivity_; ++j)
			normalization_factor += w_j[j];

		for (int j = 0; j < nodes_connectivity_; ++j)
			w_j[j] /= normalization_factor;

		new_point << 0, 0, 0;
		for (int j = 0; j < nodes_connectivity_; ++j)
			new_point += w_j[j] * (rotation_matrices_[ neighbours_nodes[j] ] * (V_.row(i).transpose() - deformation_graph->get_node( neighbours_nodes[j] ) ) 
				                                                           + deformation_graph->get_node( neighbours_nodes[j] ) + translations_[ neighbours_nodes[j] ]);

		V_deformed.row(i) = new_point;

	}
}


void embedded_deformation::show_deformation_graph()
{
	visualization::plot(*deformation_graph, "deformation graph");
}
