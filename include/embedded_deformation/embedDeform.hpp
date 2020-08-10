/*
*   embedded deformation implementation
*   by R. Falque
*   14/11/2018
*/

#ifndef EMBEDDED_DEFORMATION
#define EMBEDDED_DEFORMATION

#include <ceres/ceres.h>
#include <ceres/gradient_checker.h>

#include "../utils/options.hpp"

#include "costFunction.hpp"

#include "libGraphCpp/include/libGraphCpp/graph.hpp"
#include "libGraphCpp/include/libGraphCpp/plotGraph.hpp"

class embedded_deformation
{
public:

	// provide the mesh and the graph
	embedded_deformation(Eigen::MatrixXd V_in, 
						 Eigen::MatrixXi F_in,
						 Eigen::MatrixXd N_in, 
						 Eigen::MatrixXi E_in,
						 options opts);

	// provide the mesh
	embedded_deformation(Eigen::MatrixXd V_in, 
						 Eigen::MatrixXi F_in,
						 double grid_resolution,
						 int k);

	// provide only a point cloud
	embedded_deformation(Eigen::MatrixXd V_in, 
						 double grid_resolution,
						 int k);
	~embedded_deformation(){
	}

	void deform(Eigen::MatrixXd sources, Eigen::MatrixXd targets, Eigen::MatrixXd & V_deformed);
	void show_deformation_graph();

private:
	// by order of appearance:
	const double w_rot_ = 1;
	const double w_reg_ = 10;
	const double w_rig_ = 10;
	const double w_con_ = 100;

	// options
	bool use_knn_;
	bool use_dijkstra_;
	bool verbose_;
	int nodes_connectivity_;

	Eigen::MatrixXd V_;
	Eigen::MatrixXi F_;
	
	libgraphcpp::Graph* deformation_graph;

	// graph properties definition
	std::vector<Eigen::Matrix3d> rotation_matrices_;
	std::vector<Eigen::Vector3d> translations_;

	// specific for searching on the geodesic distance
	std::vector<int> indexes_of_deformation_graph_in_V_;
};

#endif
