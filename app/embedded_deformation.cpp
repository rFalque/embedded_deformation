/**
 * Author: R. Falque
 * 
 * main for testing the embedded deformation implementation
 * by R. Falque
 * 14/11/2018
 **/

// dependencies

#include "IO/readPLY.h"
#include "IO/writePLY.h"
#include "IO/loadCSV.hpp"



#include "embedded_deformation/embedDeform.hpp"



int main(int argc, char* argv[])
{
    int grid_resolution = 40;
    int k = 6;
    std::string path_pairwise_correspondence = "../data/M.csv";
    std::string path_input_to_deform = "../data/M.ply";
    std::string path_output = "../data/output.ply";


    /* 
     * V: vertex of the surface
     * F: faces of the surface
     * N: nodes of the deformation graph
     * E: edges of the deformation graph
     */
    
    Eigen::MatrixXd V, N;
    Eigen::MatrixXi F, E;
    embedded_deformation* non_rigid_deformation;

    igl::readPLY(path_input_to_deform ,V, F);

    if (F.rows() != 0) // use mesh
        non_rigid_deformation = new embedded_deformation(V, F, grid_resolution, k);
    else /* use knn distance */
        non_rigid_deformation = new embedded_deformation(V, grid_resolution, k);
    
    Eigen::MatrixXd old_points, new_points;
    load_csv(path_pairwise_correspondence, old_points, new_points);

    std::cout << "progress : start deformation ..." << std::endl;
    Eigen::MatrixXd V_deformed;
    non_rigid_deformation->deform(old_points, new_points, V_deformed);

    igl::writePLY(path_output, V_deformed, F);
    return 0;
}
