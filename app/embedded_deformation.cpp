/**
 * Author: R. Falque
 * 
 * main for testing the embedded deformation implementation
 * by R. Falque
 * 14/11/2018
 **/

// dependencies


#include "IO/readOBJ.h"
#include "IO/readPLY.h"
#include "IO/writePLY.h"


#include "libGraphCpp/include/libGraphCpp/readGraphOBJ.hpp"
#include "libGraphCpp/include/libGraphCpp/polyscopeWrapper.hpp"

#include "embedded_deformation/embedDeform.hpp"
#include "plotMesh.hpp"
#include "loadCSV.hpp"
#include "options.hpp"

#include <yaml-cpp/yaml.h>
#include "polyscope/polyscope.h"


int main(int argc, char* argv[])
{
    options opts;
    opts.loadYAML("../config.yaml");
std::cout << "Progress: yaml loaded\n";

    polyscope::init();
std::cout << "Progress: plyscope initialized\n";

    /* 
     * V: vertex of the surface
     * F: faces of the surface
     * N: nodes of the deformation graph
     * E: edges of the deformation graph
     */
    
    Eigen::MatrixXd V, N;
    Eigen::MatrixXi F, E;
    embedded_deformation* non_rigid_deformation;

std::cout << "Progress: load file ...";

    igl::readPLY(opts.path_input_file ,V, F);

std::cout << " done.\n";

    // check for error
    if (opts.use_geodesic and F.rows() == 0)
    {
        std::cout << "Config file error: use_geodesic = true, but nor faces were provided." <<std::endl;
        exit(-1);
    }

    if (opts.visualization)
        if (F.rows() != 0)
            visualization::plot_mesh(V,F);
        else
            visualization::plot_cloud(V);
    
    if (opts.graph_provided) /* graph is provided */
    {
        libgraphcpp::readGraphOBJ(opts.path_graph_obj ,N, E);

        // create deformation object (does not use geodesic distance)
        non_rigid_deformation = new embedded_deformation(V, F, N, E, opts);

    }
    else /* graph not provided */
    {
        if (opts.use_geodesic)
            non_rigid_deformation = new embedded_deformation(V, F, opts);
            //non_rigid_deformation = new embedded_deformation(V, F, opts.grid_resolution, opts.graph_connectivity);
        else /* use knn distance */
            non_rigid_deformation = new embedded_deformation(V, opts);
            //non_rigid_deformation = new embedded_deformation(V, opts.grid_resolution, opts.graph_connectivity);
    }

    if (opts.visualization)
        non_rigid_deformation->show_deformation_graph();
    
    Eigen::MatrixXd old_points, new_points;
    load_csv(opts.path_pairwise_correspondence, old_points, new_points);

    std::cout << "progress : start deformation ..." << std::endl;
    Eigen::MatrixXd V_deformed;
    non_rigid_deformation->deform(old_points, new_points, V_deformed);

    if (opts.visualization)
        if (F.rows() != 0)
            visualization::plot_mesh(V_deformed,F);
        else
            visualization::plot_cloud(V_deformed);

    igl::writePLY(opts.path_output_file, V_deformed, F);
    return 0;
}
