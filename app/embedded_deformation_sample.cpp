/**
 * Author: R. Falque
 * 
 * main for testing the embedded deformation implementation
 * by R. Falque
 **/

// dependencies


#include "utils/IO_libIGL/readOBJ.h"
#include "utils/IO/readPLY.h"
#include "utils/IO/writePLY.h"
#include "utils/IO/readCSV.h"



#include "libGraphCpp/readGraphOBJ.hpp"
#include "libGraphCpp/polyscopeWrapper.hpp"

#include "embedded_deformation/embedDeform.hpp"
#include "utils/visualization/plotMesh.h"
#include "utils/visualization/plotCloud.h"
#include "embedded_deformation/options.hpp"

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
    EmbeddedDeformation* non_rigid_deformation;

std::cout << "Progress: load file ...";

    readPLY(opts.path_input_file ,V, F);

std::cout << " done.\n";

    // check for error
    if (opts.use_geodesic and F.rows() == 0)
    {
        std::cout << "Config file error: use_geodesic = true, but nor faces were provided." <<std::endl;
        exit(-1);
    }

    if (opts.visualization)
        if (F.rows() != 0)
            plot_mesh(V,F);
        else
            plot_cloud(V);
    
    if (opts.graph_provided) /* graph is provided */
    {
        libgraphcpp::readGraphOBJ(opts.path_graph_obj ,N, E);

        // create deformation object (does not use geodesic distance)
        non_rigid_deformation = new EmbeddedDeformation(V, F, N, E, opts);

    }
    else /* graph not provided */
    {
        if (opts.use_geodesic)
            non_rigid_deformation = new EmbeddedDeformation(V, F, opts);
            //non_rigid_deformation = new embedded_deformation(V, F, opts.grid_resolution, opts.graph_connectivity);
        else /* use knn distance */
            non_rigid_deformation = new EmbeddedDeformation(V, opts);
            //non_rigid_deformation = new embedded_deformation(V, opts.grid_resolution, opts.graph_connectivity);
    }

    if (opts.visualization)
        non_rigid_deformation->show_deformation_graph();
    
    // read correspondences
    Eigen::MatrixXd correspondences = read_csv<Eigen::MatrixXd>(opts.path_pairwise_correspondence);
    Eigen::MatrixXd new_points = correspondences.rightCols(3).transpose();
    Eigen::MatrixXd old_points = correspondences.leftCols(3).transpose();

    std::cout << "progress : start deformation ..." << std::endl;
    Eigen::MatrixXd V_deformed;
    non_rigid_deformation->deform(old_points, new_points, V_deformed);

    if (opts.visualization)
        if (F.rows() != 0)
            plot_mesh(V_deformed,F);
        else
            plot_cloud(V_deformed);

    writePLY(opts.path_output_file, V_deformed, F, true);
    return 0;
}
