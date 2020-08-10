/**
 * Author: R. Falque
 * 
 * plot mesh in polyscope
 * by R. Falque
 * 26/09/2019
 **/

#ifndef PLOT_MESH_HPP
#define PLOT_MESH_HPP

#include <Eigen/Core>

#include <string>

#include "polyscopeWrapper.hpp"

namespace visualization {

    inline bool plot_mesh (const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces) {
        visualization::add_mesh(vertices, faces, "mesh");
        visualization::show();
        visualization::clear_structures();
        return true;
    };

    inline bool plot_mesh (const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& color) {
        visualization::add_mesh(vertices, faces, "mesh");
        visualization::add_color_to_mesh(color, "mesh", "color");
        visualization::show();
        visualization::clear_structures();
        return true;
    };

    inline bool plot_cloud (const Eigen::MatrixXd& vertices) {
        visualization::add_cloud(vertices, "cloud");
        visualization::show();
        visualization::clear_structures();
        return true;
    };

    inline bool plot_cloud (const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& color) {
        visualization::add_cloud(vertices, "cloud");
        visualization::add_color_to_cloud(color, "cloud", "color");
        visualization::show();
        visualization::clear_structures();
        return true;
    };

}

#endif
