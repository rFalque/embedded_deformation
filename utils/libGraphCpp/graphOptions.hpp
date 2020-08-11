#ifndef GRAPH_OPTIONS_HPP
#define GRAPH_OPTIONS_HPP

#include <string>
#include <iostream>

struct graphOptions
{
    // path to the graph to load
    std::string path_graph_obj;

    // general options
    bool visualization = false;
    bool verbose = false;

    // visualization parameters
    std::vector<double> nodes_color = {1.0, 0.1, 0.1};
    std::vector<double> edges_color = {0.1, 0.1, 0.1};

    void print(){
        std::cout << "list of the parameters:" << std::endl;
        std::cout << std::endl;
        std::cout << "*** IO files: ***" << std::endl;
        std::cout <<  "path_graph_obj: " << path_graph_obj << std::endl;
        std::cout << std::endl;
        std::cout << "*** General parameters ***" << std::endl;
        std::cout << "visualization: " << visualization << std::endl;
        std::cout << "verbose: " << verbose << std::endl;
        std::cout << std::endl;
        std::cout << "*** Visualization parameters ***" << std::endl;
        std::cout << "nodes_color: [" << nodes_color[0] << "," << nodes_color[1] << "," << nodes_color[2] << "]" << std::endl;
        std::cout << "edges_color: [" << edges_color[0] << "," << edges_color[1] << "," << edges_color[2] << "]" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

};

#endif
