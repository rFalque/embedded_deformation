/*
*   Greedy Search
*   by R. Falque
*   29/11/2018
*/

#ifndef GREEDY_SEARCH
#define GREEDY_SEARCH

#include <Eigen/Core>
#include <vector>

typedef int vertex_t;
typedef double weight_t;
 
struct neighbor {
    vertex_t target;
    weight_t weight;
    neighbor(vertex_t arg_target, weight_t arg_weight)
        : target(arg_target), weight(arg_weight) { }
};

typedef std::vector<std::vector<neighbor> > adjacency_list_t;

class greedy_search
{
public:
	greedy_search(Eigen::MatrixXd V, 
					Eigen::MatrixXi F,
					std::vector<int> indexes_of_sub_V_in_V);
	~greedy_search(){
	}

	std::vector< int > return_k_closest_points(int source, int k);

private:
	int k_;
	int graph_size_;
	adjacency_list_t adjacency_list_;

	Eigen::MatrixXd V_;
	Eigen::MatrixXd sub_V_;
	Eigen::MatrixXi F_;

	std::vector<int> indexes_of_sub_V_in_V_;
};

#endif
