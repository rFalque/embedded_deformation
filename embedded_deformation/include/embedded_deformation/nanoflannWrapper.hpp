/*
*   nanoflann wrapper
*   by R. Falque
*   21/11/2018
*/

#ifndef NANOFLANN_WRAPPER
#define NANOFLANN_WRAPPER

#include <iostream>
#include <memory>

#include "nanoflann.hpp"

class nanoflann_wrapper
{
private:
	Eigen::MatrixXd target_;
	std::shared_ptr < nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd > > kd_tree_index;

public:
	nanoflann_wrapper(Eigen::MatrixXd target) : target_(target)
	{
		// set up kdtree
		int leaf_size=10;
		int dimensionality=3;

		this->kd_tree_index = std::make_shared< nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd> >(dimensionality, target_, leaf_size);
		this->kd_tree_index->index->buildIndex();

	}

	~nanoflann_wrapper(){
	}

	std::vector< int > return_k_closest_points(Eigen::Vector3d query_point, int k)
	{


		// Query point:
		std::vector<double> query_pt;
		for (int d=0; d<3; d++)
			query_pt.push_back( query_point(d) );

		// set wtf vectors
		std::vector<size_t> ret_indexes(k);
		std::vector<double> out_dists_sqr(k);
		nanoflann::KNNResultSet<double> resultSet(k);
		resultSet.init( &ret_indexes.at(0), &out_dists_sqr.at(0) );

		// knn search
		this->kd_tree_index->index->findNeighbors(resultSet, &query_pt.at(0), nanoflann::SearchParams(k));

		// pack result into std::vector<int>
		std::vector< int > indexes;
		for (int i = 0; i < k; i++)
			indexes.push_back( ret_indexes.at(i) );

		return indexes;
	}

	bool radius_search(Eigen::Vector3d query_point, double max_dist, std::vector<int> & indexes, std::vector<double> & distances)
	{
		// Query point:
		std::vector<double> query_pt;
		for (int d=0; d<3; d++)
			query_pt.push_back( query_point(d) );

		//search
		nanoflann::SearchParams params;
		std::vector<std::pair<Eigen::Index, double> > matches;
		this->kd_tree_index->index->radiusSearch(&query_pt.at(0), max_dist, matches, params);

		// pack result into std::vector<int>
		for (int i = 0; i < matches.size(); i++)
		{
			indexes.push_back( matches.at(i).first );
			distances.push_back( matches.at(i).second );
		}
		
		return true;
	}

};

#endif

