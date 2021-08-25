#pragma once
// Add a backcheck on the correspondence


#include <vector>
#include <string>
#include <Eigen/Dense>
#include "nanoflannWrapper.hpp"
#include "utils/visualization/progressbar.h"

Eigen::MatrixXd sub_space(const Eigen::MatrixXd &cloud, std::vector<int> indices) {
    Eigen::MatrixXd out(3, indices.size());
    for (int i=0; i<indices.size(); i++)
        out.col(i) = cloud.col(indices.at(i));
    return out;
}

template <typename T>
T median_filter(std::vector<T> vector)
{
    std::sort(vector.begin(), vector.end());
    return vector[ std::round( vector.size()/2 ) ];
}


Eigen::Vector3d get_median_point( std::vector<int> pointIdxNKNSearch, Eigen::MatrixXd & cloud )
{
    int x=0, y=1, z=2;
    Eigen::Vector3d median_point;
    std::vector< double > x_vector, y_vector, z_vector;
    for (int i=0; i<pointIdxNKNSearch.size(); ++i)
    {
        x_vector.push_back( cloud(pointIdxNKNSearch[i], x) );
        y_vector.push_back( cloud(pointIdxNKNSearch[i], y) );
        z_vector.push_back( cloud(pointIdxNKNSearch[i], z) );
    }
    median_point(x) = median_filter(x_vector);
    median_point(y) = median_filter(y_vector);
    median_point(z) = median_filter(z_vector);
    return median_point;
}


Eigen::Vector3d point_association(std::string method, 
                                Eigen::Vector3d p_i, 
                                Eigen::Vector3d n_i, 
	                            Eigen::MatrixXd & P,
	                            Eigen::MatrixXd & N)
{
    // used for readability
    int x=0, y=1, z=2;
	Eigen::VectorXd points_dist(P.cols());
	Eigen::VectorXd normals_dist(P.cols());
	for (int i = 0; i < P.cols(); ++i)
	{
		points_dist(i) = sqrt( (p_i(x) - P.col(i)(x)) * (p_i(x) - P.col(i)(x))
			                  +(p_i(y) - P.col(i)(y)) * (p_i(y) - P.col(i)(y))
			                  +(p_i(z) - P.col(i)(z)) * (p_i(z) - P.col(i)(z)));
		normals_dist(i) = n_i(x)*N.col(i)(x) 
		                + n_i(y)*N.col(i)(y) 
		                + n_i(z)*N.col(i)(z);
	}
	Eigen::MatrixXd::Index pointIndex;
	if (method == "closest_point")
	{
		points_dist.minCoeff(&pointIndex);
	}
	if (method == "normals_distances")
	{
		normals_dist.maxCoeff(&pointIndex);
	}
	if (method == "weighted")
	{
		( points_dist/0.1 - normals_dist ).minCoeff(&pointIndex);
	}
	if (method == "hybrid")
	{
		//pointIndex = ( points_dist.array().pow(31)/normals_dist ).maxCoeff(&pointIndex);
		pointIndex = 1;
		//pointIndex = ( pow(points_dist, 31)/normals_dist ).maxCoeff(&pointIndex);
	}

    return P.row(pointIndex);
}


int find_correspondence(std::string method, 
                        Eigen::Vector3d v, 
                        Eigen::Vector3d n, 
	                    const Eigen::MatrixXd & cloud,
	                    const Eigen::MatrixXd & normals,
                        nanoflann_wrapper kd_tree,
                        int K)
{
    // downsampling of the cloud and normals with K closest points
    std::vector<int> closest_point_indices(K);
    closest_point_indices = kd_tree.return_k_closest_points(v, K);
    Eigen::MatrixXd sample_cloud = sub_space(cloud, closest_point_indices);
    Eigen::MatrixXd sample_normals = sub_space(normals, closest_point_indices);

    // build the distance between points and normals
    int x=0, y=1, z=2;
    Eigen::VectorXd points_distance(cloud.cols());
	Eigen::VectorXd normals_distance(cloud.cols());
    for (int i = 0; i < cloud.cols(); ++i)
	{
		points_distance(i) = sqrt( (v(x) - cloud.col(i)(x)) * (v(x) - cloud.col(i)(x))
			                      +(v(y) - cloud.col(i)(y)) * (v(y) - cloud.col(i)(y))
			                      +(v(z) - cloud.col(i)(z)) * (v(z) - cloud.col(i)(z)));
		normals_distance(i) = n(x)*normals.col(i)(x) 
		                    + n(y)*normals.col(i)(y) 
		                    + n(z)*normals.col(i)(z);
	}

    // return the min distance according to the selected method
	Eigen::MatrixXd::Index pointIndex;
	if (method == "closest_point")
	{
		points_distance.minCoeff(&pointIndex);
	}
	if (method == "normals_distances")
	{
		normals_distance.maxCoeff(&pointIndex);
	}
	if (method == "weighted")
	{
		( points_distance/0.1 - normals_distance ).minCoeff(&pointIndex);
	}

    return pointIndex;
}



void get_surface_association(Eigen::MatrixXd V_source, Eigen::MatrixXd N_source,
                             Eigen::MatrixXd V_target, Eigen::MatrixXd N_target,
                             Eigen::MatrixXd &source_position,
                             Eigen::MatrixXd &target_position)
{
    // this need to be moved out of hardcoded parameters into config file
    int x=0, y=1, z=2;                                  // used to access the points
    double distance_threshold = 0.05;                   // used to check that the correspondence of the correspondence is not too far
    int skip_points = 5;                               // pseudo downsampling
    int K = 50;                                         // limit the search space
    
    std::vector<Eigen::Vector3d> correspondences_on_source;
    std::vector<Eigen::Vector3d> correspondences_on_target;
    
    nanoflann_wrapper tree_target(V_target.transpose());
    nanoflann_wrapper tree_source(V_source.transpose());
    std::vector<int> closest_point_indices(K);

    // for loop
    progressbar bar(V_source.cols()/skip_points);
    for (int i = 0; i < V_source.cols(); i = i+skip_points) {
        int point_index_original = i*skip_points;

        // find the correspondence of the source on the target
        int correspondence_on_target = find_correspondence("weighted", 
                                                           V_source.col(point_index_original),
                                                           N_source.col(point_index_original), 
	                                                       V_target,
	                                                       N_target,
                                                           tree_target,
                                                           K);

        // find the correspondence of the source on the target
        int correspondence_on_source = find_correspondence("weighted", 
                                                           V_target.col(correspondence_on_target),
                                                           N_target.col(correspondence_on_target), 
	                                                       V_source,
	                                                       N_source,
                                                           tree_source,
                                                           K);

        // check the distance between the points
        if ( (V_source.col(point_index_original)-V_source.col(correspondence_on_source)).norm() < distance_threshold ) {
            correspondences_on_source.push_back(V_source.col(point_index_original));
            correspondences_on_target.push_back(V_target.col(correspondence_on_target));
        }

        bar.update();
    }

    // push back into an Eigen Matrices
    source_position.resize(3, correspondences_on_source.size());
    for (int i=0; i<correspondences_on_source.size(); i++)
        source_position.col(i) = correspondences_on_source.at(i);
    
    target_position.resize(3, correspondences_on_target.size());
    for (int i=0; i<correspondences_on_target.size(); i++)
        target_position.col(i) = correspondences_on_target.at(i);
    
    /*
    int x=0, y=1, z=2;
    //for (int i = 0; i < template_cloud->size(); i = i+10)
    //{
    //    template_cloud_downsampled->points.push_back(template_cloud->points[i]);
    //    template_normals_downsampled->points.push_back(template_normals->points[i]);
    //}
    double distance_threshold = 0.05;                   // used to check that the correspondence of the correspondence is not too far
    int skip_points = 10;                               // pseudo downsampling
    int K = 50;                                         // limit the search space
    int number_of_correspondences = int(round(V_source.cols()/skip_points));

    // generate downsampled cloud
    Eigen::MatrixXd V_source_downsampled, N_source_downsampled;
    V_source_downsampled.resize(3, number_of_correspondences);
    N_source_downsampled.resize(3, number_of_correspondences);
    for (int i = 0; i < number_of_correspondences; i++) {
        V_source_downsampled.col(i) = V_source.col(i*skip_points);
        N_source_downsampled.col(i) = N_source.col(i*skip_points);
    }

    // create kd-tree
    nanoflann_wrapper tree(V_target.transpose());
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    pointIdxNKNSearch.clear();

    Eigen::Vector3d median_point;

    source_position.resize(3, number_of_correspondences);
    target_position.resize(3, number_of_correspondences);


    for (int i = 0; i < number_of_correspondences; i++)
    {
        // WTF this is not used ????? LOOOOL
		pointIdxNKNSearch = tree.return_k_closest_points(V_source.col(i*skip_points).transpose(), K);

        //median_point = point_association("closest_point",
        //                                    V_source.col(i*skip_points),
        //                                    N_source.col(i*skip_points),
        //                                    V_target,
        //                                    N_target);
        
        median_point = V_target.col(pointIdxNKNSearch[0]);
        
        target_position.col(i) <<  median_point(x), median_point(y), median_point(z);
        source_position.col(i) << V_source.col(i*skip_points)(x),
                                  V_source.col(i*skip_points)(y),
                                  V_source.col(i*skip_points)(z);
    }
    */
}

