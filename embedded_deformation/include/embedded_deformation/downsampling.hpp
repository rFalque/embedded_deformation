/*
*   Greedy Search
*   by R. Falque
*   29/11/2018
*/

#ifndef DOWNSAMPLING
#define DOWNSAMPLING

#include <Eigen/Core>
#include <vector>
#include <limits> 

#include <cfloat>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "getMinMax.hpp"
#include "farther_sampling.hpp"
#include "nanoflannWrapper.hpp"

using namespace std;

class three_d_point{
public:
	double x;
	double y;
	double z;
};

class point_and_occurences{
public:
	double x;
	double y;
	double z;
	int occurence;
};

inline void voxel_grid_downsampling(Eigen::MatrixXd & in_cloud, double leaf_size, Eigen::MatrixXd & out_cloud)
{
	Eigen::Vector3d min_point, max_point;
	getMinMax(in_cloud, min_point, max_point);

	double inv_leaf_size;
	inv_leaf_size = 1.0/leaf_size;

	Eigen::Vector3i min_box, max_box;
	min_box << floor(min_point(0) * inv_leaf_size ), floor(min_point(1) * inv_leaf_size ), floor(min_point(2) * inv_leaf_size); 
	max_box << floor(max_point(0) * inv_leaf_size ), floor(max_point(1) * inv_leaf_size ), floor(max_point(2) * inv_leaf_size); 

    Eigen::Vector3i divb, divb_mul;
    divb << max_box(0) - min_box(0) + 1, max_box(1) - min_box(1) + 1, max_box(2) - min_box(2) + 1;
    divb_mul << 1, divb(0), divb(0) * divb(1);

	std::vector < std::vector < std::vector < point_and_occurences> > > voxels;

	voxels.resize(divb(0));
	for (int x_index = 0; x_index < voxels.size(); ++x_index)
	{
		voxels[x_index].resize(divb(1));
		for (int y_index = 0; y_index < voxels[0].size(); ++y_index)
		{
			voxels[x_index][y_index].resize(divb(2));
		}
	}

	// plus assign zeros to voxel_count
	for (int i = 0; i < in_cloud.rows(); ++i)
	{
        int x_index = static_cast<int> ( floor(in_cloud(i, 0) * inv_leaf_size) - min_box(0) );
        int y_index = static_cast<int> ( floor(in_cloud(i, 1) * inv_leaf_size) - min_box(1) );
        int z_index = static_cast<int> ( floor(in_cloud(i, 2) * inv_leaf_size) - min_box(2) );

        voxels[x_index][y_index][z_index].x += in_cloud(i,0);
        voxels[x_index][y_index][z_index].y += in_cloud(i,1);
        voxels[x_index][y_index][z_index].z += in_cloud(i,2);
        voxels[x_index][y_index][z_index].occurence ++;
	}

	std::vector< three_d_point> final_cloud;
	three_d_point temp;
	for (int x_index = 0; x_index < voxels.size(); ++x_index)
	{
		for (int y_index = 0; y_index < voxels[0].size(); ++y_index)
		{
			for (int z_index = 0; z_index < voxels[0][0].size(); ++z_index)
			{
				if (voxels[x_index][y_index][z_index].occurence!= 0)
				{
					temp.x = voxels[x_index][y_index][z_index].x / voxels[x_index][y_index][z_index].occurence;
					temp.y = voxels[x_index][y_index][z_index].y / voxels[x_index][y_index][z_index].occurence;
					temp.z = voxels[x_index][y_index][z_index].z / voxels[x_index][y_index][z_index].occurence;
					final_cloud.push_back(temp);
				}
			}
		}
	}

	out_cloud.resize(final_cloud.size(), 3);
	for (int i = 0; i < final_cloud.size(); ++i)
	{
		out_cloud.row(i) << final_cloud[i].x, final_cloud[i].y, final_cloud[i].z;
	}

};


inline void downsampling(Eigen::MatrixXd & in_cloud, 
                         Eigen::MatrixXd & out_cloud, 
						 std::vector<int> & in_cloud_samples, 
						 double grid_resolution,
						 double leaf_size, 
						 bool use_farthest_sampling, 
						 bool use_relative_grid)
{
	// overwrite the leaf_size
	if (use_relative_grid) {
		double scale;
		getScale(in_cloud, scale);
		leaf_size = scale / grid_resolution;
	}

	// downsampling
	Eigen::MatrixXd downsampled_cloud;
	if (use_farthest_sampling)
	{
		farthest_sampling_by_sphere(in_cloud, leaf_size/100, downsampled_cloud);
	}
	else
	{
		voxel_grid_downsampling(in_cloud, leaf_size, downsampled_cloud);
	}

	out_cloud.resize(downsampled_cloud.rows(), 3);

	nanoflann_wrapper tree(in_cloud);
	for (int i = 0; i < downsampled_cloud.rows(); ++i)
	{
		std::vector< int > closest_point;
		closest_point = tree.return_k_closest_points(downsampled_cloud.row(i), 1);

		out_cloud.row(i) = in_cloud.row( closest_point[0] );
		in_cloud_samples.push_back(closest_point[0]);
	}

};

#endif
