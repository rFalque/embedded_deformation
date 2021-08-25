#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>

// equation (3) in KinEtre: norm2(R'*R-I)^2
// gradient computed with http://www.matrixcalculus.org/
class RigCostFunction: public ceres::CostFunction
{
public:
	// constructor
	RigCostFunction(double cost_function_weight){
		// define residual and block_size
		set_num_residuals(1);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		block_sizes->push_back(12);
		block_sizes->push_back(12);

		cost_function_weight_ = sqrt(cost_function_weight);
	}

	~RigCostFunction(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const
    {
        Eigen::Map<const Eigen::Matrix3d> R_j(parameters[0]);
        Eigen::Map<const Eigen::Vector3d> t_j(parameters[0]+9);
        Eigen::Map<const Eigen::Matrix3d> R_k(parameters[1]);
        Eigen::Map<const Eigen::Vector3d> t_k(parameters[1]+9);
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3,3);

        Eigen::Map<Eigen::MatrixXd>res(residuals,1,3);

        double alpha_j_k = 1;

        //res(0, 0) = cost_function_weight_*( (R_j.transpose()*R_k - I).squaredNorm() + (t_k - t_j).squaredNorm() );
        res(0, 0) = cost_function_weight_*(R_j.transpose()*R_k - I).squaredNorm();
        
        if (jacobians!=NULL){
            if(jacobians[0]!=NULL){
                /*
                Eigen::Map<Eigen::MatrixXd>first_param_block_gradient(jacobians[0],1,12);
                Eigen::Matrix3d rotation_gradient;
                Eigen::Vector3d translation_gradient;

                rotation_gradient = 2 * R_k * (R_k.transpose()*R_j - I);
                translation_gradient = 2 * (t_j-t_k);
                Eigen::Map<Eigen::Matrix<double, 1, 9>, Eigen::AutoAlign> inline_matrix(rotation_gradient.data());

                first_param_block_gradient.block(0, 0, 1, 9) = inline_matrix;
                first_param_block_gradient(0,9) = translation_gradient(0);
                first_param_block_gradient(0,10) = translation_gradient(1);
                first_param_block_gradient(0,11) = translation_gradient(2);

                Eigen::Map<Eigen::MatrixXd>second_param_block_gradient(jacobians[1],1,12);

                rotation_gradient = 2 * R_j * (R_j.transpose()*R_k - I);
                translation_gradient = -2 * (t_j-t_k);
                Eigen::Map<Eigen::Matrix<double, 1, 9>, Eigen::AutoAlign> inline_matrix_2(rotation_gradient.data());

                second_param_block_gradient.block(0, 0, 1, 9) = inline_matrix_2;
                second_param_block_gradient(0,9) = translation_gradient(0);
                second_param_block_gradient(0,10) = translation_gradient(1);
                second_param_block_gradient(0,11) = translation_gradient(2);
                */

                Eigen::Map<Eigen::MatrixXd>first_param_block_gradient(jacobians[0],1,12);
                Eigen::Matrix3d rotation_gradient;
                Eigen::Vector3d translation_gradient;

                rotation_gradient = cost_function_weight_ * 2 * R_k * (R_k.transpose()*R_j - I);
                translation_gradient = 2 * (t_j-t_k);
                Eigen::Map<Eigen::Matrix<double, 1, 9>, Eigen::AutoAlign> inline_matrix(rotation_gradient.data());

                first_param_block_gradient.block(0, 0, 1, 9) = inline_matrix;
                first_param_block_gradient(0,9) = 0;
                first_param_block_gradient(0,10) = 0;
                first_param_block_gradient(0,11) = 0;

                Eigen::Map<Eigen::MatrixXd>second_param_block_gradient(jacobians[1],1,12);

                rotation_gradient = 2 * R_j * (R_j.transpose()*R_k - I);
                translation_gradient = -2 * (t_j-t_k);
                Eigen::Map<Eigen::Matrix<double, 1, 9>, Eigen::AutoAlign> inline_matrix_2(rotation_gradient.data());

                second_param_block_gradient.block(0, 0, 1, 9) = inline_matrix_2;
                second_param_block_gradient(0,9) = 0;
                second_param_block_gradient(0,10) = 0;
                second_param_block_gradient(0,11) = 0;
            }
        }
        return true;
	}

private:
	Eigen::Vector3d g_j_;
	Eigen::Vector3d g_k_;
	double cost_function_weight_;
};


// equation (5-6)
class RotCostFunction: public ceres::CostFunction
{
public:

	// constructor
	RotCostFunction(double cost_function_weight){
		// define residual and block_size
		set_num_residuals(6);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		block_sizes->push_back(12);

		cost_function_weight_ = sqrt(cost_function_weight);
	}

	~RotCostFunction(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const{

			Eigen::Map<const Eigen::Matrix3d> R_j(parameters[0]);
			Eigen::Vector3d c_1 = R_j.col(0);
			Eigen::Vector3d c_2 = R_j.col(1);
			Eigen::Vector3d c_3 = R_j.col(2);

			Eigen::Map<Eigen::MatrixXd>res(residuals,1,6);
			res(0, 0) = cost_function_weight_ * c_1.dot( c_2 );
			res(0, 1) = cost_function_weight_ * c_1.dot( c_3 );
			res(0, 2) = cost_function_weight_ * c_2.dot( c_3 );
			res(0, 3) = cost_function_weight_ * ( c_1.dot( c_1 )-1 );
			res(0, 4) = cost_function_weight_ * ( c_2.dot( c_2 )-1 );
			res(0, 5) = cost_function_weight_ * ( c_3.dot( c_3 )-1 );

			if (jacobians!=NULL){
				if(jacobians[0]!=NULL){
					// derivated block for (c_1 . c_2)^2
					jacobians[0][ 0] = cost_function_weight_ * c_2(0);		// d( (c_1.c_2)^2 ) / d( R_j(0,0) )
					jacobians[0][ 1] = cost_function_weight_ * c_2(1);		// d( (c_1.c_2)^2 ) / d( R_j(1,0) )
					jacobians[0][ 2] = cost_function_weight_ * c_2(2);		// d( (c_1.c_2)^2 ) / d( R_j(2,0) )
					jacobians[0][ 3] = cost_function_weight_ * c_1(0);		// d( (c_1.c_2)^2 ) / d( R_j(0,1) )
					jacobians[0][ 4] = cost_function_weight_ * c_1(1);		// d( (c_1.c_2)^2 ) / d( R_j(1,1) )
					jacobians[0][ 5] = cost_function_weight_ * c_1(2);;		// d( (c_1.c_2)^2 ) / d( R_j(2,1) )
					jacobians[0][ 6] = 0;									// d( (c_1.c_2)^2 ) / d( R_j(0,2) )
					jacobians[0][ 7] = 0;									// d( (c_1.c_2)^2 ) / d( R_j(1,2) )
					jacobians[0][ 8] = 0;									// d( (c_1.c_2)^2 ) / d( R_j(2,2) )
					jacobians[0][ 9] = 0;									// d( (c_1.c_2)^2 ) / d( T_j(0) )
					jacobians[0][10] = 0;									// d( (c_1.c_2)^2 ) / d( T_j(1) )
					jacobians[0][11] = 0;									// d( (c_1.c_2)^2 ) / d( T_j(2) )

					// derivated block for (c_1 . c_3)^2
					jacobians[0][ 0 + 12] = cost_function_weight_ * c_3(0);
					jacobians[0][ 1 + 12] = cost_function_weight_ * c_3(1);
					jacobians[0][ 2 + 12] = cost_function_weight_ * c_3(2);
					jacobians[0][ 3 + 12] = 0;
					jacobians[0][ 4 + 12] = 0;
					jacobians[0][ 5 + 12] = 0;
					jacobians[0][ 6 + 12] = cost_function_weight_ * c_1(0);
					jacobians[0][ 7 + 12] = cost_function_weight_ * c_1(1);
					jacobians[0][ 8 + 12] = cost_function_weight_ * c_1(2);
					jacobians[0][ 9 + 12] = 0;
					jacobians[0][10 + 12] = 0;
					jacobians[0][11 + 12] = 0;

					// derivated block for (c_2 . c_3)^2
					jacobians[0][ 0 + 12*2] = 0;
					jacobians[0][ 1 + 12*2] = 0;
					jacobians[0][ 2 + 12*2] = 0;
					jacobians[0][ 3 + 12*2] = cost_function_weight_ * c_3(0);
					jacobians[0][ 4 + 12*2] = cost_function_weight_ * c_3(1);
					jacobians[0][ 5 + 12*2] = cost_function_weight_ * c_3(2);
					jacobians[0][ 6 + 12*2] = cost_function_weight_ * c_2(0);
					jacobians[0][ 7 + 12*2] = cost_function_weight_ * c_2(1);
					jacobians[0][ 8 + 12*2] = cost_function_weight_ * c_2(2);
					jacobians[0][ 9 + 12*2] = 0;
					jacobians[0][10 + 12*2] = 0;
					jacobians[0][11 + 12*2] = 0;

					// derivated block for (C_1.C_1 - 1)^2
					jacobians[0][ 0 + 12*3] = cost_function_weight_ * 2*c_1(0);
					jacobians[0][ 1 + 12*3] = cost_function_weight_ * 2*c_1(1);
					jacobians[0][ 2 + 12*3] = cost_function_weight_ * 2*c_1(2);
					jacobians[0][ 3 + 12*3] = 0;
					jacobians[0][ 4 + 12*3] = 0;
					jacobians[0][ 5 + 12*3] = 0;
					jacobians[0][ 6 + 12*3] = 0;
					jacobians[0][ 7 + 12*3] = 0;
					jacobians[0][ 8 + 12*3] = 0;
					jacobians[0][ 9 + 12*3] = 0;
					jacobians[0][10 + 12*3] = 0;
					jacobians[0][11 + 12*3] = 0;

					// derivated block for (C_2.C_2 - 1)^2
					jacobians[0][0 + 12*4] = 0;
					jacobians[0][1 + 12*4] = 0;
					jacobians[0][2 + 12*4] = 0;
					jacobians[0][3 + 12*4] = cost_function_weight_ * 2*c_2(0);
					jacobians[0][4 + 12*4] = cost_function_weight_ * 2*c_2(1);
					jacobians[0][5 + 12*4] = cost_function_weight_ * 2*c_2(2);
					jacobians[0][6 + 12*4] = 0;
					jacobians[0][7 + 12*4] = 0;
					jacobians[0][8 + 12*4] = 0;
					jacobians[0][9 + 12*4] = 0;
					jacobians[0][10 + 12*4] = 0;
					jacobians[0][11 + 12*4] = 0;

					// derivated block for (C_3.C_3 - 1)^2
					jacobians[0][0 + 12*5] = 0;
					jacobians[0][1 + 12*5] = 0;
					jacobians[0][2 + 12*5] = 0;
					jacobians[0][3 + 12*5] = 0;
					jacobians[0][4 + 12*5] = 0;
					jacobians[0][5 + 12*5] = 0;
					jacobians[0][6 + 12*5] = cost_function_weight_ * 2*c_3(0);
					jacobians[0][7 + 12*5] = cost_function_weight_ * 2*c_3(1);
					jacobians[0][8 + 12*5] = cost_function_weight_ * 2*c_3(2);
					jacobians[0][9 + 12*5] = 0;
					jacobians[0][10 + 12*5] = 0;
					jacobians[0][11 + 12*5] = 0;
				}
			}
			return true;
	}

private:
	double cost_function_weight_;
};


// equation (3) in KinEtre: norm2(R'*R-I)^2
class RotCostFunction_2: public ceres::CostFunction
{
public:

	// constructor
	RotCostFunction_2(double cost_function_weight){
		// define residual and block_size
		set_num_residuals(9);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		block_sizes->push_back(12);

		cost_function_weight_ = sqrt(cost_function_weight);
	}

	~RotCostFunction_2(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const
    {
        Eigen::Map<const Eigen::Matrix3d> R(parameters[0]);
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3,3);

        Eigen::Map<Eigen::MatrixXd>res(residuals,3,3);
        res = cost_function_weight_*(R.transpose()*R - I);

        if (jacobians!=NULL){
            if(jacobians[0]!=NULL){

				//R'*R - I (0,0)
				jacobians[0][ 0] = cost_function_weight_*2*R(0,0);	// df/dR(0,0)
				jacobians[0][ 1] = cost_function_weight_*2*R(1,0);	// df/dR(1,0)
				jacobians[0][ 2] = cost_function_weight_*2*R(2,0);	// df/dR(2,0)
				jacobians[0][ 3] = 0;								// df/dR(0,1)
				jacobians[0][ 4] = 0;								// df/dR(1,1)
				jacobians[0][ 5] = 0;								// df/dR(2,1)
				jacobians[0][ 6] = 0;								// df/dR(0,2)
				jacobians[0][ 7] = 0;								// df/dR(1,2)
				jacobians[0][ 8] = 0;								// df/dR(2,2)
				jacobians[0][ 9] = 0;								// df/dT(0)
				jacobians[0][10] = 0;								// df/dT(1)
				jacobians[0][11] = 0;								// df/dT(2)

				// R'*R - I (0,1)
				jacobians[0][ 0 + 12] = cost_function_weight_*R(0,1);
				jacobians[0][ 1 + 12] = cost_function_weight_*R(1,1);
				jacobians[0][ 2 + 12] = cost_function_weight_*R(2,1);
				jacobians[0][ 3 + 12] = cost_function_weight_*R(0,0);
				jacobians[0][ 4 + 12] = cost_function_weight_*R(1,0);
				jacobians[0][ 5 + 12] = cost_function_weight_*R(2,0);
				jacobians[0][ 6 + 12] = 0;
				jacobians[0][ 7 + 12] = 0;
				jacobians[0][ 8 + 12] = 0;
				jacobians[0][ 9 + 12] = 0;
				jacobians[0][10 + 12] = 0;
				jacobians[0][11 + 12] = 0;

				// R'*R - I (0,2)
				jacobians[0][ 0 + 2*12] = cost_function_weight_*R(0,2);
				jacobians[0][ 1 + 2*12] = cost_function_weight_*R(1,2);
				jacobians[0][ 2 + 2*12] = cost_function_weight_*R(2,2);
				jacobians[0][ 3 + 2*12] = 0;
				jacobians[0][ 4 + 2*12] = 0;
				jacobians[0][ 5 + 2*12] = 0;
				jacobians[0][ 6 + 2*12] = cost_function_weight_*R(0,0);
				jacobians[0][ 7 + 2*12] = cost_function_weight_*R(1,0);
				jacobians[0][ 8 + 2*12] = cost_function_weight_*R(2,0);
				jacobians[0][ 9 + 2*12] = 0;
				jacobians[0][10 + 2*12] = 0;
				jacobians[0][11 + 2*12] = 0;

				// R'*R - I (1,0)
				jacobians[0][ 0 + 3*12] = cost_function_weight_*R(0,1);
				jacobians[0][ 1 + 3*12] = cost_function_weight_*R(1,1);
				jacobians[0][ 2 + 3*12] = cost_function_weight_*R(2,1);
				jacobians[0][ 3 + 3*12] = cost_function_weight_*R(0,0);
				jacobians[0][ 4 + 3*12] = cost_function_weight_*R(1,0);
				jacobians[0][ 5 + 3*12] = cost_function_weight_*R(2,0);
				jacobians[0][ 6 + 3*12] = 0;
				jacobians[0][ 7 + 3*12] = 0;
				jacobians[0][ 8 + 3*12] = 0;
				jacobians[0][ 9 + 3*12] = 0;
				jacobians[0][10 + 3*12] = 0;
				jacobians[0][11 + 3*12] = 0;

				// R'*R - I (1,1)
				jacobians[0][ 0 + 4*12] = 0;
				jacobians[0][ 1 + 4*12] = 0;
				jacobians[0][ 2 + 4*12] = 0;
				jacobians[0][ 3 + 4*12] = cost_function_weight_*2*R(0,1);
				jacobians[0][ 4 + 4*12] = cost_function_weight_*2*R(1,1);
				jacobians[0][ 5 + 4*12] = cost_function_weight_*2*R(2,1);
				jacobians[0][ 6 + 4*12] = 0;
				jacobians[0][ 7 + 4*12] = 0;
				jacobians[0][ 8 + 4*12] = 0;
				jacobians[0][ 9 + 4*12] = 0;
				jacobians[0][10 + 4*12] = 0;
				jacobians[0][11 + 4*12] = 0;

				// R'*R - I (1,2)
				jacobians[0][ 0 + 5*12] = 0;
				jacobians[0][ 1 + 5*12] = 0;
				jacobians[0][ 2 + 5*12] = 0;
				jacobians[0][ 3 + 5*12] = cost_function_weight_*R(0,2);
				jacobians[0][ 4 + 5*12] = cost_function_weight_*R(1,2);
				jacobians[0][ 5 + 5*12] = cost_function_weight_*R(2,2);
				jacobians[0][ 6 + 5*12] = cost_function_weight_*R(0,1);
				jacobians[0][ 7 + 5*12] = cost_function_weight_*R(1,1);
				jacobians[0][ 8 + 5*12] = cost_function_weight_*R(2,1);
				jacobians[0][ 9 + 5*12] = 0;
				jacobians[0][10 + 5*12] = 0;
				jacobians[0][11 + 5*12] = 0;

				// R'*R - I (2,0)
				jacobians[0][ 0 + 6*12] = cost_function_weight_*R(0,2);
				jacobians[0][ 1 + 6*12] = cost_function_weight_*R(1,2);
				jacobians[0][ 2 + 6*12] = cost_function_weight_*R(2,2);
				jacobians[0][ 3 + 6*12] = 0;
				jacobians[0][ 4 + 6*12] = 0;
				jacobians[0][ 5 + 6*12] = 0;
				jacobians[0][ 6 + 6*12] = cost_function_weight_*R(0,0);
				jacobians[0][ 7 + 6*12] = cost_function_weight_*R(1,0);
				jacobians[0][ 8 + 6*12] = cost_function_weight_*R(2,0);
				jacobians[0][ 9 + 6*12] = 0;
				jacobians[0][10 + 6*12] = 0;
				jacobians[0][11 + 6*12] = 0;

				// R'*R - I (2,1)
				jacobians[0][ 0 + 7*12] = 0;
				jacobians[0][ 1 + 7*12] = 0;
				jacobians[0][ 2 + 7*12] = 0;
				jacobians[0][ 3 + 7*12] = cost_function_weight_*R(0,2);
				jacobians[0][ 4 + 7*12] = cost_function_weight_*R(1,2);
				jacobians[0][ 5 + 7*12] = cost_function_weight_*R(2,2);
				jacobians[0][ 6 + 7*12] = cost_function_weight_*R(0,1);
				jacobians[0][ 7 + 7*12] = cost_function_weight_*R(1,1);
				jacobians[0][ 8 + 7*12] = cost_function_weight_*R(2,1);
				jacobians[0][ 9 + 7*12] = 0;
				jacobians[0][10 + 7*12] = 0;
				jacobians[0][11 + 7*12] = 0;

				// R'*R - I (2,2)
				jacobians[0][ 0 + 8*12] = 0;
				jacobians[0][ 1 + 8*12] = 0;
				jacobians[0][ 2 + 8*12] = 0;
				jacobians[0][ 3 + 8*12] = 0;
				jacobians[0][ 4 + 8*12] = 0;
				jacobians[0][ 5 + 8*12] = 0;
				jacobians[0][ 6 + 8*12] = cost_function_weight_*2*R(0,2);
				jacobians[0][ 7 + 8*12] = cost_function_weight_*2*R(1,2);
				jacobians[0][ 8 + 8*12] = cost_function_weight_*2*R(2,2);
				jacobians[0][ 9 + 8*12] = 0;
				jacobians[0][10 + 8*12] = 0;
				jacobians[0][11 + 8*12] = 0;

            }
        }
        return true;
	}

private:
	double cost_function_weight_;
};


// equation (7)
class RegCostFunction: public ceres::CostFunction
{
public:
	// constructor
	RegCostFunction(double cost_function_weight, Eigen::Vector3d g_j, Eigen::Vector3d g_k){
		// define residual and block_size
		set_num_residuals(3);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		block_sizes->push_back(12);
		block_sizes->push_back(12);

		g_j_ = g_j;
		g_k_ = g_k;
		cost_function_weight_ = sqrt(cost_function_weight);
	}

	~RegCostFunction(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const{

			Eigen::Map<const Eigen::Matrix3d> R_j(parameters[0]);
			Eigen::Map<const Eigen::Vector3d> t_j(parameters[0]+9);
			Eigen::Map<const Eigen::Matrix3d> R_k(parameters[1]);
			Eigen::Map<const Eigen::Vector3d> t_k(parameters[1]+9);

			Eigen::Map<Eigen::Vector3d>res(residuals,3);

			double alpha_j_k = 1;

			res = cost_function_weight_ * alpha_j_k * ( R_j*(g_k_ - g_j_) + g_j_ + t_j - (g_k_ + t_k) );

			if (jacobians!=NULL){
				if(jacobians[0]!=NULL){
					jacobians[0][ 0] = cost_function_weight_ * alpha_j_k * ( g_k_(0) - g_j_(0) );
					jacobians[0][ 1] = 0;
					jacobians[0][ 2] = 0;
					jacobians[0][ 3] = cost_function_weight_ * alpha_j_k * ( g_k_(1) - g_j_(1) );
					jacobians[0][ 4] = 0;
					jacobians[0][ 5] = 0;
					jacobians[0][ 6] = cost_function_weight_ * alpha_j_k * ( g_k_(2) - g_j_(2) );
					jacobians[0][ 7] = 0;
					jacobians[0][ 8] = 0;
					jacobians[0][ 9] = cost_function_weight_ * alpha_j_k;
					jacobians[0][10] = 0;
					jacobians[0][11] = 0;

					jacobians[0][ 0 + 12] = 0;
					jacobians[0][ 1 + 12] = cost_function_weight_ * alpha_j_k * ( g_k_(0) - g_j_(0) );
					jacobians[0][ 2 + 12] = 0;
					jacobians[0][ 3 + 12] = 0;
					jacobians[0][ 4 + 12] = cost_function_weight_ * alpha_j_k * ( g_k_(1) - g_j_(1) );
					jacobians[0][ 5 + 12] = 0;
					jacobians[0][ 6 + 12] = 0;
					jacobians[0][ 7 + 12] = cost_function_weight_ * alpha_j_k * ( g_k_(2) - g_j_(2) );
					jacobians[0][ 8 + 12] = 0;
					jacobians[0][ 9 + 12] = 0;
					jacobians[0][10 + 12] = cost_function_weight_ * alpha_j_k;
					jacobians[0][11 + 12] = 0;

					jacobians[0][ 0 + 12*2] = 0;
					jacobians[0][ 1 + 12*2] = 0;
					jacobians[0][ 2 + 12*2] = cost_function_weight_ * alpha_j_k * ( g_k_(0) - g_j_(0) );
					jacobians[0][ 3 + 12*2] = 0;
					jacobians[0][ 4 + 12*2] = 0;
					jacobians[0][ 5 + 12*2] = cost_function_weight_ * alpha_j_k * ( g_k_(1) - g_j_(1) );
					jacobians[0][ 6 + 12*2] = 0;
					jacobians[0][ 7 + 12*2] = 0;
					jacobians[0][ 8 + 12*2] = cost_function_weight_ * alpha_j_k * ( g_k_(2) - g_j_(2) );
					jacobians[0][ 9 + 12*2] = 0;
					jacobians[0][10 + 12*2] = 0;
					jacobians[0][11 + 12*2] = cost_function_weight_ * alpha_j_k;


					jacobians[1][ 0] = 0;
					jacobians[1][ 1] = 0;
					jacobians[1][ 2] = 0;
					jacobians[1][ 3] = 0;
					jacobians[1][ 4] = 0;
					jacobians[1][ 5] = 0;
					jacobians[1][ 6] = 0;
					jacobians[1][ 7] = 0;
					jacobians[1][ 8] = 0;
					jacobians[1][ 9] = - cost_function_weight_ * alpha_j_k;
					jacobians[1][10] = 0;
					jacobians[1][11] = 0;

					jacobians[1][ 0 + 12] = 0;
					jacobians[1][ 1 + 12] = 0;
					jacobians[1][ 2 + 12] = 0;
					jacobians[1][ 3 + 12] = 0;
					jacobians[1][ 4 + 12] = 0;
					jacobians[1][ 5 + 12] = 0;
					jacobians[1][ 6 + 12] = 0;
					jacobians[1][ 7 + 12] = 0;
					jacobians[1][ 8 + 12] = 0;
					jacobians[1][ 9 + 12] = 0;
					jacobians[1][10 + 12] = - cost_function_weight_ * alpha_j_k;
					jacobians[1][11 + 12] = 0;

					jacobians[1][ 0 + 12*2] = 0;
					jacobians[1][ 1 + 12*2] = 0;
					jacobians[1][ 2 + 12*2] = 0;
					jacobians[1][ 3 + 12*2] = 0;
					jacobians[1][ 4 + 12*2] = 0;
					jacobians[1][ 5 + 12*2] = 0;
					jacobians[1][ 6 + 12*2] = 0;
					jacobians[1][ 7 + 12*2] = 0;
					jacobians[1][ 8 + 12*2] = 0;
					jacobians[1][ 9 + 12*2] = 0;
					jacobians[1][10 + 12*2] = 0;
					jacobians[1][11 + 12*2] = - cost_function_weight_ * alpha_j_k;
				}
			}
			return true;
	}

private:
	Eigen::Vector3d g_j_;
	Eigen::Vector3d g_k_;
	double cost_function_weight_;
};


// equation (8)
class ConCostFunction: public ceres::CostFunction
{
public:
	// constructor
	ConCostFunction(double cost_function_weight, 
	                std::vector<Eigen::Vector3d> vector_g, 
					Eigen::Vector3d source, 
					Eigen::Vector3d target){
		// define residual and block_size
		source_ = source;
		target_ = target;
		vector_g_ = vector_g;
		cost_function_weight_ = sqrt(cost_function_weight);

		nodes_connectivity = vector_g_.size()-1;
		
		// set the size of the parameters input
		set_num_residuals(3);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		for (int i = 0; i < nodes_connectivity; ++i)
			block_sizes->push_back(12);
		
		// equation (3)
		w_j_.clear();
		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_.push_back( pow(1 - (source_ - vector_g_[i]).squaredNorm() / (source_ - vector_g_.back()).squaredNorm(), 2) );

		double normalization_factor = 0;
		for (int i = 0; i < nodes_connectivity; ++i)
			normalization_factor += w_j_[i];

		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_[i] /= normalization_factor;
	}

	~ConCostFunction(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const{
			Eigen::Map<Eigen::MatrixXd>res(residuals,1,3);

			// back to equation (8)
			Eigen::Vector3d new_node_position;
			new_node_position << 0, 0, 0;
			for (int i = 0; i < nodes_connectivity; ++i)
			{
				Eigen::Map<const Eigen::Matrix3d> R_j(parameters[i]);
				Eigen::Map<const Eigen::Vector3d> t_j(parameters[i]+9);
				new_node_position += w_j_[i] * (R_j * (source_ - vector_g_[i]) + vector_g_[i] + t_j);
			}

			res(0,0) = cost_function_weight_ * (new_node_position - target_)(0);
			res(0,1) = cost_function_weight_ * (new_node_position - target_)(1);
			res(0,2) = cost_function_weight_ * (new_node_position - target_)(2);


			if (jacobians!=NULL){
				if(jacobians[0]!=NULL){
					for (int i = 0; i < w_j_.size(); ++i)
					{
						jacobians[i][ 0] = cost_function_weight_ * w_j_[i]*( source_(0) - vector_g_[i](0) );
						jacobians[i][ 1] = 0;
						jacobians[i][ 2] = 0;
						jacobians[i][ 3] = cost_function_weight_ * w_j_[i]*( source_(1) - vector_g_[i](1) );
						jacobians[i][ 4] = 0;
						jacobians[i][ 5] = 0;
						jacobians[i][ 6] = cost_function_weight_ * w_j_[i]*( source_(2) - vector_g_[i](2) );
						jacobians[i][ 7] = 0;
						jacobians[i][ 8] = 0;
						jacobians[i][ 9] = cost_function_weight_ * w_j_[i];
						jacobians[i][10] = 0;
						jacobians[i][11] = 0;

						jacobians[i][ 0 + 12] = 0;
						jacobians[i][ 1 + 12] = cost_function_weight_ * w_j_[i]*( source_(0) - vector_g_[i](0) );
						jacobians[i][ 2 + 12] = 0;
						jacobians[i][ 3 + 12] = 0;
						jacobians[i][ 4 + 12] = cost_function_weight_ * w_j_[i]*( source_(1) - vector_g_[i](1) );
						jacobians[i][ 5 + 12] = 0;
						jacobians[i][ 6 + 12] = 0;
						jacobians[i][ 7 + 12] = cost_function_weight_ * w_j_[i]*( source_(2) - vector_g_[i](2) );
						jacobians[i][ 8 + 12] = 0;
						jacobians[i][ 9 + 12] = 0;
						jacobians[i][10 + 12] = cost_function_weight_ * w_j_[i];
						jacobians[i][11 + 12] = 0;

						jacobians[i][ 0 + 12*2] = 0;
						jacobians[i][ 1 + 12*2] = 0;
						jacobians[i][ 2 + 12*2] = cost_function_weight_ * w_j_[i]*( source_(0) - vector_g_[i](0) );
						jacobians[i][ 3 + 12*2] = 0;
						jacobians[i][ 4 + 12*2] = 0;
						jacobians[i][ 5 + 12*2] = cost_function_weight_ * w_j_[i]*( source_(1) - vector_g_[i](1) );
						jacobians[i][ 6 + 12*2] = 0;
						jacobians[i][ 7 + 12*2] = 0;
						jacobians[i][ 8 + 12*2] = cost_function_weight_ * w_j_[i]*( source_(2) - vector_g_[i](2) );
						jacobians[i][ 9 + 12*2] = 0;
						jacobians[i][10 + 12*2] = 0;
						jacobians[i][11 + 12*2] = cost_function_weight_ * w_j_[i];
					}
				}
			}
			return true;
	}

private:
	Eigen::Vector3d source_;
	Eigen::Vector3d target_;
	std::vector<Eigen::Vector3d> vector_g_;
	std::vector<double> w_j_;
	int nodes_connectivity;
	double cost_function_weight_;
};


// equation (20) in ElasticFusion (journal version) 
class RelCostFunction: public ceres::CostFunction
{
public:
	// constructor
	RelCostFunction(double cost_function_weight, 
					std::vector<Eigen::Vector3d> nodes_in_graph_1,
					std::vector<Eigen::Vector3d> nodes_in_graph_2,
					Eigen::Vector3d point_in_1, 
					Eigen::Vector3d point_in_2){
		// define residual and block_size
		point_in_1_ = point_in_1;
		point_in_2_ = point_in_2;
		nodes_in_graph_1_ = nodes_in_graph_1;
		nodes_in_graph_2_ = nodes_in_graph_2;
		cost_function_weight_ = sqrt(cost_function_weight);

		nodes_connectivity_1 = nodes_in_graph_1_.size()-1;
		nodes_connectivity_2 = nodes_in_graph_2_.size()-1;
		
		// set the size of the parameters input
		set_num_residuals(3);
		std::vector<int>* block_sizes =  mutable_parameter_block_sizes();
		for (int i = 0; i < nodes_connectivity_1; ++i)
			block_sizes->push_back(12);
		
		for (int i = 0; i < nodes_connectivity_2; ++i)
			block_sizes->push_back(12);

		// equation (3) tbd
		w_j_1_.clear();
		for (int i = 0; i < nodes_connectivity_1; ++i)
			w_j_1_.push_back( pow(1 - (point_in_1_ - nodes_in_graph_1_[i]).squaredNorm() / (point_in_1_ - nodes_in_graph_1_.back()).squaredNorm(), 2) ); // tbd

		double normalization_factor = 0;
		for (int i = 0; i < nodes_connectivity_1; ++i)
			normalization_factor += w_j_1_[i];

		for (int i = 0; i < nodes_connectivity_1; ++i)
			w_j_1_[i] /= normalization_factor;


		// equation (3) weight for right part of the equation
		w_j_2_.clear();
		for (int i = 0; i < nodes_connectivity_2; ++i)
			w_j_2_.push_back( pow(1 - (point_in_2_ - nodes_in_graph_2_[i]).squaredNorm() / (point_in_2_ - nodes_in_graph_2_.back()).squaredNorm(), 2) ); // tbd

		normalization_factor = 0;
		for (int i = 0; i < nodes_connectivity_2; ++i)
			normalization_factor += w_j_2_[i];

		for (int i = 0; i < nodes_connectivity_2; ++i)
			w_j_2_[i] /= normalization_factor;
	}

	~RelCostFunction(){
	}

	virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const{
			Eigen::Map<Eigen::MatrixXd>res(residuals,1,3);

			// back to equation (8)
			Eigen::Vector3d new_point_position_wrt_graph_1;
			Eigen::Vector3d new_point_position_wrt_graph_2;
			new_point_position_wrt_graph_1 << 0, 0, 0;
			for (int i = 0; i < nodes_connectivity_1; ++i)
			{
				Eigen::Map<const Eigen::Matrix3d> R_j(parameters[i]);
				Eigen::Map<const Eigen::Vector3d> t_j(parameters[i]+9);
				new_point_position_wrt_graph_1 += w_j_1_[i] * (R_j * (point_in_1_ - nodes_in_graph_1_[i])
				                                            + nodes_in_graph_1_[i] + t_j);
			}

			new_point_position_wrt_graph_2 << 0, 0, 0;
			for (int i = 0; i < nodes_connectivity_2; ++i)
			{
				Eigen::Map<const Eigen::Matrix3d> R_j(parameters[i+nodes_connectivity_1]);
				Eigen::Map<const Eigen::Vector3d> t_j(parameters[i+nodes_connectivity_1]+9);
				new_point_position_wrt_graph_2 += w_j_2_[i] * (R_j * (point_in_2_ - nodes_in_graph_2_[i]) 
				                                            + nodes_in_graph_2_[i] + t_j);
			}

			res(0,0) = cost_function_weight_ * (new_point_position_wrt_graph_1 - new_point_position_wrt_graph_2)(0);
			res(0,1) = cost_function_weight_ * (new_point_position_wrt_graph_1 - new_point_position_wrt_graph_2)(1);
			res(0,2) = cost_function_weight_ * (new_point_position_wrt_graph_1 - new_point_position_wrt_graph_2)(2);

			if (jacobians!=NULL){
				if(jacobians[0]!=NULL){
					for (int i = 0; i < nodes_connectivity_1; ++i)
					{
						jacobians[i][ 0] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(0) - nodes_in_graph_1_[i](0) );
						jacobians[i][ 1] = 0;
						jacobians[i][ 2] = 0;
						jacobians[i][ 3] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(1) - nodes_in_graph_1_[i](1) );
						jacobians[i][ 4] = 0;
						jacobians[i][ 5] = 0;
						jacobians[i][ 6] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(2) - nodes_in_graph_1_[i](2) );
						jacobians[i][ 7] = 0;
						jacobians[i][ 8] = 0;
						jacobians[i][ 9] = cost_function_weight_ * w_j_1_[i];
						jacobians[i][10] = 0;
						jacobians[i][11] = 0;

						jacobians[i][ 0 + 12] = 0;
						jacobians[i][ 1 + 12] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(0) - nodes_in_graph_1_[i](0) );
						jacobians[i][ 2 + 12] = 0;
						jacobians[i][ 3 + 12] = 0;
						jacobians[i][ 4 + 12] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(1) - nodes_in_graph_1_[i](1) );
						jacobians[i][ 5 + 12] = 0;
						jacobians[i][ 6 + 12] = 0;
						jacobians[i][ 7 + 12] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(2) - nodes_in_graph_1_[i](2) );
						jacobians[i][ 8 + 12] = 0;
						jacobians[i][ 9 + 12] = 0;
						jacobians[i][10 + 12] = cost_function_weight_ * w_j_1_[i];
						jacobians[i][11 + 12] = 0;

						jacobians[i][ 0 + 12*2] = 0;
						jacobians[i][ 1 + 12*2] = 0;
						jacobians[i][ 2 + 12*2] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(0) - nodes_in_graph_1_[i](0) );
						jacobians[i][ 3 + 12*2] = 0;
						jacobians[i][ 4 + 12*2] = 0;
						jacobians[i][ 5 + 12*2] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(1) - nodes_in_graph_1_[i](1) );
						jacobians[i][ 6 + 12*2] = 0;
						jacobians[i][ 7 + 12*2] = 0;
						jacobians[i][ 8 + 12*2] = cost_function_weight_ * w_j_1_[i]*( point_in_1_(2) - nodes_in_graph_1_[i](2) );
						jacobians[i][ 9 + 12*2] = 0;
						jacobians[i][10 + 12*2] = 0;
						jacobians[i][11 + 12*2] = cost_function_weight_ * w_j_1_[i];
					}
				}
			}


			if (jacobians!=NULL){
				if(jacobians[0]!=NULL){
					for (int i = 0; i < nodes_connectivity_2; ++i)
					{
						jacobians[i+nodes_connectivity_1][ 0] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(0) - nodes_in_graph_2_[i](0) );
						jacobians[i+nodes_connectivity_1][ 1] = 0;
						jacobians[i+nodes_connectivity_1][ 2] = 0;
						jacobians[i+nodes_connectivity_1][ 3] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(1) - nodes_in_graph_2_[i](1) );
						jacobians[i+nodes_connectivity_1][ 4] = 0;
						jacobians[i+nodes_connectivity_1][ 5] = 0;
						jacobians[i+nodes_connectivity_1][ 6] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(2) - nodes_in_graph_2_[i](2) );
						jacobians[i+nodes_connectivity_1][ 7] = 0;
						jacobians[i+nodes_connectivity_1][ 8] = 0;
						jacobians[i+nodes_connectivity_1][ 9] = -cost_function_weight_ * w_j_2_[i];
						jacobians[i+nodes_connectivity_1][10] = 0;
						jacobians[i+nodes_connectivity_1][11] = 0;

						jacobians[i+nodes_connectivity_1][ 0 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 1 + 12] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(0) - nodes_in_graph_2_[i](0) );
						jacobians[i+nodes_connectivity_1][ 2 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 3 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 4 + 12] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(1) - nodes_in_graph_2_[i](1) );
						jacobians[i+nodes_connectivity_1][ 5 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 6 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 7 + 12] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(2) - nodes_in_graph_2_[i](2) );
						jacobians[i+nodes_connectivity_1][ 8 + 12] = 0;
						jacobians[i+nodes_connectivity_1][ 9 + 12] = 0;
						jacobians[i+nodes_connectivity_1][10 + 12] = -cost_function_weight_ * w_j_2_[i];
						jacobians[i+nodes_connectivity_1][11 + 12] = 0;

						jacobians[i+nodes_connectivity_1][ 0 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 1 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 2 + 12*2] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(0) - nodes_in_graph_2_[i](0) );
						jacobians[i+nodes_connectivity_1][ 3 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 4 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 5 + 12*2] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(1) - nodes_in_graph_2_[i](1) );
						jacobians[i+nodes_connectivity_1][ 6 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 7 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][ 8 + 12*2] = -cost_function_weight_ * w_j_2_[i]*( point_in_2_(2) - nodes_in_graph_2_[i](2) );
						jacobians[i+nodes_connectivity_1][ 9 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][10 + 12*2] = 0;
						jacobians[i+nodes_connectivity_1][11 + 12*2] = -cost_function_weight_ * w_j_2_[i];
					}
				}
			}
			return true;
	}

private:
	Eigen::Vector3d point_in_1_;
	Eigen::Vector3d point_in_2_;
	std::vector<Eigen::Vector3d> nodes_in_graph_1_;
	std::vector<Eigen::Vector3d> nodes_in_graph_2_;
	std::vector<double> w_j_1_;
	std::vector<double> w_j_2_;
	int nodes_connectivity_1;
	int nodes_connectivity_2;
	double cost_function_weight_;
};
