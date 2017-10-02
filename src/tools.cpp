#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float MIN_VALUE = 0.001;
const float LEAST_VALUE = 0.1;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size

  if (estimations.size() == 0)
    return rmse;

  if (estimations.size() != ground_truth.size() )
    return rmse;

  VectorXd diff_vec;
  VectorXd rmse_vec;
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    //take the diff of the prediction vs. ground truth
    diff_vec = estimations[i] - ground_truth[i];
    //Take the sum of square values of the difference
    rmse_vec = diff_vec.array()*diff_vec.array();
    rmse = rmse + rmse_vec;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	Hj << 0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0;

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //Calculate the components required to calculate Jacobian matrix elements 
  float sq_term = std::pow(px, 2) + std::pow(py, 2);
  float sq_r_term = std::pow(sq_term, 0.5);
  float sq_three_term = std::pow(sq_term, 1.5);
  float cross_Term = vx*py - vy*px;
  float cross_Term1 = -vx*py + vy*px;
    
    
  //check division by zero
  if (sq_term == 0)
    return Hj;
  
  //Here, ele_<i><j> represents the vector element at row <i> and col <j>  
  float ele_00 = px / sq_r_term;
  float ele_01 = py / sq_r_term;
  float ele_10 = -1*py/sq_term;
  float ele_11 = px/sq_term;
  float ele_20 = (py*cross_Term)/sq_three_term;
  float ele_21 = (px*cross_Term1)/sq_three_term;
  float ele_22 = ele_00;
  float ele_23 = ele_01;

  //compute the Jacobian matrix
  Hj << ele_00, ele_01, 0, 0,
     ele_10, ele_11, 0, 0,
     ele_20, ele_21, ele_22, ele_23;

  return Hj;
}

// Helper funxtion to set min values
float Tools::SetMinValues(float in) {
	// if (pX < leastValue)
	if (fabs(in) < MIN_VALUE){
		std::cout << "LOG: min value " << in << std::endl;

		in = LEAST_VALUE;
	}
	return in;
}
