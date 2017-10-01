#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

const float PI = 3.14159265;
  
KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //Find h(x) for Radar msmts 
	//Get the predicted values of (px, py, vx, vy) from _x
	float pX = x_(0);
	float pY = x_(1);
	float vX = x_(2);
	float vY = x_(3);
	float sq_term = std::pow(pX, 2) + std::pow(pY, 2);
  float sq_r_term = std::pow(sq_term, 0.5);
	float rho = sq_r_term;
	float theta = atan2(pY, pX);
	float rho_dot = (pX*vX+pY*vY)/sq_r_term;

	VectorXd z_pred = VectorXd(3);
  z_pred << rho, theta, rho_dot;

	VectorXd y = z - z_pred;
  std::cout << "LOG: theta in y = " << y(1) << std::endl;
	//check for the range of y(1) - theta to be within -PI to PI
	while((y(1) > PI) || (y(1) < -1*PI))
	{
		if (y(1) > PI){
			y(1) -= 2*PI;
		} else if (y(1) < -1*PI){
			y(1) += 2*PI;
		}
	}
	std::cout << "LOG: theta after in y = " << y(1) << std::endl;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

	//new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
	
}
