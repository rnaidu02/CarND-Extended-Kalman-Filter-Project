#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;



/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

	// Initialize H_laser_
	H_laser_ << 1, 0, 0, 0,
							0, 1, 0, 0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;
			  
  //the initial transition matrix Q_
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  1, 0, 1, 0,
			  0, 1, 0, 1;

	//set the acceleration noise components for Q Matrix
	noise_ax = 9;
	noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
		float pX, pY, vX, vY;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
			//Get ro (range) from the measurement data
			float rho = measurement_pack.raw_measurements_[0];
			//Get theta (bearing) from the measurement data
			float phi = measurement_pack.raw_measurements_[1];
			//Get ro_dot (radial velocity) from the measurement data
			float rho_dot = measurement_pack.raw_measurements_[2];

			//Get the pX and pY values by converting from polar to cartesian space
			pX = rho * cos(phi);
			pY = rho * sin(phi);
			//Get the vX and vY values by converting from polar to cartesian space
			vX = rho_dot * cos(phi);
			vY = rho_dot * sin(phi);
			
			//Set the velocity components to be 0, as deriving velocy from 
			//the first measurement may not be perfect
			vX = vY = 0;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		  /**
		  Initialize state.
		  */
			//set the state with the initial location and zero velocity
			pX = measurement_pack.raw_measurements_[0];
			pY = measurement_pack.raw_measurements_[1];
			vX = vY = 0;
			//ekf_.x_ << , measurement_pack.raw_measurements_[1], 0, 0;

    }

		// if (pX < leastValue)
		pX = tools.SetMinValues(pX);
		pY = tools.SetMinValues(pY);
		/*
		if (fabs(pY) < MIN_VALUE){
			pY = LEAST_VALUE;
		}
		*/

		//set the state with the initial location and  velocity
		ekf_.x_ << pX, pY, vX, vY;

		//cout << "LOG: x_: " << ekf_.x_ << endl;
		//Set the time time stamp for finding the time difference between samples
		previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

	//1. Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
	//2. Set the process covariance matrix Q
	float t_sq = dt*dt;
	float t_cu = (t_sq*dt)/2;
	float t_fr = (t_sq*t_sq)/4;
	ekf_.Q_(0,0) = t_fr*noise_ax;
	ekf_.Q_(0,2) = t_cu*noise_ax;
	ekf_.Q_(1,1) = t_fr*noise_ay;
	ekf_.Q_(1,3) = t_cu*noise_ay;
	ekf_.Q_(2,0) = t_cu*noise_ax;
	ekf_.Q_(2,2) = t_sq*noise_ax;
	ekf_.Q_(3,1) = t_cu*noise_ay;
	ekf_.Q_(3,3) = t_sq*noise_ay;
	
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
		//Find the Hj from Tools::CalculateJacobian(x)

		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;
	
		//Store rho, theta, rho_top inside z
		VectorXd z(3);
		z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
		ekf_.UpdateEKF(z);

  } else {
    // Laser updates
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		
		//Store px, py inside z
		VectorXd z(2);
		z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
		ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
