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
  H_radar_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1,0,0,0,
              0,1,0,0;

  H_radar_ << 1,0,0,0,
              0,1,0,0,
              0,0,1,0;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1,0,0,0,
              0,1,0,0,
              0,0,1000,0,
              0,0,0,1000;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.P_ << 1,0,0,0,
              0,1,0,0,
              0,0,1,0,
              0,0,0,1; // some random initialization here, fixed it later.


  noise_ax = 9.0 ;
  noise_ay = 9.0;


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
  
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    //Check if  the sensor is RADAR then convert polar coordinates to cartesian.
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
         float rho = measurement_pack.raw_measurements_(0);
      float th =  measurement_pack.raw_measurements_(1);
      float drho =  measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = rho*cos(th);
      ekf_.x_(1) = rho*sin(th);
      ekf_.x_(2) = drho*cos(th); // Approximate value as dth is not known
      ekf_.x_(3) = 0; // Approximate value as dth is not known
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.x_(0) =measurement_pack.raw_measurements_(0);
      ekf_.x_(1) =measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0; // Approximate value of 0
      ekf_.x_(3) = 0; // Approximate value of 0
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;


  cout << dt << endl;
  ekf_.F_ <<  1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;

    ekf_.Q_ << (pow(dt,4)/4*noise_ax), 0, (pow(dt,3)/2*noise_ax), 0,
                  0, (pow(dt,4)/4*noise_ay), 0, (pow(dt,3)/2*noise_ay),
                  (pow(dt,3)/2*noise_ax), 0, pow(dt,2)*noise_ax, 0,
0, (pow(dt,3)/2*noise_ay), 0, pow(dt,2)*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/



  //Perform update according to sensor type
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
     H_radar_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = H_radar_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
