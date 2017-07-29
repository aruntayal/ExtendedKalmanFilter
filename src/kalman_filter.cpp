#include "kalman_filter.h"
#include<iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::ComputeEstimateAndUpdate(const VectorXd& y)
{
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);

  //Update Process covarience matrix
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::Update(const VectorXd &z) {

  //Get Measurement from sensor
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  ComputeEstimateAndUpdate(y);  
}


MatrixXd GetRadarMeas(const VectorXd& x_state) {

    VectorXd z_radar;
    z_radar = VectorXd(3);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  

  //check division by zero
    z_radar(0) = sqrt(pow(px,2)+pow(py,2));

    if (sqrt(pow(px,2)+pow(py,2))<.0001){
        z_radar(2) = 0;
    }
    else{
        z_radar(2) = (px*vx+py*vy)/sqrt(pow(px,2)+pow(py,2));
    }


  z_radar(1) = atan2(py,px);

  return z_radar;




}
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  const float PI = 3.14;
  VectorXd z_pred = GetRadarMeas(x_); //H_ * x_;
  VectorXd y = z - z_pred;
  float angle = y(1);

  if ((angle > PI) || (angle < -1 * PI))
  { 
    while(angle > PI)
    {
      angle -= 2 * PI ;
    }

    while(angle < (-1 * PI))
    {
      angle += 2 * PI ;
    }

    y(1) = angle;
  }

  ComputeEstimateAndUpdate(y); 

}


