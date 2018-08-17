#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  // Declare and set transpose of F
  Eigen::MatrixXd F_t = F_.transpose();

  // KK Predict Funcitons
  x_ = F_ * x_;
  P_ = F_ * P_ * F_t + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  // Declare and set transpose of H
  MatrixXd H_t = H_.transpose();

  // Declare and set I (4x4 Identity Matrix)
  Eigen::MatrixXd I = MatrixXd::Identity(4,4);

  // KM Update Equations
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_t + R_;
  MatrixXd K = P_ * H_t * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * H_)*P_;
}

VectorXd KalmanFilter::Cart2Pol(const VectorXd &x_vec){

  // Extract Position and Velocity values from the state vector
  double px = x_vec[0];
  double py = x_vec[1];
  double vx = x_vec[2];
  double vy = x_vec[3];

  // Covert the cartesian values to polar values
  double rho = sqrt(px*px + py*py);
  double phi = atan2(py,px);

  // Check rho to avoid divide by zero
  if (rho < 0.00001){
    rho = 0.00001;
  }

  double rho_dot = (px*vx + py*vy) / rho;

  //Initialize and set polar state vector
  VectorXd x_pol = VectorXd(3);
  x_pol << rho, phi, rho_dot;

  return x_pol;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Convert state vector to polar coordinates
  VectorXd z_pol = Cart2Pol(x_);
  VectorXd y = z - z_pol;

  // Adjust rho to keep value between -Pi and Pi
  if(y(1) > M_PI){
    y(1) -= 2.0*M_PI;
  } else if (y(1) < -M_PI){
    y(1) += 2.0*M_PI;
  }

  // Declare and set transpose of H
  MatrixXd H_t = H_.transpose();

  // Declare and set I (4x4 Identity Matrix)
  Eigen::MatrixXd I = MatrixXd::Identity(4,4);

  // KM Update Equations
  MatrixXd S = H_ * P_ * H_t + R_;
  MatrixXd K = P_ * H_t * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * H_)*P_;
}
