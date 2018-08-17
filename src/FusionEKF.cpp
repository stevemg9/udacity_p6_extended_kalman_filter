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

  // Initializing Matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Measurement Noise Covariance Matrix - LiDAR
  R_laser_ << 0.0225,       0,
                   0,  0.0225;

  // Measurement Noise Covariance Matrix - Radar
  R_radar_ << 0.09,       0,     0,
                 0,  0.0009,     0,
                 0,       0,  0.09;

  // Measurement Matrix - LiDAR
  H_laser_ <<  1, 0, 0, 0,
               0, 1, 0, 0;

  // State Covariance Matrix
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ <<  1,  0,     0,     0,
              0,  1,     0,     0,
              0,  0,  1000,     0,
              0,  0,     0,  1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  cout << "FUSION EKF PROCESS MEASUREMENT" << std::endl;
  if (!is_initialized_) {

    // Initialize x_ for ekf
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      // Read in initial Radar measurements
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      // Convert rho, phi and rho_dot to x, y vx and vy
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      //Set Initial State
      ekf_.x_ <<   x,
                   y,
                  vx,
                  vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      //For LiDAR measurements, directly set x and y and set vx and vy to zero
      double x = measurement_pack.raw_measurements_[0];
      double y = measurement_pack.raw_measurements_[1];

      //Set Initial State
      ekf_.x_ <<  x,
                  y,
                  0,
                  0;
    }

    // Update Timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // Set Initialization flag
    is_initialized_ = true;

    // Return without Update or Predict
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Calculate delta t (time between last and current measurement)
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  //Initialize and set State Transition Matrix F
  Eigen::MatrixXd F = MatrixXd(4,4);
  F <<  1,  0, dt,  0,
        0,  1,  0, dt,
        0,  0,  1,  0,
        0,  0,  0,  1;

  //Calculate delta t values for use in Process Noise Covariance Matrix
  double dt_2 = dt*dt;
  double dt_3d2 = (dt_2*dt)/2.0;
  double dt_4d4 = (dt_2*dt*dt)/4.0;
  double noise_ax = 9;
  double noise_ay = 9;

  // Initialize and set Process Noise Covariance Matrix Q
  Eigen::MatrixXd Q = MatrixXd(4,4);
  Q <<  dt_4d4*noise_ax,                 0,   dt_3d2*noise_ax,                 0,
                      0,   dt_4d4*noise_ay,                 0,   dt_3d2*noise_ay,
        dt_3d2*noise_ax,                 0,     dt_2*noise_ax,                 0,
                      0,   dt_3d2*noise_ay,                 0,     dt_2*noise_ay;

  // Assign F and Q to ekf member variables
  ekf_.F_ = F;
  ekf_.Q_ = Q;
  ekf_.Predict();



  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar Update------------------------------------------------
    
    // Calculate Jacobian
    Hj_ = tools.CalculateJacobian(ekf_.x_);

    // Assign Jacobian and Measurement Noise Covariance Matrices to ekf
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    // Call EKF Update Function
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // LiDAR Update-----------------------------------------------

    // Assign Measurement and Measurement Noise Covariance Matrices to ekf
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    //Call Update Function
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // Print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
