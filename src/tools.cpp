#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // Declare RMSE and set initial values to 0
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Check to make sure estimation and ground truth vectors are of same length
  // If not print error and return rmse with value 0, 0, 0, 0.
  if(estimations.size() != ground_truth.size()){
    std::cout << "ERROR in CalculateRMSE(): Estimation and ground truth size mismatch!" << std::endl;
    return rmse;
  }


  for(unsigned int i=0; i<estimations.size(); i++){

    // Difference between estimation and ground truth
    VectorXd diff = estimations[i] - ground_truth[i];

    // Perform element wise squaring of the difference
    diff = diff.array()*diff.array();

    // Add the difference to rmse (summation)
    rmse += diff;
  }

  // Divide summation by number of elements and take element-wise sqrt
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // Extract position and velocity values from x_state;
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Check to make sure px and py are not both zero
  if (abs(px) < 0.00001 && abs(py) < 0.00001){
    // Set px to small value if close to zero
    if(px >= 0){
      px = 0.00001;
    } else {
      px = -0.00001;
    }
    
    // Set py to small value if close to zero
    if(py >= 0){
      py = 0.00001;
    } else {
      py = -0.00001;
    }

    // Print Warning Message
    std::cout << "WARNING in CalculateJacobian: px and py values adjusted to avoid 'Divide By Zero'" << std::endl;
  }

  // Declare Jacobian
  MatrixXd H_j(3,4);

  // Calculate the donominators for the jacobian matrix calculation
  double denom_1 = sqrt(px*px + py*py);
  double denom_2 = px*px + py*py;
  double denom_3 = denom_1 * denom_2;

  // Calculate the Jacobian
  H_j <<                 px/denom_1,                  py/denom_1,           0,           0,
                      -(py/denom_2),                  px/denom_2,           0,           0,
         py*(vx*py - vy*px)/denom_3,  px*(vy*px - vx*py)/denom_3,  px/denom_1,  py/denom_1;

  return H_j;
}


