#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  /**
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estiate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  VectorXd z_pred_polar = GetPredictedInPolarForm();
  
  VectorXd y = z - z_pred_polar;
  
  while ( y(1) > M_PI )
  {
    y(1) -= (2 * M_PI);
  }
  
  while ( y(1) < -M_PI )
  {
    y(1) += (2 * M_PI);
  }
  
  MatrixXd Hj_t = H_.transpose();
  MatrixXd S = H_ * P_ * Hj_t + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Hj_t;
  MatrixXd K = PHt * Si;
  
  //new estiate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::GetPredictedInPolarForm()
{
  VectorXd z_pred_polar = VectorXd(3);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho     = 0;
  float phi     = 0;
  float rho_dot = 0;
  
  //calculate rho
  rho = sqrt(pow(px,2) + pow(py,2));
  
  //calculate phi
  phi = atan2(py,px);
  
  //calculate phi dot
  float rho_dot_temp = (px * vx) + (py * vy);

  if(rho < 0.0001)
  {
    rho_dot = rho_dot_temp / 0.0001;
  }
  else
  {
    rho_dot = rho_dot_temp / rho;
  }
  
  z_pred_polar << rho, phi, rho_dot;
  
  return z_pred_polar;
}
