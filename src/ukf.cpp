#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ <<    0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.45;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.45;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.1;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.1;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.2;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.2;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.2;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // number of states
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // Xsig Prediction
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Xaug
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Initialize lambda
  lambda_ = 3 - n_x_;

  // initialize weights
  weights_= VectorXd(2*n_aug_+1);

  // Initialize NIS
  NIS_radar_ = 0;
  NIS_laser_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double x_loc = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      double y_loc = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
      double v = meas_package.raw_measurements_[2];
      double vx = meas_package.raw_measurements_[2]*cos(meas_package.raw_measurements_[1]);
      double vy = meas_package.raw_measurements_[2]*sin(meas_package.raw_measurements_[1]);
      double psi = atan(vy/vx);
      double psi_dot = 0;

      if (x_loc == 0 && y_loc == 0) {
        x_loc = 1e-100;
        y_loc = 1e-100;
      }

      x_ << x_loc, y_loc, v, psi, psi_dot;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double x_loc = meas_package.raw_measurements_[0];
      double y_loc = meas_package.raw_measurements_[1];

      if (x_loc == 0 && y_loc == 0) {
        x_loc = 1e-100;
        y_loc = 1e-100;
      }

      x_ << x_loc, y_loc, 0, 0, 0;
    }
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction Step
  float dt = (meas_package.timestamp_ - previous_timestamp_)/ 1000000.0;
  Prediction(dt);

  previous_timestamp_ = meas_package.timestamp_;

  // Update Step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)  {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean stateVectorXd x_aug = VectorXd(7);
  VectorXd x_aug = VectorXd(7);

  x_aug << x_, 0,0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug_.col(0) = x_aug;

  for (int i=0; i<n_aug_; i++) {
      Xsig_aug_.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
      Xsig_aug_.col(n_aug_+i+1) = x_aug - sqrt(lambda_+n_aug_)*A.col(i);
  }

  //predict sigma points
  for (int i=0; i < n_aug_*2+1; i++) {

      double px = Xsig_aug_.col(i)(0);
      double py = Xsig_aug_.col(i)(1);
      double v = Xsig_aug_.col(i)(2);
      double theta = Xsig_aug_.col(i)(3);
      double theta_dot = Xsig_aug_.col(i)(4);
      double sig_a = Xsig_aug_.col(i)(5);
      double sig_theta = Xsig_aug_.col(i)(6);

      if (theta_dot > 0) {
          Xsig_pred_.col(i)(0) = px + (v/theta_dot)*(sin(theta+theta_dot*delta_t)-sin(theta))+0.5*delta_t*delta_t*cos(theta)*sig_a;
          Xsig_pred_.col(i)(1) = py + (v/theta_dot)*(-cos(theta+theta_dot*delta_t)+cos(theta))+0.5*delta_t*delta_t*sin(theta)*sig_a;
      } else {
          Xsig_pred_.col(i)(0) = px + (v)*(cos(theta)*delta_t)+0.5*delta_t*delta_t*cos(theta)*sig_a;
          Xsig_pred_.col(i)(1) = py + (v)*(sin(theta)*delta_t)+0.5*delta_t*delta_t*sin(theta)*sig_a;
      }
      Xsig_pred_.col(i)(2) = v+ delta_t*sig_a;
      Xsig_pred_.col(i)(3) =theta + theta_dot*delta_t+0.5*delta_t*delta_t*sig_theta;
      Xsig_pred_.col(i)(4) = theta_dot + sig_theta*delta_t;
  }

  //predict mean and covariance
  //set weights
  weights_.setConstant(1/(2*(lambda_+n_aug_)));
  weights_(0) = lambda_/(lambda_+n_aug_);

  //predict state mean
  x_.setZero();
  for (int i = 0; i<Xsig_pred_.cols(); i++) {
      x_ += weights_(i)*Xsig_pred_.col(i);
  }
  P_.setZero();
  for (int i = 0; i<Xsig_pred_.cols(); i++) {
      P_ += weights_(i)*(Xsig_pred_.col(i)-x_)*(Xsig_pred_.col(i)-x_).transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z = 2;
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  double px;
  double py;

  Zsig.setZero();
  z_pred.setZero();
  MatrixXd R(n_z,n_z);
  R.setZero();
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  //transform sigma points into measurement space
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      px = Xsig_pred_.col(i)(0);
      py = Xsig_pred_.col(i)(1);

      Zsig.col(i) << px, py;
  }

  //calculate mean predicted measurement
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      z_pred += weights_(i)*Zsig.col(i);
  }
  MatrixXd diff_sig_pred;
  MatrixXd diff_sig_pred_t;
  S.setZero();

  //calculate measurement covariance matrix S
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      diff_sig_pred = Zsig.col(i) - z_pred;
      diff_sig_pred_t = diff_sig_pred.transpose();
      S += weights_(i)*diff_sig_pred*diff_sig_pred_t;
  }

  S += R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();

  //calculate cross correlation matrix
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    Tc += weights_(i)*(x_diff)*(z_diff).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  x_+= K*(z-z_pred);
  P_ -= K*S*K.transpose();
  NIS_laser_ = (z-z_pred).transpose()*S.inverse()*(z-z_pred);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  double rho;
  double phi;
  double rho_dot;
  double px;
  double py;
  double v;
  double psi;

  Zsig.setZero();
  z_pred.setZero();
  MatrixXd R(n_z,n_z);
  R.setZero();
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  //transform sigma points into measurement space
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      px = Xsig_pred_.col(i)(0);
      py = Xsig_pred_.col(i)(1);
      v = Xsig_pred_.col(i)(2);
      psi = Xsig_pred_.col(i)(3);

      rho = sqrt(px*px+py*py);
      phi = atan(py/px);
      rho_dot = (px*cos(psi)*v+py*sin(psi)*v)/(rho);
      Zsig.col(i) << rho, phi, rho_dot;
  }

  //calculate mean predicted measurement
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      z_pred += weights_(i)*Zsig.col(i);
  }
  MatrixXd diff_sig_pred;
  MatrixXd diff_sig_pred_t;
  S.setZero();

  //calculate measurement covariance matrix S
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      diff_sig_pred = Zsig.col(i) - z_pred;
      diff_sig_pred_t = diff_sig_pred.transpose();
      S += weights_(i)*diff_sig_pred*diff_sig_pred_t;// + R;
  }
  S += R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();

  //calculate cross correlation matrix
  for (int i=0; i<Xsig_pred_.cols(); i++) {
      //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      Tc += weights_(i)*(x_diff)*(z_diff).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  x_+= K*(z-z_pred);
  P_ -= K*S*K.transpose();

  NIS_radar_ = (z-z_pred).transpose()*S.inverse()*(z-z_pred);
}
