#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

static inline void NormalizeAngle(double &angle) {
  while (angle > M_PI) angle -= 2. * M_PI;
  while (angle < -M_PI) angle += 2. * M_PI;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();

 
  std_a_ = 1.5;      // m/s^2
  std_yawdd_ = 0.5;  // rad/s^2

  

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

 

  is_initialized_ = false;
  time_us_ = 0;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // Weights for sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); ++i) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  //  initialize NIS
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // Initialization
 
  if (!is_initialized_) {
    // Initialize state x_ from first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // lidar provides px, py
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      x_ << px, py, 0.0, 0.0, 0.0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // radar provides rho, phi, rho_dot
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      double px = rho * std::cos(phi);
      double py = rho * std::sin(phi);

      // We can initialize speed magnitude with rho_dot (not perfect, but helps)
      double v = std::fabs(rho_dot);

      x_ << px, py, v, 0.0, 0.0;
    }

    // Initialize covariance (reasonable)
    P_.setIdentity();
    P_(0,0) = 1.0;
    P_(1,1) = 1.0;
    P_(2,2) = 1.0;
    P_(3,3) = 1.0;
    P_(4,4) = 1.0;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction
 
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Safety: ignore extremely small dt
  if (dt > 0.0) {
    Prediction(dt);
  }

 
  // Update
  
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0.0;       // nu_a
  x_aug(n_x_ + 1) = 0.0;   // nu_yawdd

  // Create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // Square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;

  double sqrt_lambda_n_aug = std::sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_lambda_n_aug * L.col(i);
  }

  // Predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v   = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd= Xsig_aug(4, i);

    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p;

    // CTRV model
    if (std::fabs(yawd) > 1e-6) {
      px_p = p_x + (v / yawd) * (std::sin(yaw + yawd * delta_t) - std::sin(yaw));
      py_p = p_y + (v / yawd) * (-std::cos(yaw + yawd * delta_t) + std::cos(yaw));
    } else {
      px_p = p_x + v * delta_t * std::cos(yaw);
      py_p = p_y + v * delta_t * std::sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add process noise
    px_p  += 0.5 * nu_a * delta_t * delta_t * std::cos(yaw);
    py_p  += 0.5 * nu_a * delta_t * delta_t * std::sin(yaw);
    v_p   += nu_a * delta_t;

    yaw_p  += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));  // yaw angle
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Lidar measurement dimension
  int n_z = 2;

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    Zsig(0, i) = Xsig_pred_(0, i); // px
    Zsig(1, i) = Xsig_pred_(1, i); // py
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  S += R;

  // Cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // Update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // NIS (optional)
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Radar measurement dimension
  int n_z = 3;

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px  = Xsig_pred_(0, i);
    double py  = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double vx = v * std::cos(yaw);
    double vy = v * std::sin(yaw);

    double rho = std::sqrt(px * px + py * py);
    if (rho < 1e-6) rho = 1e-6;

    Zsig(0, i) = rho;
    Zsig(1, i) = std::atan2(py, px);
    Zsig(2, i) = (px * vx + py * vy) / rho;
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Measurement covariance S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  S += R;

  // Cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));

  // Update state and covariance
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // NIS (optional)
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
