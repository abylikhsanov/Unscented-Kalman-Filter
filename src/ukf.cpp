#include "ukf.h"
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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);

  // H matrix for Lidar update (measurement), as we are using Standard KF for this
  H_ = MatrixXd(2,5);
  H_.fill(0.0);
  H_(0,0) = 1;
  H_(1,1) = 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.1;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/12;

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

  n_aug_ = 7;

  n_x_ = 5;

  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  Xsig = MatrixXd(n_x_,2*n_aug_+1);
  Xsig.fill(0);

  lambda_ = 3-n_aug_;

  time_us_ = 0;

  is_initialized_ = false;

    NIS_laser_ = 0.0;
    NIS_radar_ = 0.0;





}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if(!is_initialized_){
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            // Converting polar to cartesian
            double ro = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            while(phi<M_PI) phi -= 2.*M_PI;
            while(phi>-M_PI) phi += 2.*M_PI;
            double ro_dot = meas_package.raw_measurements_(2);
            double vx = ro_dot * cos(phi);
            double vy = ro_dot * sin(phi);

            x_(0) = ro * cos(phi);
            x_(1) = ro * sin(phi);
            x_(2) = sqrt(vx*vx + vy*vy);
            P_(0,0) = std_radr_;
            P_(1,1) = std_radr_;
        }

        else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
            P_(0,0) = std_laspx_;
            P_(1,1) = std_laspy_;
        }
        P_(2,2) = 20;
        P_(3,3) = 2;
        P_(4,4) = 1;
        cout<<x_<<endl;
        cout<<P_<<endl;
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }


    float delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
      cout<<"Updating Radar"<<endl;
    UpdateRadar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      cout<<"Updating Laser"<<endl;
    UpdateLidar(meas_package);
  }

  cout<<"NIS LASER IS "<<NIS_laser_<<endl;
  cout<<"NIS RADAR IS "<<NIS_radar_<<endl;
  cout<<endl;

}


void UKF::Prediction(double delta_t) {

    cout<<"Prediction"<<endl;


  // Now adding noise vector to the sigma points (Augmented State and Cov. Matrix)
    MatrixXd Xsig_aug(n_aug_,2*n_aug_+1);
    Xsig_aug.fill(0.0);
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.fill(0.0);
    MatrixXd Paug(7,7);
    Paug.fill(0.0);

    x_aug.head(5) = x_;
    Xsig_aug.col(0) = x_aug;
    Paug.topLeftCorner(5,5) = P_;
    Paug(5,5) = pow(std_a_,2);
    Paug(6,6) = pow(std_yawdd_,2);
    MatrixXd L = Paug.llt().matrixL();
    for(int i = 0; i<n_aug_;++i){
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*L.col(i);
        Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_+n_aug_)*L.col(i);
    }


    // Next, we shall predict the sigma points.
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double px = Xsig_aug(0,i);
        double py = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        double temp = yaw+yawd*delta_t;

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.002) {
            px_p = px + v/yawd * ( sin (temp) - sin(yaw));
            py_p = py + v/yawd * ( cos(yaw) - cos(temp) );
        }
        else {
            px_p = px + v*delta_t*cos(yaw);
            py_p = py + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    // Now, we need the Mean State and Covariance matrix, in other words, mean of all sigma points (predicted)
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

    x_.fill(0.0);
    P_.fill(0.0);
    for(int i = 0; i<2*n_aug_+1; ++i){
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
    }


    for(int i = 0; i<2*n_aug_+1; ++i){
        VectorXd diff = Xsig_pred_.col(i)-x_;
        // Need to normalise the angle
        while(diff(3)>M_PI) diff(3)-=2.*M_PI;
        while(diff(3)<-M_PI) diff(3)+=2.*M_PI;
        P_ += weights_(i)*diff*(diff.transpose());
    }




    // After we predicted, now it is turn for an update:


}


void UKF::UpdateLidar(MeasurementPackage meas_package) {

    // Implementing Standard Kalman Filter for the Lidar measurement update;
    VectorXd z = VectorXd(2);
    z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);

    VectorXd z_pred = H_*x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd R(2,2);
    R.fill(0.0);
    R(0,0) = pow(std_laspx_,2);
    R(1,1) = pow(std_laspy_,2);
    MatrixXd S = H_*P_*Ht + R;
    MatrixXd K = P_*Ht*(S.inverse());
    x_ = x_ + K*y;

    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size,x_size);
    P_ = (I-K*H_)*P_;

    NIS_laser_ = (meas_package.raw_measurements_-z_pred).transpose()*S.inverse()*(meas_package.raw_measurements_-z_pred);

}


void UKF::UpdateRadar(MeasurementPackage meas_package) {
// FIRST PART: Converting all the predicted Sigma points to Radar dimensional space


    // Generating pred sigma points (converting Xsig_pred to Zsig_pred)
    MatrixXd Zsig_pred(3,2*n_aug_+1);
    Zsig_pred.fill(0);
    for(int i = 0; i<2*n_aug_+1; ++i){
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v = Xsig_pred_.col(i)(2);
        double yaw = Xsig_pred_.col(i)(3);
        double yaw_d = Xsig_pred_.col(i)(4);
        Zsig_pred.col(i)(0) = sqrt(pow(px,2)+pow(py,2));

        if(Zsig_pred.col(i)(0) > 0.01){
            Zsig_pred.col(i)(1) = atan2(py,px);
            Zsig_pred.col(i)(2) = ((px*cos(yaw))*v + py*sin(yaw)*v)/sqrt(pow(px,2)+pow(py,2));
        }
    }

    // Calculating the z mean (predicted)
    VectorXd z_mean(3);
    z_mean.fill(0.0);
    for(int i = 0; i<2*n_aug_+1;++i){
        z_mean += weights_(i)*Zsig_pred.col(i);
    }

    // Calculating the Predicted covariance in radar space
    MatrixXd S(3,3);
    S.fill(0.0);
    for(int i = 0; i<2*n_aug_+1;++i){
        VectorXd diff(3);
        diff = Zsig_pred.col(i) - z_mean;
        while(diff(1)>M_PI) diff(1)-=2.*M_PI;
        while(diff(1)<-M_PI) diff(1)+=2.*M_PI;
        S += weights_(i)*diff*(diff.transpose());
    }

    // Adding noise R to the predicted State covariance for the Radar
    MatrixXd R(3,3);
    R.fill(0.0);
    R(0,0) = pow(std_radr_,2);
    R(1,1) = pow(std_radphi_,2);
    R(2,2) = pow(std_radrd_,2);

    S += R;

// SECOND PART: Now, the actual update step, after conversion

    // Cross-correlation Matrix.

    MatrixXd T(n_x_,3);
    T.fill(0);
    for(int i = 0; i<2*n_aug_+1; ++i){
        VectorXd diff_z(3);
        diff_z = Zsig_pred.col(i)-z_mean;
        while(diff_z(1)>M_PI) diff_z(1)-=2.*M_PI;
        while(diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;

        VectorXd diff_x = Xsig_pred_.col(i)-x_;
        while(diff_x(3)>M_PI) diff_x(3)-=2.*M_PI;
        while(diff_x(3)<-M_PI) diff_x(3)+=2.*M_PI;

        T += weights_(i)*diff_x*(diff_z.transpose());
    }

    // Kalman Gain Matrix

    MatrixXd K(n_x_,n_x_);
    K = T*(S.inverse());

    // Now it is turn to update the Actual X (state)
    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_diff = z - z_mean;

    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    x_ += K*z_diff;

    P_ = P_ - K*S*(K.transpose());
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}
