#ifndef _Kalman_H_
#define _Kalman_H_
typedef struct {
    float x;  // state parameters of system
    float A;  // x(n)=A*x(n-1)+u(n),u(n)~N(0,q)
    float H;  // z(n)=H*x(n)+w(n),w(n)~N(0,r)
    float q;  // estimation noise covariance
    float r;  // measurement noise covariance
    float p;  // Estimation error covariance
    float gain;// kalman gain
} KalmanStructTypedef;
void kalmanFilter_init(KalmanStructTypedef *kalmanFilter, float init_x, float init_p,float predict_q,float newMeasured_q);
float kalmanFilter_filter(KalmanStructTypedef *kalmanFilter, float newMeasured);
#endif