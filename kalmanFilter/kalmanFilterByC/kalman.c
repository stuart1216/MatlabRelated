#include "Kalman.h"
/**
 *@init_x：initial value of measurement
 *@init_p：The initial value of the variance of the error of the estimated value of the posterior state
 */
void kalmanFilter_init(KalmanStructTypedef *kalmanFilter, float init_x, float init_p,float predict_q,float newMeasured_q)
{
    kalmanFilter->x = init_x;// middle value
    kalmanFilter->p = init_p;// not equals to 0
    kalmanFilter->A = 1;
    kalmanFilter->H = 1;
    kalmanFilter->q = predict_q;//Prediction (process) noise variance, realated to convergence rate
    kalmanFilter->r = newMeasured_q;//measurement noise variance, can be got from experiment
}

/**
 *@function: - 
 *@kalmanFilter:
 *@newMeasured
 */
float kalmanFilter_filter(KalmanStructTypedef *kalmanFilter, float newMeasured)
{
    /* Predict */
    kalmanFilter->x = kalmanFilter->A * kalmanFilter->x;
    kalmanFilter->p = kalmanFilter->A * kalmanFilter->A * kalmanFilter->p + kalmanFilter->q;  

    /* Correct */
    kalmanFilter->gain = kalmanFilter->p * kalmanFilter->H / (kalmanFilter->p * kalmanFilter->H * kalmanFilter->H + kalmanFilter->r);
    kalmanFilter->x = kalmanFilter->x + kalmanFilter->gain * (newMeasured - kalmanFilter->H * kalmanFilter->x);
    kalmanFilter->p = (1 - kalmanFilter->gain * kalmanFilter->H) * kalmanFilter->p;

    return kalmanFilter->x;
}