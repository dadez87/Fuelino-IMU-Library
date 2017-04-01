#ifndef IMU_Fuelino_cal_h
#define IMU_Fuelino_cal_h

// Filter type (1=Cavaliere, 0=Madgwick)
#define CAVALIERE_FILTER_ENABLE 0

// Constant values which should not be touched. 
// Do not touch these defines, unless you know what you are doing.
#define RAD_TO_DEG (double)57.295777951 // convert radiants to degree
#define DEG_TO_RAD (double)1.0/RAD_TO_DEG // convert radiants to degree
#define M_PI (double)3.14159  // pi greco value
#define GYRO_SENSITIVITY (double)250.0/(double)32768.0 // 250 deg/s is the full scale, which corresponds to 2^15 = 32768
#define GYRO_TO_RAD (double)DEG_TO_RAD*GYRO_SENSITIVITY // needed to tranform sensor raw data into rad per second (3.14/180)*(1/32.8). Gyro range is 32.8 LSB/(deg/s) (Arduino range = 1000)
#define VERBOSE_MODE 0 // Flag to allow prints messages on screen

// In case of MPU-6050, each signal (3x acceleration, 3x gyroscope) corresponds to a 16bit signed integer (signed short)

// Accelerometer Calibration Data (360 degrees rotation test)
#define ACC_X_MIN (int)-15836 //-16384
#define ACC_X_MAX (int)16848 //16384
#define ACC_Y_MIN (int)-16429 //-16384
#define ACC_Y_MAX (int)16399 //16384
#define ACC_Z_MIN (int)-17617 //-16384
#define ACC_Z_MAX (int)15737 //16384
#define ACC_X_OFFSET (ACC_X_MAX + ACC_X_MIN)/2 // X acc offset
#define ACC_Y_OFFSET (ACC_Y_MAX + ACC_Y_MIN)/2 // Y acc offset
#define ACC_Z_OFFSET (ACC_Z_MAX + ACC_Z_MIN)/2 // Z acc offset
#define ACC_X_ONEG (ACC_X_MAX - ACC_X_MIN)/2 // X acc value for 1g = 9.81m/s2
#define ACC_Y_ONEG (ACC_Y_MAX - ACC_Y_MIN)/2 // Y acc value for 1g = 9.81m/s2
#define ACC_Z_ONEG (ACC_Z_MAX - ACC_Z_MIN)/2 // Z acc value for 1g = 9.81m/s2

// Gyroscope Calibration Data (steady state output signal)
#define GYR_X_OFFSET (int)-144
#define GYR_Y_OFFSET (int)239
#define GYR_Z_OFFSET (int)-12

// Sensor Orientation Data (steady state) - Sensor offset when in horizontal position
#define PITCH_OFFSET 0//(double)-14.75*DEG_TO_RAD
#define ROLL_OFFSET 0//(double)9.8*DEG_TO_RAD

// Offset compensation in realtime, for gyroscope signals
#define G_OFFSET_CORR_ACT_FLAG 0 // Enables Gyroscope raw data offset correction algorithm
#define G_STEADY_STATE_SIGMA_THR (double)50.0 // Threshold for gyroscope value: division between sigma and average

// Acceleration Pitch and Roll angles: filtering and saturations
#define A_F_COR_IIR_FILTER_COEFF (double)0.1 // filter coefficient for IIR filter
#define MAX_PITCH_VAR_ACC (double)250.0*(DEG_TO_RAD) // maximum variation acceptable in 1 second (acceleration sensor filtering)
#define MAX_ROLL_VAR_ACC (double)250.0*(DEG_TO_RAD) // maximum variation acceptable in 1 second (acceleration sensor filtering)
#define MAX_PITCH (double)85.0*(DEG_TO_RAD) // maximum value (acceleration sensor filtering)
#define MAX_ROLL (double)85.0*(DEG_TO_RAD) // maximum value (acceleration sensor filtering)

// Fusion constants, used to fuse the pitch and roll angles calculated from Accelerometer and Gyroscope data.
#define FUSION_FILTER_MAX (double)0.1 // when acceleration pitch/roll has high priority

#endif