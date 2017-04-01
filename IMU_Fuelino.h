// IMU Library. Last update: 27 March 2017 
// Davide Cavaliere www.monocilindro.com dadez87-at-gmail.com

#ifndef IMU_Fuelino_h
#define IMU_Fuelino_h

#include <tchar.h>
#include <string>
#include <conio.h>
#include <time.h>
#include <math.h>

// History samples for calculating statistics (average and sigma) and real-time corrections
#define A_HIST_NUM (unsigned int)10 // number of samples for acceleration sensor history
#define G_HIST_NUM (unsigned int)10 // number of samples for gyroscope sensor history

// IMU Class
class IMU_Fuelino_class{

public:
	IMU_Fuelino_class(); // constructor, to initialize data
	void calculation_task(); // This function does the job. It should be called from outside.
	void IMU_Fuelino_class::IMU_init(); // initializes the variables
	
	int evaluate_yaw_rotation_turns(double yaw_now, double* yaw_prev_ptr, int turns_prev);
	double integration_time; // Integration time [s]
	int ax_raw, ay_raw, az_raw; // Acceleration raw data [LSB]
	int gx_raw, gy_raw, gz_raw; // Gyroscope raw data [LSB]
	double a_g_norm; // Acceleration G norm, after sensor offset and gain correction
	double pitch_a; // pitch calculated from accelerometer data
	double roll_a; // roll calculated from accelerometer data
	double a_f_delta_norm; // acceleration vector delta between sample and average
	double filter_req_pitchroll; // filter constant for yaw and pitch (depends on how acceleration data is reliable)
	double q_f_out[4]; // Output Quaternion
	double yaw_out; // Output [rad]
	double pitch_out; // Output [rad]
	double roll_out; // Output [rad]
	int turns_clockwise; // how many turns clockwise have been done
	
private:
	void process_raw_data(); // fixes the sign and scale, and arranges data into float format
	void calculate_rotation_matrix_elements(double yaw_in, double pitch_in, double roll_in, double M[3][3]);
	void calculate_quaternion_elements(double R_a[3][3], double* q_a);
	void calculate_pitchroll_from_acceleration();
	void calculate_yawpitchroll_from_quaternions(double* q_in, double* yaw_out, double* pitch_out, double* roll_out);
	void integrate_quaternion_g(double* q_in, double* q_out);
	void normalize_quaternion(double* q_in, double* q_out);
	void statistics_calculation();
	void fuse_yawpitchroll_block();
	void quaternion_product(double* q_in1, double* q_in2, double* q_out);
	void Madgwick_updateIMU(double* q_in, double ax, double ay, double az, double gx, double gy, double gz, double int_time, double madg_beta, double* q_out);
	void vector_rotation_correction(double *vet_in, double *vet_out, bool correction); // orientation correction

	double a_f_cor[3]; // Acceleration [g] after correction
	double g_f_cor[3]; // Gyroscope [rad/s]
	double a_f_cor_fil[3]; // Acceleration [g] after correction, filtered
	double yaw_old; // used to estimate rotation turns
	double yaw_fuse, pitch_fuse, roll_fuse; // yaw, pitch, roll calculated from fusion
	double M_rot[3][3]; // rotation matrix calculated from yaw, pitch, roll
	double M_cor[3][3]; // rotation matrix for input signal correction
	double q_M_cor[4]; // quaternion or M correction
	double q_int_in[4]; // quaternion (before integration)
	double q_int_out[4]; // quaternion (after integration)

	// Offsets estimation variables for gyroscope data
	double gx_f_offset;
	double gy_f_offset;
	double gz_f_offset;

	// History variables
	double a_mod_samples_hist[A_HIST_NUM];
	double gx_samples_hist[G_HIST_NUM];
	double gy_samples_hist[G_HIST_NUM];
	double gz_samples_hist[G_HIST_NUM];
	unsigned int a_samples_index;
	unsigned int g_samples_index;

};

extern IMU_Fuelino_class IMU_Fuelino;

#endif