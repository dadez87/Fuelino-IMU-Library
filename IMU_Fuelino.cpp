// IMU Library. Last update: 27 March 2017 
// Davide Cavaliere www.monocilindro.com dadez87-at-gmail.com

#include "stdafx.h"
#include "IMU_Fuelino.h" // Class and types definitions
#include "IMU_Fuelino_cal.h" // Calibrateable values

IMU_Fuelino_class IMU_Fuelino;

// Returns the sign of the input
double segno(double x) {
	return (x >= 0.0) ? +1.0 : -1.0;
}

// Vectors product.
double vectors_dot_prod(double *x, double *y){
	double res = 0.0;
	for (unsigned int i = 0; i < 3; i++)
	{
		res += x[i] * y[i];
	}
	return res;
}

// Matrix*Vector product.
void IMU_Fuelino_class::vector_rotation_correction(double *vet_in, double *vet_out, bool correction){
	for (unsigned int i = 0; i < 3; i++){
		if (correction == true) { // correction requested
			double vec_tmp[3];
			vec_tmp[0] = M_cor[i][0];
			vec_tmp[1] = M_cor[i][1];
			vec_tmp[2] = M_cor[i][2];
			vet_out[i] = vectors_dot_prod(vec_tmp, vet_in);
		}
		else { // just bypass
			vet_out[i] = vet_in[i]; // bypass
		}
	}
}

// Madgwick filter modified algorithm. Uses double (64 bit) precision calculation. There should be no problem on fast computers.
void IMU_Fuelino_class::Madgwick_updateIMU(double* q_in, double ax, double ay, double az, double gx, double gy, double gz, double int_time, double madg_beta, double* q_out) {

	// Buffering (Input)
	double q0 = q_in[0];
	double q1 = q_in[1];
	double q2 = q_in[2];
	double q3 = q_in[3];

	// Calculation
	double Norm;
	double s0, s1, s2, s3;
	double qDot1, qDot2, qDot3, qDot4;
	double _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2, _8q1, _8q2, q0q0, q1q1, q2q2, q3q3;
	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx);
	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if (!((ax == 0.0) && (ay == 0.0) && (az == 0.0))) {
		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0 * q0;
		_2q1 = 2.0 * q1;
		_2q2 = 2.0 * q2;
		_2q3 = 2.0 * q3;
		_4q0 = 4.0 * q0;
		_4q1 = 4.0 * q1;
		_4q2 = 4.0 * q2;
		_8q1 = 8.0 * q1;
		_8q2 = 8.0 * q2;
		q0q0 = q0 * q0;
		q1q1 = q1 * q1;
		q2q2 = q2 * q2;
		q3q3 = q3 * q3;
		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0 * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0 * q1q1 * q3 - _2q1 * ax + 4.0 * q2q2 * q3 - _2q2 * ay;
		Norm = sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 /= Norm;
		s1 /= Norm;
		s2 /= Norm;
		s3 /= Norm;
		// Apply feedback step
		qDot1 -= madg_beta * s0;
		qDot2 -= madg_beta * s1;
		qDot3 -= madg_beta * s2;
		qDot4 -= madg_beta * s3;
	}
	// Integrate rate of change of quaternion to yield quaternion
	q0 += qDot1 * int_time;
	q1 += qDot2 * int_time;
	q2 += qDot3 * int_time;
	q3 += qDot4 * int_time;
	// Normalise quaternion
	Norm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 /= Norm;
	q1 /= Norm;
	q2 /= Norm;
	q3 /= Norm;

	// Posting (Output)
	q_out[0] = q0;
	q_out[1] = q1;
	q_out[2] = q2;
	q_out[3] = q3;
}

// Calculates the quaternion product: q_out = q_in1 * q_in2
void IMU_Fuelino_class::quaternion_product(double* q_in1_tmp, double* q_in2, double* q_out) {
	double q_in1[4]; // buffer for input
	q_in1[0] = q_in1_tmp[0];
	q_in1[1] = q_in1_tmp[1];
	q_in1[2] = q_in1_tmp[2];
	q_in1[3] = q_in1_tmp[3];
	q_out[0] = q_in1[0] * q_in2[0] - q_in1[1] * q_in2[1] - q_in1[2] * q_in2[2] - q_in1[3] * q_in2[3];
	q_out[1] = q_in1[0] * q_in2[1] + q_in1[1] * q_in2[0] + q_in1[2] * q_in2[3] - q_in1[3] * q_in2[2];
	q_out[2] = q_in1[0] * q_in2[2] - q_in1[1] * q_in2[3] + q_in1[2] * q_in2[0] + q_in1[3] * q_in2[1];
	q_out[3] = q_in1[0] * q_in2[3] + q_in1[1] * q_in2[2] - q_in1[2] * q_in2[1] + q_in1[3] * q_in2[0];
	normalize_quaternion(q_out, q_out);
}

// Calculates average and sigma of acceleration and gyroscope raw data
void IMU_Fuelino_class::statistics_calculation() {

	// Saves current acceleration modulus and gyroscope samples in the history arrays
	a_mod_samples_hist[a_samples_index] = a_g_norm; // g
	a_samples_index++; // increase array number
	if (a_samples_index >= A_HIST_NUM) a_samples_index = 0;
	gx_samples_hist[g_samples_index] = g_f_cor[0]; // rad/s
	gy_samples_hist[g_samples_index] = g_f_cor[1]; // rad/s
	gz_samples_hist[g_samples_index] = g_f_cor[2]; // rad/s
	g_samples_index++; // increase array number
	if (g_samples_index >= G_HIST_NUM) g_samples_index = 0;

	// Calculates acceleration array average and sigma
	double average_acc_hist = 0.0; // average
	double sigma_acc_hist = 0.0; // sigma
	for (unsigned int i = 0; i < A_HIST_NUM; i++) {
		average_acc_hist += a_mod_samples_hist[i];
	}
	average_acc_hist = average_acc_hist / A_HIST_NUM; // average calculated
	for (unsigned int i = 0; i < A_HIST_NUM; i++) {
		sigma_acc_hist += pow((a_mod_samples_hist[i] - average_acc_hist), 2);
	}
	sigma_acc_hist = sigma_acc_hist / A_HIST_NUM;
	sigma_acc_hist = sqrt(sigma_acc_hist); // sigma calculated

	// Calculates gyroscope arrays average and sigma
	double average_gx_hist = 0.0;
	double average_gy_hist = 0.0;
	double average_gz_hist = 0.0;
	double sigma_gx_hist = 0.0;
	double sigma_gy_hist = 0.0;
	double sigma_gz_hist = 0.0;
	for (unsigned int i = 0; i < G_HIST_NUM; i++) {
		average_gx_hist += gx_samples_hist[i];
		average_gy_hist += gy_samples_hist[i];
		average_gz_hist += gz_samples_hist[i];
	}
	average_gx_hist /= G_HIST_NUM;
	average_gy_hist /= G_HIST_NUM;
	average_gz_hist /= G_HIST_NUM;
	for (unsigned int i = 0; i < G_HIST_NUM; i++) {
		sigma_gx_hist += pow(gx_samples_hist[i] - average_gx_hist, 2);
		sigma_gy_hist += pow(gy_samples_hist[i] - average_gy_hist, 2);
		sigma_gz_hist += pow(gz_samples_hist[i] - average_gz_hist, 2);
	}
	sigma_gx_hist = sigma_gx_hist / G_HIST_NUM; // average calculated
	sigma_gy_hist = sigma_gy_hist / G_HIST_NUM;
	sigma_gz_hist = sigma_gz_hist / G_HIST_NUM;
	sigma_gx_hist = sqrt(sigma_gx_hist); // sigma calculated
	sigma_gy_hist = sqrt(sigma_gy_hist);
	sigma_gz_hist = sqrt(sigma_gz_hist);

	// Evaluate if the sensor is moving or not
	//double sigma_gx_relative, sigma_gy_relative, sigma_gz_relative;
	if ((average_gx_hist != 0.0) && (average_gy_hist != 0.0) && (average_gz_hist != 0.0)) { // checks to avoid division by zero
		if ((abs(sigma_gz_hist) <= G_STEADY_STATE_SIGMA_THR) && (abs(sigma_gz_hist) <= G_STEADY_STATE_SIGMA_THR) && (abs(sigma_gz_hist) <= G_STEADY_STATE_SIGMA_THR)){
			gx_f_offset = GYRO_TO_RAD*average_gx_hist;
			gy_f_offset = GYRO_TO_RAD*average_gy_hist;
			gz_f_offset = GYRO_TO_RAD*average_gz_hist;
		}
	}

	if (VERBOSE_MODE) printf("Acceleration modulo. Average= %f. Sigma= %f. Now= %f. Filter= %f.\n", average_acc_hist, sigma_acc_hist, a_g_norm, filter_req_pitchroll);
	if (VERBOSE_MODE) printf("gx_off=%f gy_off=%f gz_off=%f gx_s=%f gy_s=%f gz_s=%f\n", gx_f_offset, gy_f_offset, gz_f_offset, sigma_gx_hist, sigma_gy_hist, sigma_gz_hist);

}

// Normalizes (module=1) a quaternion
void IMU_Fuelino_class::normalize_quaternion(double* q_in, double* q_out) {
	double Norm = 1 / sqrt(q_in[0] * q_in[0] + q_in[1] * q_in[1] + q_in[2] * q_in[2] + q_in[3] * q_in[3]);
	q_out[0] = q_in[0] / Norm;
	q_out[1] = q_in[1] / Norm;
	q_out[2] = q_in[2] / Norm;
	q_out[3] = q_in[3] / Norm;
}

// Integrates the quaternion based on gyroscope angle rate signals
void IMU_Fuelino_class::integrate_quaternion_g(double* q_in, double* q_out) {
	double wx = g_f_cor[0]; // speed around X axis
	double wy = g_f_cor[1]; // speed around Y axis
	double wz = g_f_cor[2]; // speed around Z axis
	double wm = sqrt(wx*wx + wy*wy + wz*wz); // modulus of speed vector [rad/s]
	double k = (double)1.0; // initialization (worst case)
	double w0 = (double)1.0; // initialization (worst case)
	double qOmega[4]; // The following quaternion is the present rotation
	if (wm != (double)0.0) { // to avoid division by zero in case of no rotation
		double angolo = (wm * integration_time) / (double)2.0; // speed integrated in time gives rotation angle [rad]. Need to divide by 2 (quaternion)
		w0 = cos(angolo); // real part of the quaternion
		double sin_angolo_meno = -sin(angolo); // sinus of rotation angle, negated
		wx *= sin_angolo_meno; // immaginary part (i)
		wy *= sin_angolo_meno; // immaginary part (j)
		wz *= sin_angolo_meno; // immaginary part (k)
		k = sqrt(((double)1.0 - w0*w0) / (wx*wx + wy*wy + wz*wz)); // needed to normalize the quaternion Search k to have modulus = 1. q = w0 + k*wx + k*wy + k*wz
	}
	qOmega[0] = w0;
	qOmega[1] = wx*k;
	qOmega[2] = wy*k;
	qOmega[3] = wz*k;
	quaternion_product(qOmega, q_in, q_out); // Rotates the "q_in" of "qOmega", and outputs "q_out"
}

// Calculates rotation matrix elements values, starting from Tait-Bryan angles (yaw, pitch, roll)
void IMU_Fuelino_class::calculate_rotation_matrix_elements(double yaw_in, double pitch_in, double roll_in, double R_a[3][3]) {

	double sin_psi_a = sin(yaw_in); // yaw (psi)
	double cos_psi_a = cos(yaw_in);
	double sin_teta_a = sin(pitch_in); // pitch (teta)
	double cos_teta_a = cos(pitch_in);
	double sin_phi_a = sin(roll_in); // roll (phi)
	double cos_phi_a = cos(roll_in);

	// Calculates third column elements (Z axis), which require just the knowledge of pitch and roll angles
	R_a[0][2] = -sin_teta_a;
	R_a[1][2] = sin_phi_a*cos_teta_a;
	R_a[2][2] = cos_phi_a*cos_teta_a;

	// Calculates first and second column elements, which require the knowledge of yaw angle (psi)
	R_a[0][0] = cos_teta_a*cos_psi_a;
	R_a[0][1] = cos_teta_a*sin_psi_a;
	R_a[1][0] = -cos_phi_a*sin_psi_a + sin_phi_a*sin_teta_a*cos_psi_a;
	R_a[1][1] = cos_phi_a*cos_psi_a + sin_phi_a*sin_teta_a*sin_psi_a;
	R_a[2][0] = sin_phi_a*sin_psi_a + cos_phi_a*sin_teta_a*cos_psi_a;
	R_a[2][1] = -sin_phi_a*cos_psi_a + cos_phi_a*sin_teta_a*sin_psi_a;
}

// Calculates quaternion elements from a rotation matrix (Standard and simple method), only if the division is possible (!=0.0)
void IMU_Fuelino_class::calculate_quaternion_elements(double R_a[3][3], double* q_a) {
	double temp_q_a = (double)0.5 * sqrt((double)1.0 + R_a[0][0] + R_a[1][1] + R_a[2][2]);
	if (temp_q_a != (double)0.0) { // this value is used for division
		q_a[0] = temp_q_a;
		q_a[1] = (R_a[2][1] - R_a[1][2]) / ((double)4.0 * q_a[0]);
		q_a[2] = (R_a[0][2] - R_a[2][0]) / ((double)4.0 * q_a[0]);
		q_a[3] = (R_a[1][0] - R_a[0][1]) / ((double)4.0 * q_a[0]);
		normalize_quaternion(q_a, q_a); // modulus becomes 1
	}
	//printf("%f %f %f %f\n\n", q_a[0], q_a[1], q_a[2], q_a[3]);
}

// Calculates Yaw Pitch Roll angles from quaternion
void IMU_Fuelino_class::calculate_yawpitchroll_from_quaternions(double* q_in, double* yaw_out, double* pitch_out, double* roll_out) {
	*yaw_out = atan2((2.0 * q_in[1] * q_in[2] - 2.0 * q_in[0] * q_in[3]), (2.0 * q_in[0] * q_in[0] + 2.0 * q_in[1] * q_in[1]) - 1.0);
	*pitch_out = (-1.0) * asin((2.0 * q_in[1] * q_in[3] + 2.0 * q_in[0] * q_in[2]));
	*roll_out = atan2((2.0 * q_in[2] * q_in[3] - 2.0 * q_in[0] * q_in[1]), (2.0 * q_in[0] * q_in[0] + 2.0 * q_in[3] * q_in[3]) - 1.0);
	if (VERBOSE_MODE) printf("Yawt=%f Pitch= %f Roll=%f\n", RAD_TO_DEG * (*yaw_out), RAD_TO_DEG * (*pitch_out), RAD_TO_DEG * (*roll_out));
}

// Correct signs, scales inputs, and transforms data into floating point format
void IMU_Fuelino_class::process_raw_data() {

	// Correct acceleration and gyroscope offsets and gain (measured during offline calibration)
	double ax_g = ((double)(ax_raw - ACC_X_OFFSET)) / ((double)ACC_X_ONEG); // g, corrected value
	double ay_g = ((double)(ay_raw - ACC_Y_OFFSET)) / ((double)ACC_Y_ONEG); // g, corrected value
	double az_g = ((double)(az_raw - ACC_Z_OFFSET)) / ((double)ACC_Z_ONEG); // g, corrected value
	a_g_norm = sqrt(ax_g*ax_g + ay_g*ay_g + az_g*az_g); // g (total acceleration, norm)
	
	// Normalizes acceleration data inside [-1, +1] range
	double a_f[3]; // acceleration vector before correction, normalized
	a_f[0] = ax_g / a_g_norm; // normalized (0.0..1.0)
	a_f[1] = ay_g / a_g_norm; // normalized (0.0..1.0)
	a_f[2] = az_g / a_g_norm; // normalized (0.0..1.0)

	// Transforms gyroscope data from raw data into rad/s
	int gx = gx_raw - GYR_X_OFFSET; // offset removal
	int gy = gy_raw - GYR_Y_OFFSET; // offset removal
	int gz = gz_raw - GYR_Z_OFFSET; // offset removal
	double g_f[3]; // gyroscope vector before correction
	g_f[0] = GYRO_TO_RAD*(double)gx; // physical value (rad/s)
	g_f[1] = GYRO_TO_RAD*(double)gy; // physical value (rad/s)
	g_f[2] = GYRO_TO_RAD*(double)gz; // physical value (rad/s)

	// Correct signs (some axis have negative sign compared on standard convention)
	//ax = -ax;
	//gy = -gy;
	//gz = -gz;
	a_f[2] = -a_f[2]; // AZ - Fuelino Proto3 on CBR125R - MPU6050
	g_f[0] = -g_f[0]; // GX - Fuelino Proto3 on CBR125R - MPU6050
	g_f[1] = -g_f[1]; // GY - Fuelino Proto3 on CBR125R - MPU6050

	// Correct sensor alignment
	vector_rotation_correction(a_f, a_f_cor, false); // accelerometer correction (to correct sensor position rotation)
	vector_rotation_correction(g_f, g_f_cor, false); // gyroscope correction (to correct sensor position rotation)

	// Correct gyroscope offset (based on estimated offset values)
	if (G_OFFSET_CORR_ACT_FLAG) {
		g_f_cor[0] -= gx_f_offset; // Standard value is 0.0
		g_f_cor[1] -= gy_f_offset; // Standard value is 0.0
		g_f_cor[2] -= gz_f_offset; // Standard value is 0.0
	}
}

// Calculates pitch and roll angles from acceleration data, and determines the acceleration data reliability coefficient for fusion
void IMU_Fuelino_class::calculate_pitchroll_from_acceleration() {
	
	// Filter acceleration signal (1st order IIR filter)
	for (unsigned int i = 0; i < 3; i++) { // 3 axis acceleration
		a_f_cor_fil[i] = A_F_COR_IIR_FILTER_COEFF*a_f_cor[i] + ((double)1.0 - A_F_COR_IIR_FILTER_COEFF)*a_f_cor_fil[i]; // IIR filter
	}

	double pitch_old = pitch_a; // store old value
	double roll_old = roll_a; // store old value
	
	double sin_phi = (double)0.0; // tentative value
	double sin_teta = -a_f_cor_fil[0]; // X
	double cos_teta = sqrt((double)1.0 - sin_teta*sin_teta);
	if (cos_teta != (double)0.0) {
		sin_phi = a_f_cor_fil[1] / cos_teta; // Y
	}
	pitch_a = asin(sin_teta); // Pitch Acceleration (rad)
	roll_a = asin(sin_phi); // Roll Acceleration (rad)

	// Slope saturations (rationality limits) and Minimum and maximum saturations (rationality limits)
	double max_pitch_change_accepted = MAX_PITCH_VAR_ACC * integration_time; // max variation acceptable
	double max_roll_change_accepted = MAX_ROLL_VAR_ACC * integration_time; // max variation acceptable
	if ((abs(pitch_a - pitch_old) > max_pitch_change_accepted) || (abs(roll_a - roll_old) > max_roll_change_accepted) || (abs(pitch_a) > MAX_PITCH) || (abs(roll_a) > MAX_ROLL)) {
		if (pitch_a < (pitch_old - max_pitch_change_accepted)) pitch_a = pitch_old - max_pitch_change_accepted;
		if (pitch_a > (pitch_old + max_pitch_change_accepted)) pitch_a = pitch_old + max_pitch_change_accepted;
		if (roll_a < (roll_old - max_roll_change_accepted)) roll_a = roll_old - max_roll_change_accepted;
		if (roll_a > (roll_old + max_roll_change_accepted)) roll_a = roll_old + max_roll_change_accepted;
		if (pitch_a < -MAX_PITCH) pitch_a = -MAX_PITCH;
		if (pitch_a > MAX_PITCH) pitch_a = MAX_PITCH;
		if (roll_a < -MAX_ROLL) roll_a = -MAX_ROLL;
		if (roll_a > MAX_ROLL) roll_a = MAX_ROLL;
	}

	double a_f_delta[3]; // delta between present acceleration vector, and average vector
	a_f_delta_norm = 0.0;
	for (unsigned int i = 0; i < 3; i++) { // 3 axis acceleration delta (present sample - filtered sample)
		a_f_delta[i] = a_f_cor[i] - a_f_cor_fil[i];
		a_f_delta_norm += pow(a_f_delta[i], 2);
	}
	a_f_delta_norm = sqrt(a_f_delta_norm);

	// Filtering coefficient for fusion
	filter_req_pitchroll = FUSION_FILTER_MAX*exp(-a_f_delta_norm*(double)100.0);
	
	printf("%f %f \n", (float)a_f_delta_norm, (float)filter_req_pitchroll);
	//printf("%f %f \n", (float)pitch_a, (float)roll_a);
}

// This block fuses yaw, pitch and roll angles, depending on the reliability of the source
void IMU_Fuelino_class::fuse_yawpitchroll_block() {
	yaw_fuse = yaw_out; // use gyroscope previously integrated data
	pitch_fuse = filter_req_pitchroll*pitch_a + ((double)1.0 - filter_req_pitchroll)*pitch_out; // IIR filter
	roll_fuse = filter_req_pitchroll*roll_a + ((double)1.0 - filter_req_pitchroll)*roll_out; // IIR filter
}

// Evaluates if a complete turn has been done in clockwise / anticlockwise direction, by checking the difference between 2 consecutive angles
int IMU_Fuelino_class::evaluate_yaw_rotation_turns(double yaw_now, double* yaw_prev_ptr, int turns_prev) {
	int turns_new_tmp = turns_prev; // Starting value
	double angle_difference = abs(yaw_now - *yaw_prev_ptr);
	if (angle_difference > M_PI) { // difference higher than 180 deg
		if ((segno(yaw_now) == (double)-1.0) && (segno(*yaw_prev_ptr) == (double)1.0)) {
			turns_new_tmp = turns_prev + 1; // rotation was in clockwise direction
		}
		else if ((segno(yaw_now) == (double)1.0) && (segno(*yaw_prev_ptr) == (double)-1.0)) {
			turns_new_tmp = turns_prev - 1; // rotation was in anti-clockwise direction
		}
	}
	*yaw_prev_ptr = yaw_now; // stores the old value for next calculation
	return turns_new_tmp;
}

// Calculation task assumes that raw data (ax, ay, az, gx, gy, gz) has already been acquired
// This function is the main block. It performs all calculations, starting from data correction to yaw, pitch, roll calculation
void IMU_Fuelino_class::calculation_task() {
		
	process_raw_data(); // corrects sign and scale of raw data, corrects offsets, and transforms into floating point format
	calculate_pitchroll_from_acceleration(); // calculates pitch and roll angles from acceleration data

#if CAVALIERE_FILTER_ENABLE == 1
	// CAVALIERE FILTER: in case Cavaliere filter is used (monocilindro.com)		
	//statistics_calculation(); // updates filter and offsets values
	fuse_yawpitchroll_block(); // fuses sensor data depending on their reliability
	calculate_rotation_matrix_elements(yaw_fuse, pitch_fuse, roll_fuse, M_rot); // calculates rotation matrix elements
	calculate_quaternion_elements(M_rot, q_int_in); // calculates quaternion starting from rotation matrix
	integrate_quaternion_g(q_int_in, q_int_out); // calculates "gyroscope" quaternion by integrating the gyroscope data, based on input quaternion
#else
	// MADGWICK FILTER: in case Madgwick filter is used
	// Madgwick original algoritm works using an inverted quaternion, therefore the components have to be inverted before using the algorithm
	q_int_in[0] = q_int_out[0]; // sign inversion
	q_int_in[1] = -q_int_out[1];
	q_int_in[2] = -q_int_out[2];
	q_int_in[3] = -q_int_out[3];
	Madgwick_updateIMU(q_int_in, a_f_cor_fil[0], a_f_cor_fil[1], a_f_cor_fil[2], g_f_cor[0], g_f_cor[1], g_f_cor[2], integration_time, filter_req_pitchroll, q_int_out);
	q_int_out[1] = -q_int_out[1]; // sign inversion
	q_int_out[2] = -q_int_out[2];
	q_int_out[3] = -q_int_out[3];
#endif

	/*for (unsigned int i = 0; i < 4; i++) {
		q_f_out[i] = q_int_out[i];
	}*/	

	quaternion_product(q_M_cor, q_int_out, q_f_out); // corrects the rotation taking into account of sensor position rotation with respect to the main reference system
	calculate_yawpitchroll_from_quaternions(q_f_out, &yaw_out, &pitch_out, &roll_out); // calculates yaw pitch roll angles from quaternion
	turns_clockwise = evaluate_yaw_rotation_turns(yaw_out, &yaw_old, turns_clockwise); // calculate if a complete turn has been done (clockwise / anticlockwise) based on the yaw angle

}

// Initialization function
void IMU_Fuelino_class::IMU_init() {

	integration_time = 0.0; // Integration time [s]
	yaw_out = 0.0; // Output [rad]
	pitch_out = 0.0; // Output [rad]
	roll_out = 0.0; // Output [rad]
	turns_clockwise = 0; // how many turns clockwise have been done
	pitch_a = 0.0; // pitch calculated from accelerometer data
	roll_a = 0.0; // roll calculated from accelerometer data
	yaw_old = 0.0; // used to estimate rotation turns
	a_f_cor_fil[0] = 0.0; // X
	a_f_cor_fil[1] = 0.0; // Y
	a_f_cor_fil[2] = 1.0; // Z

	// Output quaternion starting values (0 degrees rotation)
	q_int_out[0] = 1.0;
	q_int_out[1] = 0.0;
	q_int_out[2] = 0.0;
	q_int_out[3] = 0.0;

	// Starting values of offsets for Gyroscope raw signal correction
	gx_f_offset = 0.0;
	gy_f_offset = 0.0;
	gz_f_offset = 0.0;

	// Resets acceleration sensors history
	a_samples_index = 0;
	for (unsigned int i = 0; i < A_HIST_NUM; i++) {
		a_mod_samples_hist[i] = 0.0;
	}

	// Resets gyroscope sensors history
	g_samples_index = 0;
	for (unsigned int i = 0; i < G_HIST_NUM; i++) {
		gx_samples_hist[i] = 0.0;
		gy_samples_hist[i] = 0.0;
		gz_samples_hist[i] = 0.0;
	}

	// Filtering
	filter_req_pitchroll = FUSION_FILTER_MAX; // filter constant for yaw and pitch (depends on how acceleration data is reliable)
	
	// Calculate correction matric elements and quaternion
	double s_teta = sin(PITCH_OFFSET);
	double c_teta = cos(PITCH_OFFSET);
	double s_phi = sin(ROLL_OFFSET);
	double c_phi = cos(ROLL_OFFSET);
	M_cor[0][0] = c_teta;
	M_cor[0][1] = s_phi * s_teta;
	M_cor[0][2] = c_phi * s_teta;
	M_cor[1][0] = 0.0;
	M_cor[1][1] = c_phi;
	M_cor[1][2] = -s_phi;
	M_cor[2][0] = -s_teta;
	M_cor[2][1] = s_phi * c_teta;
	M_cor[2][2] = c_phi * c_teta;
	calculate_quaternion_elements(M_cor, q_M_cor);
}

// Initializes the object
IMU_Fuelino_class::IMU_Fuelino_class() {
	IMU_init();
}