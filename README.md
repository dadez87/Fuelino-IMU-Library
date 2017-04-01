# Fuelino-IMU-Library
IMU Library, based on Madgwick's filter, for Fuelino fuel injection and motorcycle data logger.

Davide Cavaliere
www.monocilindro.com/Fuelino
dadez87@gmail.com
1 April 2017

In case you don't need it, please remove "#include "stdafx.h".

This library successfully works with raw signal of MPU6050.
To run the IMU, all you need to do is to set the raw signals of accelerometer and gyroscope. And the integration time (s).

IMU_Fuelino->ax_raw
IMU_Fuelino->ay_raw
IMU_Fuelino->az_raw
IMU_Fuelino->gx_raw
IMU_Fuelino->gy_raw
IMU_Fuelino->gz_raw
IMU_Fuelino->integration_time

Then launch the Calculation Task:

IMU_Fuelino->calculation_task();

The calculated Yaw, Pitch, Roll, are available as:

IMU_Fuelino->yaw_out
IMU_Fuelino->pitch_out
IMU_Fuelino->roll_out

By modifying the file "IMU_Fuelino_cal.h", you can calibrate:
1) Sensors offset and gain
2) Filters constants
3) Rationality limits

