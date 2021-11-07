// Dem.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include"..\include\Dem.h"
using namespace cv;
inline bool return_check(int ret, const char* detail_info, const char* error_head)
{
	if (ret < 0)
	{
		fprintf(stderr, "%s %s\n\n", error_head, detail_info);
		return true;
	}
	else
	{
		return false;
	}
}
inline bool parallel_check(volatile bool parallel_flag, const char* detail_info, const char* parallel_error_head)
{
	if (!parallel_flag)
	{
		fprintf(stderr, "%s %s\n\n", parallel_error_head, detail_info);
		return true;
	}
	else
	{
		return false;
	}
}
int Dem::phase2dem_newton_iter(
	Mat unwrapped_phase,
	Mat flat_phase,
	Mat& DEM_height,
	Mat auxi_m,
	Mat auxi_s,
	Mat orbit_m,
	Mat orbit_s,
	Mat doppler_frequency,
	Mat gcps,
	Mat regis_out,
	double delta_m,
	double delta_s,
	int multilook_times,
	int mode,
	int iters = 30
)
{
	if (unwrapped_phase.rows < 1 ||
		unwrapped_phase.cols < 1 ||
		unwrapped_phase.type() != CV_64F ||
		unwrapped_phase.channels() != 1 ||
		flat_phase.rows != unwrapped_phase.rows ||
		flat_phase.cols != unwrapped_phase.cols ||
		flat_phase.type() != CV_64F ||
		flat_phase.channels() != 1 ||
		auxi_m.rows < 1 ||
		auxi_m.cols != 5 ||
		auxi_m.type() != CV_64F ||
		auxi_m.channels() != 1 ||
		auxi_s.rows != auxi_m.rows ||
		auxi_s.cols != auxi_m.cols ||
		auxi_s.type() != CV_64F ||
		auxi_s.channels() != 1 ||
		orbit_m.rows < 7 ||
		orbit_m.cols != 7 ||
		orbit_m.type() != CV_64F ||
		orbit_m.channels() != 1 ||
		orbit_s.rows < 7 ||
		orbit_s.cols != 7 ||
		orbit_s.type() != CV_64F ||
		orbit_s.channels() != 1 ||
		doppler_frequency.rows != 1 ||
		doppler_frequency.cols < 1 ||
		doppler_frequency.type() != CV_64F ||
		doppler_frequency.channels() != 1 ||
		gcps.rows < 1 ||
		gcps.cols != 5 ||
		gcps.type() != CV_64F ||
		gcps.channels() != 1 ||
		regis_out.rows != 1 ||
		regis_out.cols != 4 ||
		regis_out.type() != CV_32S||
		regis_out.channels() != 1||
		delta_m <= 0.0 ||
		delta_s <= 0.0 ||
		multilook_times < 1||
		mode < 1||
		mode > 2
		)
	{
		fprintf(stderr, "phase2dem_newton_iter(): input check failed!\n\n");
		return -1;
	}
	int move_r = regis_out.at<int>(0, 0);
	int move_c = regis_out.at<int>(0, 1);
	int nr = regis_out.at<int>(0, 2);
	int nc = regis_out.at<int>(0, 3);
	if (nr < 1 || nc < 1 || fabs(move_c) >= nc || fabs(move_r) >= nr)
	{
		fprintf(stderr, "phase2dem_newton_iter(): input check failed!\n\n");
		return -1;
	}
	double C = 4 * 3.1415926535;
	if (mode == 2)
	{
		C = 2 * 3.1415926535;
	}
	Deflat deflat;
	Mat coef_m, coef_s;
	int ret;
	ret = deflat.Orbit_Polyfit(orbit_m, coef_m);
	if (return_check(ret, "deflat.Orbit_Polyfit(*, *)", error_head)) return -1;
	ret = deflat.Orbit_Polyfit(orbit_s, coef_s);
	if (return_check(ret, "deflat.Orbit_Polyfit(*, *)", error_head)) return -1;
	Mat image_time_m = Mat::zeros(1, nr, CV_64F);
	Mat image_time_s = Mat::zeros(1, nr, CV_64F);
	Mat temp;
	Mat S_position_m, S_position_s, S_velocity_m;
	S_position_m = Mat::zeros(nr, 3, CV_64F);
	S_position_s = Mat::zeros(nr, 3, CV_64F);
	S_velocity_m = Mat::zeros(nr, 3, CV_64F);
	for (int i = 0; i < nr; i++)
	{
		image_time_m.at<double>(0, i) = double(i) * delta_m + auxi_m.at<double>(0, 2);
		image_time_s.at<double>(0, i) = double(i) * delta_s + auxi_s.at<double>(0, 2);
	}
	for (int i = 0; i < nr; i++)
	{
		ret = deflat.get_xyz(image_time_m.at<double>(0, i), coef_m, temp);
		if (return_check(ret, "deflat.get_xyz(*, *, *)", error_head)) return -1;
		S_position_m.at<double>(i, 0) = temp.at<double>(0, 0);
		S_position_m.at<double>(i, 1) = temp.at<double>(0, 1);
		S_position_m.at<double>(i, 2) = temp.at<double>(0, 2);
		ret = deflat.get_xyz(image_time_s.at<double>(0, i), coef_s, temp);
		if (return_check(ret, "deflat.get_xyz(*, *, *)", error_head)) return -1;
		S_position_s.at<double>(i, 0) = temp.at<double>(0, 0);
		S_position_s.at<double>(i, 1) = temp.at<double>(0, 1);
		S_position_s.at<double>(i, 2) = temp.at<double>(0, 2);
		ret = deflat.get_vel(image_time_m.at<double>(0, i), coef_m, temp);
		if (return_check(ret, "deflat.get_vel(*, *, *)", error_head)) return -1;
		S_velocity_m.at<double>(i, 0) = temp.at<double>(0, 0);
		S_velocity_m.at<double>(i, 1) = temp.at<double>(0, 1);
		S_velocity_m.at<double>(i, 2) = temp.at<double>(0, 2);
	}
	int rows_start_m, rows_end_m, rows_start_s, rows_end_s, cols_start_m, cols_start_s, cols_end_m, cols_end_s;
	if (move_r >= 0)
	{
		rows_start_m = 0;
		rows_end_m = nr - move_r;
		rows_start_s = move_r;
		rows_end_s = nr;
	}
	else
	{
		rows_start_m = -move_r;
		rows_end_m = nr;
		rows_start_s = 0;
		rows_end_s = nr + move_r;
	}
	if (move_c >= 0)
	{
		cols_start_m = 0;
		cols_end_m = nc - move_c;
		cols_start_s = move_c;
		cols_end_s = nc;
	}
	else
	{
		cols_start_m = -move_c;
		cols_end_m = nc;
		cols_start_s = 0;
		cols_end_s = nc + move_c;
	}
	if ((rows_end_m - rows_start_m) != unwrapped_phase.rows ||
		(rows_end_s - rows_start_s) != unwrapped_phase.rows ||
		(cols_end_m - cols_start_m) != unwrapped_phase.cols ||
		(cols_end_s - cols_start_s) != unwrapped_phase.cols)
	{
		fprintf(stderr, "offset cutsize and unwrapped_phase size mismatch!\n\n");
		return -1;
	}
	S_position_m = S_position_m(Range(rows_start_m, rows_end_m), Range(0, 3));
	S_position_s = S_position_s(Range(rows_start_s, rows_end_s), Range(0, 3));
	S_velocity_m = S_velocity_m(Range(rows_start_m, rows_end_m), Range(0, 3));
	doppler_frequency = doppler_frequency(Range(0, 1), Range(cols_start_m, cols_end_m));
	unwrapped_phase = unwrapped_phase + flat_phase;
	
	int row = gcps.at<double>(2, 0);
	int col = gcps.at<double>(2, 1);
	if (row < 1 ||
		row > S_position_m.rows ||
		row < 1 ||
		row > S_position_s.rows)
	{
		fprintf(stderr, "ground control point out of image !\n\n");
		return -1;
	}

	///////////////////////////////校正至绝对相位///////////////////////////////////
	Mat Control_Point_Position = gcps(Range(2, 3), Range(2, 5));
	Mat tmp = Control_Point_Position - 
		S_position_m(Range(row - 1, row), Range(0, 3));
	double distance_r_m = 2 * cv::norm(tmp, NORM_L2);
	tmp = Control_Point_Position - 
		S_position_m(Range(row - 1, row), Range(0, 3));
	Mat tmp1 = Control_Point_Position - 
		S_position_s(Range(row - 1, row), Range(0, 3));
	double distance_r_s;
	if (mode == 1)
	{
		distance_r_s = 2 * cv::norm(tmp1, NORM_L2);
	}
	else
	{
		distance_r_s = cv::norm(tmp, NORM_L2) + cv::norm(tmp1, NORM_L2);
	}
	if (fabs(auxi_m.at<double>(0, 4)) < 1e-14)
	{
		auxi_m.at<double>(0, 4) = 1e-14;
	}
	double lambda = 3e8 / auxi_m.at<double>(0, 4);
	double real_phase = C / lambda * (distance_r_m - distance_r_s);
	double K = round((real_phase - unwrapped_phase.at<double>(row - 1, col - 1)) / C);
	unwrapped_phase = unwrapped_phase + K * C;


	double R0_m = auxi_m.at<double>(0, 3) * 3e8 / 2;//最近斜距
	double spc_m = auxi_m.at<double>(0, 1);//距离向采样间隔
	Mat Rm = Mat::zeros(1, nc, CV_64F);
	for (int i = 0; i < nc; i++)
	{
		Rm.at<double>(0, i) = double(i) * spc_m + R0_m;
	}
	Rm = Rm(Range(0, 1), Range(cols_start_m, cols_end_m));
	Mat ones = Mat::ones(unwrapped_phase.rows, 1, CV_64F);
	Mat R_M = ones * Rm;
	Mat R_F = ones * Rm * 2.0 - lambda * unwrapped_phase / C;
	Mat Satellite_M_T_Position = S_position_m;//主星发射位置
	Mat Satellite_S_T_Position;
	if (mode == 1)
	{
		Satellite_S_T_Position = S_position_s;//辅星发射位置
	}
	else
	{
		Satellite_S_T_Position = S_position_m;//辅星发射位置
	}
	Mat Satellite_M_R_Position = S_position_m;//主星接收位置
	Mat Satellite_S_R_Position = S_position_s;//辅星接收位置
	Mat Satellite_M = (Satellite_M_R_Position + Satellite_M_T_Position) / 2;
	Mat Vs = S_velocity_m;
	ones = Mat::ones(unwrapped_phase.rows, unwrapped_phase.cols, CV_64F);
	Mat P1 = ones * Control_Point_Position.at<double>(0, 0);
	Mat P2 = ones * Control_Point_Position.at<double>(0, 1);
	Mat P3 = ones * Control_Point_Position.at<double>(0, 2);
	int fine_size_rows = unwrapped_phase.rows;
	int fine_size_cols = unwrapped_phase.cols;
	Mat M_T, M_R, S_T, S_R;
	Mat f1, f2, f3;
	Mat det_Df, Df_ni11, Df_ni12, Df_ni13, Df_ni21, Df_ni22, Df_ni23;
	Mat Df11, Df12, Df13, Df21, Df22, Df23, Df31, Df32, Df33, Df_ni31, Df_ni32, Df_ni33;
	Mat delta_Rt1, delta_Rt2, delta_Rt3;
	ones = Mat::ones(1, fine_size_cols, CV_64F);
	Mat temp_var, temp_var1;
	Mat fd = doppler_frequency;
	for (int i = 0; i < iters; i++)
	{
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(M_T);
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		M_T = M_T + temp_var;
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		M_T = M_T + temp_var;
		cv::sqrt(M_T, f1);
		f1 = f1 * 2 - 2 * R_M;

		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(S_T);
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		S_T = S_T + temp_var;
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		S_T = S_T + temp_var;

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(S_R);
		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		S_R = S_R + temp_var;
		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		S_R = S_R + temp_var;

		cv::sqrt(S_T, f2);
		cv::sqrt(S_R, temp_var);
		f2 = f2 + temp_var - R_F;


		temp_var = Vs(Range(0, Vs.rows), Range(0, 1)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(0, 1)) * ones - P1;
		f3 = temp_var.mul(temp_var1);
		temp_var = Vs(Range(0, Vs.rows), Range(1, 2)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(1, 2)) * ones - P2;
		f3 = f3 + temp_var.mul(temp_var1);
		temp_var = Vs(Range(0, Vs.rows), Range(2, 3)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(2, 3)) * ones - P3;
		f3 = f3 + temp_var.mul(temp_var1);
		ones = Mat::ones(fine_size_rows, 1, CV_64F);
		temp_var = ones * fd;
		temp_var1 = R_M * lambda / 2.0;
		f3 = f3 + temp_var.mul(temp_var1);

		//Dff
	    //第一行：f(1)的x，y，z的导数
		ones = Mat::ones(1, fine_size_cols, CV_64F);
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df11 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df11 = Df11 + temp_var.mul(temp_var1);

		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df12 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df12 = Df12 + temp_var.mul(temp_var1);

		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df13 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df13 = Df13 + temp_var.mul(temp_var1);


		//第二行：f(2)的x，y，z的导数
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df21 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df21 = Df21 + temp_var.mul(temp_var1);


		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df22 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df22 = Df22 + temp_var.mul(temp_var1);


		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df23 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df23 = Df23 + temp_var.mul(temp_var1);

		//第三行：f(3)的x，y，z的导数
		Df31 = -Vs(Range(0, Vs.rows), Range(0, 1)) * ones;
		Df32 = -Vs(Range(0, Vs.rows), Range(1, 2)) * ones;
		Df33 = -Vs(Range(0, Vs.rows), Range(2, 3)) * ones;

		temp_var = Df11.mul(Df22);
		temp_var = temp_var.mul(Df33);
		temp_var.copyTo(det_Df);

		temp_var = Df12.mul(Df23);
		temp_var = temp_var.mul(Df31);
		det_Df = det_Df + temp_var;

		temp_var = Df13.mul(Df21);
		temp_var = temp_var.mul(Df32);
		det_Df = det_Df + temp_var;

		temp_var = Df31.mul(Df22);
		temp_var = temp_var.mul(Df13);
		det_Df = det_Df - temp_var;

		temp_var = Df32.mul(Df23);
		temp_var = temp_var.mul(Df11);
		det_Df = det_Df - temp_var;

		temp_var = Df33.mul(Df21);
		temp_var = temp_var.mul(Df12);
		det_Df = det_Df - temp_var;

		Df_ni11 = (Df22.mul(Df33) - Df32.mul(Df23)) / det_Df;
		Df_ni12 = -(Df12.mul(Df33) - Df32.mul(Df13)) / det_Df;
		Df_ni13 = (Df12.mul(Df23) - Df22.mul(Df13)) / det_Df;
		delta_Rt1 = Df_ni11.mul(f1) + Df_ni12.mul(f2) + Df_ni13.mul(f3);


		Df_ni21 = -(Df21.mul(Df33) - Df31.mul(Df23)) / det_Df;
		Df_ni22 = (Df11.mul(Df33) - Df31.mul(Df13)) / det_Df;
		Df_ni23 = -(Df11.mul(Df23) - Df21.mul(Df13)) / det_Df;
		delta_Rt2 = Df_ni21.mul(f1) + Df_ni22.mul(f2) + Df_ni23.mul(f3);

		Df_ni31 = (Df21.mul(Df32) - Df31.mul(Df22)) / det_Df;
		Df_ni32 = -(Df11.mul(Df32) - Df31.mul(Df12)) / det_Df;
		Df_ni33 = (Df22.mul(Df11) - Df21.mul(Df12)) / det_Df;
		delta_Rt3 = Df_ni31.mul(f1) + Df_ni32.mul(f2) + Df_ni33.mul(f3);

		P1 = P1 - delta_Rt1;
		P2 = P2 - delta_Rt2;
		P3 = P3 - delta_Rt3;
	}
	volatile bool parallel_flag = true;
	DEM_height = Mat::zeros(fine_size_rows, fine_size_cols, CV_64F);
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (int i = 0; i < fine_size_rows; i++)
	{
		if (!parallel_flag) continue;
		Utils util;
		for (int j = 0; j < fine_size_cols; j++)
		{
			if (!parallel_flag) continue;
			Mat xyz, llh;
			xyz = Mat::zeros(1, 3, CV_64F);
			xyz.at<double>(0, 0) = P1.at<double>(i, j);
			xyz.at<double>(0, 1) = P2.at<double>(i, j);
			xyz.at<double>(0, 2) = P3.at<double>(i, j);
			
			ret = util.xyz2ell(xyz, llh);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			DEM_height.at<double>(i, j) = llh.at<double>(0, 2);
		}
	}
	if (parallel_check(parallel_flag, "phase2dem_newton_iter()", parallel_error_head)) return -1;
	Mat xyz, llh;
	Utils util;
	ret = util.xyz2ell(Control_Point_Position, llh);
	if (return_check(ret, "util.xyz2ell(*, *)", error_head)) return -1;
	DEM_height = DEM_height + llh.at<double>(0, 2) - DEM_height.at<double>(row - 1, col - 1);
	return 0;
}

