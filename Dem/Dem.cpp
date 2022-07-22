// Dem.cpp : 定义 DLL 应用程序的导出函数。
//
#include "stdafx.h"
#include"..\include\Dem.h"
#include"..\include\tinyxml.h"
#include"..\include\FormatConversion.h"
#include"..\include\ComplexMat.h"
#ifdef _DEBUG
#pragma comment(lib, "FormatConversion_d.lib")
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "Deflat_d.lib")
#pragma comment(lib, "ComplexMat_d.lib")
#else
#pragma comment(lib, "FormatConversion.lib")
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "Deflat.lib")
#pragma comment(lib, "ComplexMat.lib")
#endif // _DEBUG
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

Dem::Dem()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "DEM_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "DEM_DLL_ERROR: error happens when using parallel computing in function: ");
}

Dem::~Dem()
{
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

int Dem::dem_newton_iter(const char* unwrapped_phase_file, Mat& dem, const char* project_path, int iter_times, int mode)
{
	if (unwrapped_phase_file == NULL ||
		project_path == NULL ||
		iter_times < 1 ||
		mode < 1 || mode > 2)
	{
		fprintf(stderr, "dem_newton_iter(): input check failed!\n");
		return -1;
	}

	/*
	* 校正至绝对相位
	*/

	FormatConversion conversion; Utils util;
	int nr, nc, ret, offset_row, offset_col;
	double time_interval1, time_interval2;
	string source_1, source_2, tmp;
	string project(project_path);
	Mat unwrapped_phase, flat_phase_coefficient, gcps, temp, range_spacing,
		stateVec1, stateVec2, lat_coefficient, lon_coefficient, prf1, prf2, carrier_frequency;
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "phase", unwrapped_phase);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nr = unwrapped_phase.rows; nc = unwrapped_phase.cols;
	if (nr < 1 || nc < 1)
	{
		fprintf(stderr, "dem_newton_iter(): invalide unwrapped_phase !\n");
		return -1;
	}
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "flat_phase_coefficient", flat_phase_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_1", source_1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_2", source_2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	tmp = project;
	source_1 = tmp + source_1;
	source_2 = tmp + source_2;
	ret = conversion.read_array_from_h5(source_1.c_str(), "gcps", gcps);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_row", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_row = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_col", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_col = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "range_spacing", range_spacing);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "state_vec", stateVec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_2.c_str(), "state_vec", stateVec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "prf", prf1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval1 = 1.0 / (prf1.at<double>(0, 0) + 1e-10);
	ret = conversion.read_array_from_h5(source_2.c_str(), "prf", prf2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval2 = 1.0 / (prf2.at<double>(0, 0) + 1e-10);
	ret = conversion.read_array_from_h5(source_1.c_str(), "lat_coefficient", lat_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "lon_coefficient", lon_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "carrier_frequency", carrier_frequency);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	//寻找图像范围内的控制点信息

	int num_gcps = gcps.rows;
	int row, col, i = 0;
	bool b_gcp = false;
	for (i = 0; i < num_gcps; i++)
	{
		row = (int)gcps.at<double>(i, 3);
		col = (int)gcps.at<double>(i, 4);
		if ((row - offset_row) >= 0 && (row - offset_row) < nr && (col - offset_col) >= 0 && (col - offset_col) < nc)
		{
			b_gcp = true;
			break;
		}
	}
	
	Mat llh(1, 3, CV_64F), xyz_ground(1, 3, CV_64F);
	Mat row_coord(1, 1, CV_64F), col_coord(1, 1, CV_64F), lat, lon;
	if (b_gcp)
	{
		llh.at<double>(0, 0) = gcps.at<double>(i, 1);
		llh.at<double>(0, 1) = gcps.at<double>(i, 0);
		llh.at<double>(0, 2) = gcps.at<double>(i, 2);
		row = row - offset_row; col = col - offset_col;
	}
	else
	{
		row_coord.at<double>(0, 0) = offset_row + double(nr) / 2.0;
		col_coord.at<double>(0, 0) = offset_col + double(nc) / 2.0;
		ret = util.coord_conversion(lat_coefficient, row_coord, col_coord, lat);
		if (return_check(ret, "coord_conversion", error_head)) return -1;
		ret = util.coord_conversion(lon_coefficient, row_coord, col_coord, lon);
		if (return_check(ret, "coord_conversion", error_head)) return -1;
		llh.at<double>(0, 0) = lat.at<double>(0, 0);
		llh.at<double>(0, 1) = lon.at<double>(0, 0);
		llh.at<double>(0, 2) = 0.0;
		row = (int)nr / 2; col = (int)nc / 2;
	}
	ret = util.ell2xyz(llh, xyz_ground);
	if (return_check(ret, "ell2xyz", error_head)) return -1;


	/*
	* 轨道插值
	*/
	ret = util.stateVec_interp(stateVec1, time_interval1, stateVec1);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	ret = util.stateVec_interp(stateVec2, time_interval2, stateVec2);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	/*
	* 寻找图像左上角成像卫星位置
	*/
	Mat sate1_xyz, sate2_xyz, sate1_v, sate2_v;
	Mat sate1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat satev1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat sate2 = Mat::zeros(nr, 3, CV_64F);//辅星位置
	stateVec1(cv::Range(0, stateVec1.rows), cv::Range(1, 4)).copyTo(sate1_xyz);
	stateVec1(cv::Range(0, stateVec1.rows), cv::Range(4, 7)).copyTo(sate1_v);
	stateVec2(cv::Range(0, stateVec2.rows), cv::Range(1, 4)).copyTo(sate2_xyz);
	stateVec2(cv::Range(0, stateVec2.rows), cv::Range(4, 7)).copyTo(sate2_v);
	//找到零多普勒位置
	Mat xyz, llh_upperleft(1, 3, CV_64F);
	row_coord.at<double>(0, 0) = offset_row;
	col_coord.at<double>(0, 0) = offset_col + double(nc) / 2.0;
	ret = util.coord_conversion(lat_coefficient, row_coord, col_coord, lat);
	if (return_check(ret, "coord_conversion", error_head)) return -1;
	ret = util.coord_conversion(lon_coefficient, row_coord, col_coord, lon);
	if (return_check(ret, "coord_conversion", error_head)) return -1;
	llh_upperleft.at<double>(0, 0) = lat.at<double>(0, 0);
	llh_upperleft.at<double>(0, 1) = lon.at<double>(0, 0);
	llh_upperleft.at<double>(0, 2) = 0.0;
	ret = util.ell2xyz(llh_upperleft, xyz);
	if (return_check(ret, "ell2xyz", error_head)) return -1;
	Mat dop = Mat::zeros(sate1_xyz.rows, 1, CV_64F);
	Mat r;
	for (int j = 0; j < sate1_xyz.rows; j++)
	{
		r = xyz - sate1_xyz(Range(j, j + 1), Range(0, 3));
		dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate1_v(Range(j, j + 1), Range(0, 3))))[0]);
	}
	Point peak_loc;
	cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
	int xxxx;
	for (int j = 0; j < nr; j++)
	{
		xxxx = (peak_loc.y + j) > (sate1_xyz.rows - 1) ? (sate1_xyz.rows - 1) : (peak_loc.y + j);
		sate1_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate1(Range(j, j + 1), Range(0, 3)));
		sate1_v(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(satev1(Range(j, j + 1), Range(0, 3)));
	}
	//卫星2
	dop = Mat::zeros(sate2_xyz.rows, 1, CV_64F);
	for (int j = 0; j < sate2_xyz.rows; j++)
	{
		r = xyz - sate2_xyz(Range(j, j + 1), Range(0, 3));
		dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate2_v(Range(j, j + 1), Range(0, 3))))[0]);
	}
	cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
	for (int j = 0; j < nr; j++)
	{
		xxxx = (peak_loc.y + j) > (sate2_xyz.rows - 1) ? (sate2_xyz.rows - 1) : (peak_loc.y + j);
		sate2_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate2(Range(j, j + 1), Range(0, 3)));
	}
	//加回平地相位
	//Mat coef;
	//cv::transpose(flat_phase_coefficient, coef);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		Mat temp(1, 6, CV_64F);
		for (int j = 0; j < nc; j++)
		{
			temp.at<double>(0, 0) = 1.0;
			temp.at<double>(0, 1) = i;
			temp.at<double>(0, 2) = j;
			temp.at<double>(0, 3) = i * j;
			temp.at<double>(0, 4) = i * i;
			temp.at<double>(0, 5) = j * j;
			unwrapped_phase.at<double>(i, j) = unwrapped_phase.at<double>(i, j) + sum(temp.mul(flat_phase_coefficient))[0];
		}
	}
	//控制点绝对相位计算
	double r_main = sqrt(sum((sate1(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground).mul(sate1(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground))[0]);
	double r_slave = sqrt(sum((sate2(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground).mul(sate2(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground))[0]);
	double C = mode == 1 ? 4 * PI : 2 * PI;
	double lambda = 3e8 / (carrier_frequency.at<double>(0, 0) + 1e-10);
	double phase_real = (r_slave - r_main) / lambda * C;
	double K = round((phase_real - unwrapped_phase.at<double>(row, col)) / (2 * PI));
	unwrapped_phase = unwrapped_phase + K * 2 * PI;//相位校正
	
	//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\unwrapped_phase_abs.bin", unwrapped_phase);

	/*
	* 反演高程
	*/

	Mat R_M(1, nc, CV_64F);
	for (int i = 0; i < nc; i++)
	{
		R_M.at<double>(0, i) = r_main + range_spacing.at<double>(0, 0) * (double(i) - col);
	}
	Mat ones = Mat::ones(nr, 1, CV_64F);
	R_M = ones * R_M;
	Mat R_F = R_M * 2.0 + lambda * unwrapped_phase / (2 * PI);
	Mat Satellite_M_T_Position = sate1;//主星发射位置
	Mat Satellite_S_T_Position;
	if (mode == 1)
	{
		Satellite_S_T_Position = sate2;//辅星发射位置
	}
	else
	{
		Satellite_S_T_Position = sate1;//辅星发射位置
	}
	Mat Satellite_M_R_Position = sate1;//主星接收位置
	Mat Satellite_S_R_Position = sate2;//辅星接收位置
	Mat Satellite_M = (Satellite_M_R_Position + Satellite_M_T_Position) / 2;
	Mat Vs = satev1;
	ones = Mat::ones(nr, nc, CV_64F);
	Mat P1 = ones * xyz_ground.at<double>(0, 0);
	Mat P2 = ones * xyz_ground.at<double>(0, 1);
	Mat P3 = ones * xyz_ground.at<double>(0, 2);
	Mat M_T, M_R, S_T, S_R;
	Mat f1, f2, f3;
	Mat det_Df, Df_ni11, Df_ni12, Df_ni13, Df_ni21, Df_ni22, Df_ni23;
	Mat Df11, Df12, Df13, Df21, Df22, Df23, Df31, Df32, Df33, Df_ni31, Df_ni32, Df_ni33;
	Mat delta_Rt1, delta_Rt2, delta_Rt3;
	ones = Mat::ones(1, nc, CV_64F);
	Mat temp_var, temp_var1;
	Mat fd = Mat::zeros(1, nc, CV_64F);
	for (int i = 0; i < iter_times; i++)
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
		ones = Mat::ones(nr, 1, CV_64F);
		temp_var = ones * fd;
		temp_var1 = R_M * lambda / 2.0;
		f3 = f3 + temp_var.mul(temp_var1);

		//Dff
		//第一行：f(1)的x，y，z的导数
		ones = Mat::ones(1, nc, CV_64F);
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
	delta_Rt1.release(); delta_Rt2.release(); delta_Rt3.release(); Df_ni31.release(); Df_ni32.release();
	Df_ni33.release(); Df_ni21.release(); Df_ni22.release(); Df_ni23.release(); Df_ni11.release(); Df_ni12.release();
	Df_ni13.release();
	volatile bool parallel_flag = true;
	dem.create(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (int i = 0; i < nr; i++)
	{
		if (!parallel_flag) continue;
		Utils util;
		for (int j = 0; j < nc; j++)
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
			dem.at<double>(i, j) = llh.at<double>(0, 2);
		}
	}
	if (parallel_check(parallel_flag, "dem_newton_iter()", parallel_error_head)) return -1;
	dem = dem + llh.at<double>(0, 2) - dem.at<double>(row, col);
	return 0;
}

int Dem::dem_newton_iter_test(const char* unwrapped_phase_file, Mat& dem, const char* project_path, int iter_times, int mode)
{
	if (unwrapped_phase_file == NULL ||
		project_path == NULL ||
		iter_times < 1 ||
		mode < 1 || mode > 2)
	{
		fprintf(stderr, "dem_newton_iter(): input check failed!\n");
		return -1;
	}

	/*
	* 校正至绝对相位
	*/

	FormatConversion conversion; Utils util;
	int nr, nc, ret, offset_row, offset_col;
	double time_interval1, time_interval2, acquisitionStartTime1, acquisitionStartTime2, acquisitionStopTime1,
		acquisitionStopTime2, wavelength, nearRange;
	string source_1, source_2, tmp, start_time, end_time;
	string project(project_path);
	Mat unwrapped_phase, flat_phase_coefficient, gcps, temp, range_spacing,
		stateVec1, stateVec2, lat_coefficient, lon_coefficient, prf1, prf2, carrier_frequency;
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "phase", unwrapped_phase);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nr = unwrapped_phase.rows; nc = unwrapped_phase.cols;
	if (nr < 1 || nc < 1)
	{
		fprintf(stderr, "dem_newton_iter(): invalide unwrapped_phase !\n");
		return -1;
	}
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "flat_phase_coefficient", flat_phase_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_1", source_1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_2", source_2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	tmp = project;
	source_1 = tmp + source_1;
	source_2 = tmp + source_2;
	ret = conversion.read_array_from_h5(source_1.c_str(), "GCP", gcps);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_row", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_row = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_col", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_col = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "range_spacing", range_spacing);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "state_vec", stateVec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_2.c_str(), "state_vec", stateVec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "prf", prf1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval1 = 1.0 / (prf1.at<double>(0, 0) + 1e-10);
	ret = conversion.read_array_from_h5(source_2.c_str(), "prf", prf2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval2 = 1.0 / (prf2.at<double>(0, 0) + 1e-10);
	//ret = conversion.read_array_from_h5(source_1.c_str(), "lat_coefficient", lat_coefficient);
	//if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	//ret = conversion.read_array_from_h5(source_1.c_str(), "lon_coefficient", lon_coefficient);
	//if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "carrier_frequency", carrier_frequency);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	wavelength = 3e8 / carrier_frequency.at<double>(0, 0);

	ret = conversion.read_str_from_h5(source_1.c_str(), "acquisition_start_time", start_time);
	ret = conversion.utc2gps(start_time.c_str(), &acquisitionStartTime1);
	ret = conversion.read_str_from_h5(source_1.c_str(), "acquisition_stop_time", end_time);
	ret = conversion.utc2gps(end_time.c_str(), &acquisitionStopTime1);
	ret = conversion.read_str_from_h5(source_2.c_str(), "acquisition_start_time", start_time);
	ret = conversion.utc2gps(start_time.c_str(), &acquisitionStartTime2);
	ret = conversion.read_str_from_h5(source_2.c_str(), "acquisition_stop_time", end_time);
	ret = conversion.utc2gps(end_time.c_str(), &acquisitionStopTime2);

	ret = conversion.read_double_from_h5(source_1.c_str(), "slant_range_first_pixel", &nearRange);

	//寻找图像范围内的控制点信息

	int num_gcps = gcps.rows;
	int row, col, i_gcp = 0, count = 0;
	bool b_gcp = false;
	vector<int> valid_row;
	for (i_gcp = 0; i_gcp < num_gcps; i_gcp++)
	{
		row = (int)gcps.at<double>(i_gcp, 0);
		col = (int)gcps.at<double>(i_gcp, 1);
		if ((row - offset_row - 1) >= 0 && (row - offset_row - 1) < nr && (col - offset_col - 1) >= 0 && (col - offset_col - 1) < nc)
		{
			b_gcp = true;
			valid_row.push_back(i_gcp);
			//break;
		}
	}
	Mat llh(1, 3, CV_64F), xyz_ground(1, 3, CV_64F);
	if (b_gcp)
	{
		i_gcp = valid_row[0];
		row = (int)gcps.at<double>(i_gcp, 0);
		col = (int)gcps.at<double>(i_gcp, 1);
		llh.at<double>(0, 0) = gcps.at<double>(i_gcp, 3);
		llh.at<double>(0, 1) = gcps.at<double>(i_gcp, 2);
		llh.at<double>(0, 2) = gcps.at<double>(i_gcp, 4);
		row = row - offset_row; col = col - offset_col;
	}
	else
	{
		return -1;
	}
	ret = util.ell2xyz(llh, xyz_ground);
	if (return_check(ret, "ell2xyz", error_head)) return -1;


	/*
	* 轨道插值
	*/
	orbitStateVectors stateVectors1(stateVec1, acquisitionStartTime1, acquisitionStopTime1);
	stateVectors1.applyOrbit();
	orbitStateVectors stateVectors2(stateVec2, acquisitionStartTime2, acquisitionStopTime2);
	stateVectors2.applyOrbit();
	/*
	* 寻找图像左上角成像卫星位置
	*/
	Mat sate1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat satev1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat sate2 = Mat::zeros(nr, 3, CV_64F);//辅星位置
	//找到零多普勒位置
	Position pos;Velocity vel;
	for (int i = 0; i < nr; i++)
	{
		stateVectors1.getPosition(acquisitionStartTime1 + double(offset_row + i) * time_interval1, pos);
		stateVectors1.getVelocity(acquisitionStartTime1 + double(offset_row + i) * time_interval1, vel);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
		satev1.at<double>(i, 0) = vel.vx;
		satev1.at<double>(i, 1) = vel.vy;
		satev1.at<double>(i, 2) = vel.vz;
	}

	//卫星2
	Position groundPosition;
	double latitude, longitude, height, dopplerFrequency = 0.0;
	latitude = gcps.at<double>(i_gcp, 3);
	longitude = gcps.at<double>(i_gcp, 2);
	longitude = longitude > 180.0 ? (longitude - 360.0) : longitude;
	height = gcps.at<double>(i_gcp, 4);
	Utils::ell2xyz(longitude, latitude, height, groundPosition);
	int numOrbitVec = stateVectors2.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int iii = 0; iii < numOrbitVec; iii++) {
		Position orb_pos(stateVectors2.newStateVectors.at<double>(iii, 1), stateVectors2.newStateVectors.at<double>(iii, 2),
			stateVectors2.newStateVectors.at<double>(iii, 3));
		Velocity orb_vel(stateVectors2.newStateVectors.at<double>(iii, 4), stateVectors2.newStateVectors.at<double>(iii, 5),
			stateVectors2.newStateVectors.at<double>(iii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = time_interval2;

	int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	int numIterations = 0; 
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors2.getPosition(midTime, pos);
		stateVectors2.getVelocity(midTime, vel);
		xdiff = groundPosition.x - pos.x;
		ydiff = groundPosition.y - pos.y;
		zdiff = groundPosition.z - pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength * distance);
		if ((midFreq - dopplerFrequency) * (lowerBoundFreq - dopplerFrequency) > 0.0) {
			lowerBoundTime = midTime;
			lowerBoundFreq = midFreq;
		}
		else if ((midFreq - dopplerFrequency) * (upperBoundFreq - dopplerFrequency) > 0.0) {
			upperBoundTime = midTime;
			upperBoundFreq = midFreq;
		}
		else if (fabs(midFreq - dopplerFrequency) < 0.01) {
			zeroDopplerTime = midTime;
			break;
		}

		diffTime = fabs(upperBoundTime - lowerBoundTime);
		numIterations++;
	}
	zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
	acquisitionStartTime2 = zeroDopplerTime - (row - 1) * time_interval2;
	for (int i = 0; i < nr; i++)
	{
		stateVectors2.getPosition(acquisitionStartTime2 + double(i) * time_interval2, pos);
		sate2.at<double>(i, 0) = pos.x;
		sate2.at<double>(i, 1) = pos.y;
		sate2.at<double>(i, 2) = pos.z;
	}

	//加回平地相位
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		Mat temp(1, 6, CV_64F);
		for (int j = 0; j < nc; j++)
		{
			temp.at<double>(0, 0) = 1.0;
			temp.at<double>(0, 1) = i;
			temp.at<double>(0, 2) = j;
			temp.at<double>(0, 3) = i * j;
			temp.at<double>(0, 4) = i * i;
			temp.at<double>(0, 5) = j * j;
			unwrapped_phase.at<double>(i, j) = unwrapped_phase.at<double>(i, j) + sum(temp.mul(flat_phase_coefficient))[0];
		}
	}
	//控制点绝对相位计算
	double lambda = 3e8 / (carrier_frequency.at<double>(0, 0) + 1e-10);
	double K = 0.0;
	for (int i = 0; i < valid_row.size(); i++)
	{
		int rrr = (int)gcps.at<double>(valid_row[i], 0) - offset_row;
		int ccc = (int)gcps.at<double>(valid_row[i], 1) - offset_col;
		Mat ground(1, 3, CV_64F), llh_temp(1, 3, CV_64F);
		llh_temp.at<double>(0, 0) = gcps.at<double>(valid_row[i], 3);
		llh_temp.at<double>(0, 1) = gcps.at<double>(valid_row[i], 2);
		llh_temp.at<double>(0, 2) = gcps.at<double>(valid_row[i], 4);
		util.ell2xyz(llh_temp, ground);
		double r_main = sqrt(sum((sate1(cv::Range(rrr - 1, rrr), cv::Range(0, 3)) - ground).mul(sate1(cv::Range(rrr - 1, rrr), cv::Range(0, 3)) - ground))[0]);
		double r_slave = sqrt(sum((sate2(cv::Range(rrr - 1, rrr), cv::Range(0, 3)) - ground).mul(sate2(cv::Range(rrr - 1, rrr), cv::Range(0, 3)) - ground))[0]);
		double C = mode == 1 ? 4 * PI : 2 * PI;
		double phase_real = (r_slave - r_main) / lambda * C;
		K += ((phase_real - unwrapped_phase.at<double>(rrr - 1, ccc - 1)) / (2 * PI));
	}
	K /= (double)valid_row.size();
	unwrapped_phase = unwrapped_phase + K * 2 * PI;//相位校正

	/*
	* 反演高程
	*/

	Mat R_M(1, nc, CV_64F);
	for (int i = 0; i < nc; i++)
	{
		R_M.at<double>(0, i) = nearRange + range_spacing.at<double>(0, 0) * double(i + offset_col);
	}
	Mat ones = Mat::ones(nr, 1, CV_64F);
	R_M = ones * R_M;
	Mat R_F = R_M * 2.0 + lambda * unwrapped_phase / (2 * PI);
	Mat Satellite_M_T_Position = sate1;//主星发射位置
	Mat Satellite_S_T_Position;
	if (mode == 1)
	{
		Satellite_S_T_Position = sate2;//辅星发射位置
	}
	else
	{
		Satellite_S_T_Position = sate1;//辅星发射位置
	}
	Mat Satellite_M_R_Position = sate1;//主星接收位置
	Mat Satellite_S_R_Position = sate2;//辅星接收位置
	Mat Satellite_M = (Satellite_M_R_Position + Satellite_M_T_Position) / 2;
	Mat Vs = satev1;
	ones = Mat::ones(nr, nc, CV_64F);
	Mat P1 = ones * xyz_ground.at<double>(0, 0);
	Mat P2 = ones * xyz_ground.at<double>(0, 1);
	Mat P3 = ones * xyz_ground.at<double>(0, 2);
	Mat M_T, M_R, S_T, S_R;
	Mat f1, f2, f3;
	Mat det_Df, Df_ni11, Df_ni12, Df_ni13, Df_ni21, Df_ni22, Df_ni23;
	Mat Df11, Df12, Df13, Df21, Df22, Df23, Df31, Df32, Df33, Df_ni31, Df_ni32, Df_ni33;
	Mat delta_Rt1, delta_Rt2, delta_Rt3;
	ones = Mat::ones(1, nc, CV_64F);
	Mat temp_var, temp_var1;
	Mat fd = Mat::zeros(1, nc, CV_64F);
	for (int i = 0; i < iter_times; i++)
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
		ones = Mat::ones(nr, 1, CV_64F);
		temp_var = ones * fd;
		temp_var1 = R_M * lambda / 2.0;
		f3 = f3 + temp_var.mul(temp_var1);

		//Dff
		//第一行：f(1)的x，y，z的导数
		ones = Mat::ones(1, nc, CV_64F);
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
	delta_Rt1.release(); delta_Rt2.release(); delta_Rt3.release(); Df_ni31.release(); Df_ni32.release();
	Df_ni33.release(); Df_ni21.release(); Df_ni22.release(); Df_ni23.release(); Df_ni11.release(); Df_ni12.release();
	Df_ni13.release();
	volatile bool parallel_flag = true;
	dem.create(nr, nc, CV_64F);
	Mat lon, lat;
	//lon.create(nr, nc, CV_64F); lat.create(nr, nc, CV_64F); 
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (int i = 0; i < nr; i++)
	{
		Utils util;
		for (int j = 0; j < nc; j++)
		{
			Mat xyz, llh;
			xyz = Mat::zeros(1, 3, CV_64F);
			xyz.at<double>(0, 0) = P1.at<double>(i, j);
			xyz.at<double>(0, 1) = P2.at<double>(i, j);
			xyz.at<double>(0, 2) = P3.at<double>(i, j);
			ret = util.xyz2ell(xyz, llh);
			dem.at<double>(i, j) = llh.at<double>(0, 2);
			//lat.at<double>(i, j) = llh.at<double>(0, 0);
			//lon.at<double>(i, j) = llh.at<double>(0, 1);
		}
	}
	if (parallel_check(parallel_flag, "dem_newton_iter()", parallel_error_head)) return -1;
	//dem = dem + llh.at<double>(0, 2) - dem.at<double>(row - 1, col - 1);
	Mat error(valid_row.size(), 3, CV_64F);

	for (int i = 0; i < valid_row.size(); i++)
	{
		int r, c;
		Utils::ell2xyz(gcps.at<double>(valid_row[i], 2), gcps.at<double>(valid_row[i], 3), gcps.at<double>(valid_row[i], 4), pos);
		r = gcps.at<double>(valid_row[i], 0) - offset_row;
		c = gcps.at<double>(valid_row[i], 1) - offset_col;
		error.at<double>(i, 0) = P1.at<double>(r - 1, c - 1) - pos.x;
		error.at<double>(i, 1) = P2.at<double>(r - 1, c - 1) - pos.y;
		error.at<double>(i, 2) = P3.at<double>(r - 1, c - 1) - pos.z;
	}
	util.cvmat2bin("G:\\tmp\\error.bin", error);
	return 0;
}

