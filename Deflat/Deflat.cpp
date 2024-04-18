// Deflat.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include"..\include\Deflat.h"
#include"..\include\FormatConversion.h"
#include<direct.h>
#include<Windows.h>
#include<SensAPI.h>
#include<urlmon.h>
#pragma comment(lib,"URlmon")
#pragma comment(lib, "Sensapi.lib")
#ifdef _DEBUG
#pragma comment(lib, "FormatConversion_d.lib")
#pragma comment(lib, "ComplexMat_d.lib")
#pragma comment(lib, "Utils_d.lib")
#else
#pragma comment(lib, "FormatConversion.lib")
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Utils.lib")
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

inline bool parallel_check(volatile bool parallel_flag, const char* detail_info,
	const char* parallel_error_head)
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

inline bool parallel_flag_change(volatile bool parallel_flag, int ret)
{
	if (ret < 0)
	{
		parallel_flag = false;
		return true;
	}
	else
	{
		return false;
	}
}

Deflat::Deflat()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "DEFLAT_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "DEFLAT_DLL_ERROR: error happens when using parallel computing in function: ");
	latSpacing = 5.0 / 6000.0;
	lonSpacing = 5.0 / 6000.0;
}

Deflat::~Deflat()
{
}

int Deflat::get_xyz(double aztime, Mat& coef, Mat& pos_xyz)
{
	if (aztime < 0 ||
		coef.rows != 6 ||
		coef.cols != 6 ||
		coef.type() != CV_64F ||
		coef.channels() != 1)
	{
		fprintf(stderr, "get_xyz(): input check failed!\n\n");
		return -1;
	}
	Mat tmp = Mat::ones(1, 6, CV_64F);
	tmp.at<double>(0, 1) = aztime;
	tmp.at<double>(0, 2) = aztime * aztime;
	tmp.at<double>(0, 3) = aztime * aztime * aztime;
	tmp.at<double>(0, 4) = aztime * aztime * aztime * aztime;
	tmp.at<double>(0, 5) = aztime * aztime * aztime * aztime * aztime;
	cv::transpose(tmp, tmp);
	pos_xyz = Mat::zeros(1, 3, CV_64F);
	pos_xyz.at<double>(0, 0) = tmp.dot(coef(Range(0, coef.rows), Range(0, 1)));
	pos_xyz.at<double>(0, 1) = tmp.dot(coef(Range(0, coef.rows), Range(1, 2)));
	pos_xyz.at<double>(0, 2) = tmp.dot(coef(Range(0, coef.rows), Range(2, 3)));
	return 0;
}

int Deflat::get_vel(double aztime, Mat& coef, Mat& vel_xyz)
{
	if (aztime < 0 ||
		coef.rows != 6 ||
		coef.cols != 6 ||
		coef.type() != CV_64F ||
		coef.channels() != 1)
	{
		fprintf(stderr, "get_vel(): input check failed!\n\n");
		return -1;
	}
	Mat tmp = Mat::ones(1, 6, CV_64F);
	tmp.at<double>(0, 1) = aztime;
	tmp.at<double>(0, 2) = aztime * aztime;
	tmp.at<double>(0, 3) = aztime * aztime * aztime;
	tmp.at<double>(0, 4) = aztime * aztime * aztime * aztime;
	tmp.at<double>(0, 5) = aztime * aztime * aztime * aztime * aztime;
	cv::transpose(tmp, tmp);
	vel_xyz = Mat::zeros(1, 3, CV_64F);
	vel_xyz.at<double>(0, 0) = tmp.dot(coef(Range(0, coef.rows), Range(3, 4)));
	vel_xyz.at<double>(0, 1) = tmp.dot(coef(Range(0, coef.rows), Range(4, 5)));
	vel_xyz.at<double>(0, 2) = tmp.dot(coef(Range(0, coef.rows), Range(5, 6)));
	return 0;
}

int Deflat::get_acc(double aztime, Mat& coef, Mat& acc_xyz)
{
	if (aztime < 0 ||
		coef.rows != 6 ||
		coef.cols != 6 ||
		coef.type() != CV_64F ||
		coef.channels() != 1)
	{
		fprintf(stderr, "get_acc(): input check failed!\n\n");
		return -1;
	}
	Mat tmp = Mat::ones(1, 5, CV_64F);
	acc_xyz = Mat::zeros(1, 3, CV_64F);
	tmp.at<double>(0, 0) = 1.0;
	tmp.at<double>(0, 1) = 2.0 * aztime;
	tmp.at<double>(0, 2) = 3.0 * aztime * aztime;
	tmp.at<double>(0, 3) = 4.0 * aztime * aztime * aztime;
	tmp.at<double>(0, 4) = 5.0 * aztime * aztime * aztime * aztime;
	cv::transpose(tmp, tmp);
	acc_xyz.at<double>(0, 0) = tmp.dot(coef(Range(1, 6), Range(3, 4)));
	acc_xyz.at<double>(0, 1) = tmp.dot(coef(Range(1, 6), Range(4, 5)));
	acc_xyz.at<double>(0, 2) = tmp.dot(coef(Range(1, 6), Range(5, 6)));
	return 0;
}

int Deflat::get_satellite_aztime_NEWTON(double center_time, Mat& coef, Mat pos_xyz, double start_time, double* aztime)
{
	if (center_time < 0.0 ||
		coef.rows != 6 ||
		coef.cols != 6 ||
		coef.type() != CV_64F ||
		coef.channels() != 1 ||
		pos_xyz.type() != CV_64F ||
		pos_xyz.rows != 1 ||
		pos_xyz.cols != 3 ||
		start_time > center_time||
		aztime == NULL)
	{
		fprintf(stderr, "get_satellite_aztime_NEWTON(): input check failed!\n\n");
		return -1;
	}
	double init_time = (center_time - start_time) * 100000;
	double sol = 0;
	*aztime = init_time;
	Mat S_xyz, V_xyz, A_xyz, D_xyz;
	int ret;
	for (int iter = 0; iter < 15; iter++)
	{
		ret = get_xyz(*aztime, coef, S_xyz);
		if (return_check(ret, "get_xyz(*, *, *)", error_head)) return -1;
		ret = get_vel(*aztime, coef, V_xyz);
		if (return_check(ret, "get_vel(*, *, *)", error_head)) return -1;
		ret = get_acc(*aztime, coef, A_xyz);
		if (return_check(ret, "get_acc(*, *, *)", error_head)) return -1;
		D_xyz = pos_xyz - S_xyz;
		sol = -(D_xyz.dot(V_xyz)) / (A_xyz.dot(D_xyz) - V_xyz.dot(V_xyz) + 0.000000001);
		*aztime = *aztime + sol;
		if (fabs(sol) < 0.0000454)
		{
			break;
		}
		
	}
	return 0;
}

int Deflat::Orbit_Polyfit(Mat& Orbit, Mat& coef)
{
	if (Orbit.rows < 7 ||
		Orbit.cols != 7 ||
		Orbit.type() != CV_64F ||
		Orbit.channels() != 1)
	{
		fprintf(stderr, "Orbit_Polyfit(): input check failed!\n\n");
		return -1;
	}
	int nr = Orbit.rows;
	coef = Mat::zeros(6, 6, CV_64F);
	double tmp;
	Mat t = (Orbit(Range(0, Orbit.rows), Range(0, 1)) - Orbit.at<double>(0, 0))*100000;
	Mat A = Mat::zeros(nr, 6, CV_64F);
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			tmp = t.at<double>(i, 0);
			A.at<double>(i, j) = pow(tmp, j);
		}
	}
	Mat A_t, _coef;
	cv::transpose(A, A_t);
	double ret = invert(A_t * A, A);
	if (fabs(ret) < 0.0)
	{
		fprintf(stderr, "matrix is singular!\n\n");
		return -1;
	}
	for (int j = 1; j < 7; j++)
	{
		_coef = A * A_t * Orbit(Range(0, nr), Range(j, j + 1));
		for (int i = 0; i < 6; i++)
		{
			coef.at<double>(i, j - 1) = _coef.at<double>(i, 0);
		}
	}
	
	return 0;
}

int Deflat::deflat(
	Mat& phase,
	Mat& phase_deflat,
	Mat& flat_phase,
	Mat auxi,
	Mat gcps,
	Mat orbit_main,
	Mat orbit_slave,
	int mode,
	int multilook_times
)
{
	if (phase.cols < 2 ||
		phase.rows < 2 ||
		phase.type() != CV_64F ||
		phase.channels() != 1 ||
		auxi.rows != 1 ||
		auxi.cols != 5 ||
		auxi.type() != CV_64F ||
		auxi.channels() != 1 ||
		(mode == 1 || mode == 2) == false ||
		gcps.rows < 2 ||
		gcps.cols != 5 ||
		orbit_main.rows < 1 ||
		orbit_main.cols != 7 ||
		orbit_slave.cols != 7 ||
		orbit_slave.rows < 1 ||
		orbit_main.type() != CV_64F ||
		orbit_main.channels() != 1 ||
		orbit_slave.type() != CV_64F ||
		orbit_slave.channels() != 1||
		multilook_times < 1)
	{
		fprintf(stderr, "deflat(): input check failed!\n\n");
		return -1;
	}
	int ret;
	double C = 2 * 3.1415926535;
	double lambda = 300000000.0 / (auxi.at<double>(0, 4));
	if (mode == 1) C = 4 * 3.1415926535;
	else
	{
		C = 2.0 * 3.1415926535;
	}
	if (multilook_times > 1)
	{
		for (int i = 0; i < gcps.rows; i++)
		{
			gcps.at<double>(i, 0) = gcps.at<double>(i, 0) / multilook_times + 1;
			gcps.at<double>(i, 1) = gcps.at<double>(i, 1) / multilook_times + 1;
		}
	}
	///////////////////////////////////机载无轨道数据情况//////////////////////////////////////////
	if (orbit_main.rows == 1 && orbit_slave.rows == 1 && gcps.rows == 2)
	{
		double delta_r1 = cv::norm(orbit_main(Range(0, 1), Range(1, 4)) - gcps(Range(0, 1), Range(2, 5))) - 
			cv::norm(orbit_slave(Range(0, 1), Range(1, 4)) - gcps(Range(0, 1), Range(2, 5)));
		double delta_r2 = cv::norm(orbit_main(Range(0, 1), Range(1, 4)) - gcps(Range(1, 2), Range(2, 5))) -
			cv::norm(orbit_slave(Range(0, 1), Range(1, 4)) - gcps(Range(1, 2), Range(2, 5)));
		double a = C * (delta_r2 - delta_r1) / lambda / (gcps.at<double>(1, 1) - gcps.at<double>(0, 1));
		double b = C * delta_r1 / lambda - a * gcps.at<double>(0, 1);
		flat_phase = Mat::zeros(phase.rows, phase.cols, CV_64F);
		int nr = phase.rows;
		int nc = phase.cols;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				flat_phase.at<double>(i, j) = a * double(i) + b;
			}
		}
		Utils util;
		phase_deflat = phase - flat_phase;
		ret = util.wrap(phase_deflat, phase_deflat);
		if (return_check(ret, "util.wrap(*, *)", error_head)) return -1;
		//ret = util.wrap(flat_phase, flat_phase);
		if (return_check(ret, "util.wrap(*, *)", error_head)) return -1;

	}
	///////////////////////////////////机载无轨道数据情况//////////////////////////////////////////

	///////////////////////////////////星载有轨道数据情况//////////////////////////////////////////
	else
	{
		Mat coef_m, coef_s;
		ret = Orbit_Polyfit(orbit_main, coef_m);
		if (return_check(ret, "Orbit_Polyfit(*, *)", error_head)) return -1;
		ret = Orbit_Polyfit(orbit_slave, coef_s);
		if (return_check(ret, "Orbit_Polyfit(*, *)", error_head)) return -1;
		int N_gcps = 2;
		Mat Satemain_xyz = Mat::zeros(N_gcps, 3, CV_64F);
		Mat Sateslave_xyz = Mat::zeros(N_gcps, 3, CV_64F);
		Mat tmp_xyz = Mat::zeros(1, 3, CV_64F);
		double aztime = 0.0;
		double temp1,temp2;
		Mat tmp;
		for (int i = 0; i < N_gcps; i++)
		{
			temp1 = orbit_main.at<double>(int(orbit_main.rows / 2), 0);
			temp2 = orbit_main.at<double>(0, 0);
			tmp = gcps(Range(i, i + 1), Range(2, 5));
			ret = get_satellite_aztime_NEWTON(temp1,coef_m, tmp, temp2, &aztime);
			if (return_check(ret, "get_satellite_aztime_NEWTON(*, *, *, *, *)", error_head)) return -1;
			ret = get_xyz(aztime, coef_m, tmp_xyz);
			if (return_check(ret, "get_xyz(*, *, *)", error_head)) return -1;
			tmp_xyz.copyTo(Satemain_xyz(Range(i, i + 1), Range(0, 3)));

			ret = get_satellite_aztime_NEWTON(orbit_slave.at<double>(int(orbit_slave.rows / 2), 0),coef_s, gcps(Range(i, i + 1), Range(2, 5)), orbit_slave.at<double>(0, 0), &aztime);
			if (return_check(ret, "get_satellite_aztime_NEWTON(*, *, *, *, *)", error_head)) return -1;
			ret = get_xyz(aztime, coef_s, tmp_xyz);
			if (return_check(ret, "get_xyz(*, *, *)", error_head)) return -1;
			tmp_xyz.copyTo(Sateslave_xyz(Range(i, i + 1), Range(0, 3)));
		}
		double delta_r1 = cv::norm(Satemain_xyz(Range(0, 1), Range(0, 3)) - gcps(Range(0, 1), Range(2, 5))) -
			cv::norm(Sateslave_xyz(Range(0, 1), Range(0, 3)) - gcps(Range(0, 1), Range(2, 5)));
		double delta_r2 = cv::norm(Satemain_xyz(Range(1, 2), Range(0, 3)) - gcps(Range(1, 2), Range(2, 5))) -
			cv::norm(Sateslave_xyz(Range(1, 2), Range(0, 3)) - gcps(Range(1, 2), Range(2, 5)));
		double a = C * (delta_r2 - delta_r1) / lambda / (gcps.at<double>(1, 1) - gcps.at<double>(0, 1));
		double b = C * delta_r1 / lambda - a * gcps.at<double>(0, 1);
		flat_phase = Mat::zeros(phase.rows, phase.cols, CV_64F);
		int nr = phase.rows;
		int nc = phase.cols;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				flat_phase.at<double>(i, j) = a * double(j) + b;
			}
		}
		Utils util;
		phase_deflat = phase - flat_phase;
		ret = util.wrap(phase_deflat, phase_deflat);
		if (return_check(ret, "util.wrap(*, *)", error_head)) return -1;
		//ret = util.wrap(flat_phase, flat_phase);
		if (return_check(ret, "util.wrap(*, *)", error_head)) return -1;

	}
	///////////////////////////////////星载有轨道数据情况//////////////////////////////////////////
	return 0;
}

int Deflat::deflat(
	const Mat& stateVec1,
	const Mat& stateVec2,
	const Mat& lon_coef,
	const Mat& lat_coef,
	const Mat& phase, 
	int offset_row,
	int offset_col, 
	double height,
	double time_interval, 
	double time_interval2,
	int mode,
	double wave_length,
	Mat& phase_deflated,
	Mat& flat_phase_coef
)
{
	if (stateVec1.cols != 7 ||
		stateVec2.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec2.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.type() != CV_64F ||
		lon_coef.cols != 32 ||
		lat_coef.cols != 32 ||
		lon_coef.rows != 1 ||
		lat_coef.rows != 1 ||
		lon_coef.type() != CV_64F ||
		lat_coef.type() != CV_64F ||
		phase.empty() ||
		phase.type() != CV_64F ||
		phase.cols * phase.rows < 500||
		height < 0.0||
		time_interval < 0.0||
		time_interval2 < 0.0 ||
		mode > 2||
		mode < 1||
		wave_length <= 0.0
		)
	{
		fprintf(stderr, "deflat(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion;
	Utils util;
	int ret;
	int rows = phase.rows; int cols = phase.cols;
	/*
	* 轨道插值
	*/
	Mat state_vec1, state_vec2;
	stateVec1.copyTo(state_vec1);
	stateVec2.copyTo(state_vec2);
	ret = util.stateVec_interp(state_vec1, time_interval, state_vec1);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	ret = util.stateVec_interp(state_vec2, time_interval2, state_vec2);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;

	Mat sate1_xyz, sate2_xyz, sate1_v, sate2_v;
	state_vec1(cv::Range(0, state_vec1.rows), cv::Range(1, 4)).copyTo(sate1_xyz);
	state_vec1(cv::Range(0, state_vec1.rows), cv::Range(4, 7)).copyTo(sate1_v);
	state_vec2(cv::Range(0, state_vec2.rows), cv::Range(1, 4)).copyTo(sate2_xyz);
	state_vec2(cv::Range(0, state_vec2.rows), cv::Range(4, 7)).copyTo(sate2_v);
	/*
	* 图像坐标转经纬坐标
	*/
	Mat row, col;
	row.create(rows, cols, CV_64F); col.create(rows, cols, CV_64F);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			row.at<double>(i, j) = i + offset_row;
			col.at<double>(i, j) = j + offset_col;
		}
	}
	Mat lon, lat, lon_coefficient, lat_coefficient;
	lon_coef.copyTo(lon_coefficient);
	lat_coef.copyTo(lat_coefficient);
	ret = util.coord_conversion(lon_coefficient, row, col, lon);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	ret = util.coord_conversion(lat_coefficient, row, col, lat);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;


	/*
	* 图像1成像点位置计算
	*/

	Mat sate1 = Mat::zeros(rows, 3, CV_64F);
	Mat sate2 = Mat::zeros(rows, 3, CV_64F);

	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = height;
		util.ell2xyz(tmp, xyz);
		
		//找到零多普勒位置
		Mat dop = Mat::zeros(sate1_xyz.rows, 1, CV_64F);
		Mat r;
		for (int j = 0; j < sate1_xyz.rows; j++)
		{
			r = xyz - sate1_xyz(Range(j, j + 1), Range(0, 3));
			dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate1_v(Range(j, j + 1), Range(0, 3))))[0]);
		}
		Point peak_loc;
		cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
		//sate1_xyz(Range(peak_loc.y, peak_loc.y + 1), Range(0, 3)).copyTo(sate1(Range(i, i + 1), Range(0, 3)));
		int xxxx;
		for (int j = 0; j < rows; j++)
		{
			xxxx = (peak_loc.y + j) > (sate1_xyz.rows - 1) ? (sate1_xyz.rows - 1) : (peak_loc.y + j);
			sate1_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate1(Range(j, j + 1), Range(0, 3)));
		}
	}

	/*
	* 图像2成像点位置计算
	*/

	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = height;
		util.ell2xyz(tmp, xyz);

		//找到零多普勒位置
		Mat dop = Mat::zeros(sate2_xyz.rows, 1, CV_64F);
		Mat r;
		for (int j = 0; j < sate2_xyz.rows; j++)
		{
			r = xyz - sate2_xyz(Range(j, j + 1), Range(0, 3));
			dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate2_v(Range(j, j + 1), Range(0, 3))))[0]);
		}
		Point peak_loc;
		cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
		//sate2_xyz(Range(peak_loc.y, peak_loc.y + 1), Range(0, 3)).copyTo(sate2(Range(i, i + 1), Range(0, 3)));
		int xxxx;
		for (int j = 0; j < rows; j++)
		{
			xxxx = (peak_loc.y + j) > (sate2_xyz.rows - 1) ? (sate2_xyz.rows - 1) : (peak_loc.y + j);
			sate2_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate2(Range(j, j + 1), Range(0, 3)));
		}
	}


	/*
	* 计算斜距和平地相位
	*/

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat r, xyz; double r1, r2;
		for (int j = 0; j < cols; j++)
		{
			tmp.at<double>(0, 0) = lat.at<double>(i, j);
			tmp.at<double>(0, 1) = lon.at<double>(i, j);
			tmp.at<double>(0, 2) = height;
			util.ell2xyz(tmp, xyz);
			r = xyz - sate1(Range(i, i + 1), Range(0, 3));
			r1 = sqrt(sum(r.mul(r))[0]);//主星斜距
			r = xyz - sate2(Range(i, i + 1), Range(0, 3));
			r2 = sqrt(sum(r.mul(r))[0]);//辅星斜距
			lat.at<double>(i, j) = (r2 - r1) / wave_length * (1 / (double)mode) * 4 * PI;
		}
	}

	/*
	* 拟合平地相位（二阶拟合），取1/20进行拟合
	*/
	int rows1 = (int)floor(rows * cols / 20);
	Mat A = Mat::ones(rows1, 6, CV_64F); Mat ro(rows1, 1, CV_64F); Mat co(rows1, 1, CV_64F); Mat delta_phi(rows1, 1, CV_64F);
	lat = lat.reshape(0, cols* rows);
	row = row - offset_row; col = col - offset_col;
	row = row.reshape(0, cols * rows); col = col.reshape(0, cols * rows);
	for (int i = 0; i < rows1; i++)
	{
		ro.at<double>(i, 0) = row.at<double>(i * 20, 0);
		co.at<double>(i, 0) = col.at<double>(i * 20, 0);
		delta_phi.at<double>(i, 0) = lat.at<double>(i * 20, 0);
	}
	row.release(); col.release();
	ro.copyTo(A(Range(0, rows1), Range(1, 2)));
	co.copyTo(A(Range(0, rows1), Range(2, 3)));
	Mat tmp = ro.mul(co);
	tmp.copyTo(A(Range(0, rows1), Range(3, 4)));
	tmp = ro.mul(ro);
	tmp.copyTo(A(Range(0, rows1), Range(4, 5)));
	tmp = co.mul(co);
	tmp.copyTo(A(Range(0, rows1), Range(5, 6)));

	Mat coef, A_t, b;
	cv::transpose(A, A_t);
	b = A_t * delta_phi;
	A = A_t * A;
	if (!cv::solve(A, b, coef, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "deflat(): matrix deficiency!");
		return -1;
	}

	//求平地相位
	cv::transpose(coef, coef);
	coef.copyTo(flat_phase_coef);
	lat = lat.reshape(0, rows);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		Mat temp(1, 6, CV_64F);
		for (int j = 0; j < cols; j++)
		{
			temp.at<double>(0, 0) = 1.0;
			temp.at<double>(0, 1) = i;
			temp.at<double>(0, 2) = j;
			temp.at<double>(0, 3) = i * j;
			temp.at<double>(0, 4) = i * i;
			temp.at<double>(0, 5) = j * j;
			lat.at<double>(i, j) = sum(temp.mul(coef))[0];
		}
	}
	lat = phase - lat;
	util.wrap(lat, lat);
	lat.copyTo(phase_deflated);
	return 0;
}


int Deflat::topo_removal(
	const Mat& phase,
	const Mat& dem, 
	const Mat& dem_range_lon,
	const Mat& dem_range_lat, 
	const Mat& stateVec1, 
	const Mat& stateVec2, 
	const Mat& lon_coef,
	const Mat& lat_coef,
	const Mat& inc_coef,
	int offset_row, 
	int offset_col,
	double interp_interval1,
	double interp_interval2,
	int mode,
	double wavelength,
	double B_effect,
	Mat& phase_detopo
)
{
	if (phase.empty() ||
		phase.type() != CV_64F ||
		dem.empty() ||
		dem_range_lon.empty() ||
		dem_range_lat.empty() ||
		dem_range_lat.rows != 1 || dem_range_lat.cols != 2 || dem_range_lon.rows != 1 || dem_range_lon.cols != 2 ||
		dem_range_lat.type() != CV_64F || dem_range_lon.type() != CV_64F ||
		stateVec1.cols != 7 || stateVec1.rows < 7 || stateVec1.type() != CV_64F || stateVec2.cols != 7 || stateVec2.rows < 7 || stateVec2.type() != CV_64F ||
		lon_coef.rows != 1 || lon_coef.cols != 32 || lon_coef.type() != CV_64F || lat_coef.rows != 1 || lat_coef.cols != 32 || lat_coef.type() != CV_64F ||
		interp_interval1 < 0.0 || interp_interval2 < 0.0 ||
		mode > 2 || mode < 1 || wavelength < 0.0 || inc_coef.type() != CV_64F||inc_coef.rows != 1|| inc_coef.cols != 11
		)
	{
		fprintf(stderr, "topo_removal(): input check failed!\n");
		return -1;
	}

	/*
	* 轨道插值
	*/

	FormatConversion conversion;
	Utils util;
	int ret;
	int rows = phase.rows; int cols = phase.cols;
	/*
	* 轨道插值
	*/
	Mat state_vec1, state_vec2;
	stateVec1.copyTo(state_vec1);
	stateVec2.copyTo(state_vec2);
	ret = util.stateVec_interp(state_vec1, interp_interval1, state_vec1);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	ret = util.stateVec_interp(state_vec2, interp_interval2, state_vec2);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;

	Mat sate1_xyz, sate2_xyz, sate1_v, sate2_v;
	state_vec1(cv::Range(0, state_vec1.rows), cv::Range(1, 4)).copyTo(sate1_xyz);
	state_vec1(cv::Range(0, state_vec1.rows), cv::Range(4, 7)).copyTo(sate1_v);
	state_vec2(cv::Range(0, state_vec2.rows), cv::Range(1, 4)).copyTo(sate2_xyz);
	state_vec2(cv::Range(0, state_vec2.rows), cv::Range(4, 7)).copyTo(sate2_v);
	/*
	* 图像坐标转经纬坐标
	*/
	Mat row, col;
	row.create(rows, cols, CV_64F); col.create(rows, cols, CV_64F);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			row.at<double>(i, j) = i + offset_row;
			col.at<double>(i, j) = j + offset_col;
		}
	}
	Mat lon, lat, lon_coefficient, lat_coefficient;
	lon_coef.copyTo(lon_coefficient);
	lat_coef.copyTo(lat_coefficient);
	ret = util.coord_conversion(lon_coefficient, row, col, lon);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	ret = util.coord_conversion(lat_coefficient, row, col, lat);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;

	/*
	* 图像1成像点位置计算
	*/

	Mat sate1 = Mat::zeros(rows, 3, CV_64F);
	Mat sate2 = Mat::zeros(rows, 3, CV_64F);
	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = 0;
		util.ell2xyz(tmp, xyz);

		//找到零多普勒位置
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
		for (int j = 0; j < rows; j++)
		{
			xxxx = (peak_loc.y + j) > (sate1_xyz.rows - 1) ? (sate1_xyz.rows - 1) : (peak_loc.y + j);
			sate1_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate1(Range(j, j + 1), Range(0, 3)));
		}
	}

	/*
	* 图像2成像点位置计算
	*/

	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = 0;
		util.ell2xyz(tmp, xyz);

		//找到零多普勒位置
		Mat dop = Mat::zeros(sate2_xyz.rows, 1, CV_64F);
		Mat r;
		for (int j = 0; j < sate2_xyz.rows; j++)
		{
			r = xyz - sate2_xyz(Range(j, j + 1), Range(0, 3));
			dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate2_v(Range(j, j + 1), Range(0, 3))))[0]);
		}
		Point peak_loc;
		cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
		int xxxx;
		for (int j = 0; j < rows; j++)
		{
			xxxx = (peak_loc.y + j) > (sate2_xyz.rows - 1) ? (sate2_xyz.rows - 1) : (peak_loc.y + j);
			sate2_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate2(Range(j, j + 1), Range(0, 3)));
		}
	}

	/*
	* 根据外部DEM数据计算相应的地形相位
	*/
	double lon_left, lon_right, lat_top, lat_bottom, delta_lon_dem, delta_lat_dem,
		offset_inc, scale_inc, offset_y, scale_y, a0, a1, a2, a3, a4, a5;
	offset_inc = inc_coef.at<double>(0, 0);
	scale_inc = inc_coef.at<double>(0, 1);
	offset_y = inc_coef.at<double>(0, 2);
	scale_y = inc_coef.at<double>(0, 3);
	a0 = inc_coef.at<double>(0, 4);
	a1 = inc_coef.at<double>(0, 5);
	a2 = inc_coef.at<double>(0, 6);
	a3 = inc_coef.at<double>(0, 7);
	a4 = inc_coef.at<double>(0, 8);
	a5 = inc_coef.at<double>(0, 9);
	int dem_rows = dem.rows; int dem_cols = dem.cols;
	lon_left = dem_range_lon.at<double>(0, 0); lon_right = dem_range_lon.at<double>(0, 1);
	lat_top = dem_range_lat.at<double>(0, 1); lat_bottom = dem_range_lat.at<double>(0, 0);
	delta_lon_dem = fabs(lon_right - lon_left) > 180.0 ? (lon_right - lon_left + 360.0) / (double)dem_cols : (lon_right - lon_left) / (double)dem_cols;
	delta_lat_dem = (lat_top - lat_bottom) / (double)dem_rows;
	Mat DEM; DEM.create(lat.size(), CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat r, xyz; 
		double r1, r2, latitude, longtitude, inc, jj, delta_x, delta_y, upper, lower, dem_interp;
		int m, n, m_p1, n_p1;
		for (int j = 0; j < cols; j++)
		{
			jj = (double)j + (double)offset_col;
			jj = (jj - offset_y) / scale_y;
			inc = a0 + a1 * jj + a2 * jj * jj + a3 * jj * jj * jj + a4 * jj * jj * jj * jj + a5 * jj * jj * jj * jj * jj;
			inc = inc * scale_inc + offset_inc;
			latitude = lat.at<double>(i, j);
			longtitude = lon.at<double>(i, j);
			n = fabs(longtitude - lon_left) > 180.0 ? (int)floor((longtitude - lon_left + 360.0) / delta_lon_dem) : (int)floor((longtitude - lon_left) / delta_lon_dem);
			m = (int)floor((lat_top - latitude) / delta_lat_dem);
			m = m > (dem_rows - 1) ? (dem_rows - 1) : m;
			n = n > (dem_cols - 1) ? (dem_cols - 1) : n;
			m_p1 = m + 1; n_p1 = n + 1;
			m_p1 = m_p1 > (dem_rows - 1) ? (dem_rows - 1) : m_p1;
			n_p1 = n_p1 > (dem_cols - 1) ? (dem_cols - 1) : n_p1;
			delta_x = fabs(delta_lon_dem * (double)n + lon_left) > 180.0 ? (delta_lon_dem * (double)n + lon_left + 360.0) : (delta_lon_dem * (double)n + lon_left);
			delta_x = longtitude - delta_x;
			delta_y = lat_top - (delta_lat_dem * (double)m + latitude);
			upper = (double)dem.at<short>(m, n) + delta_x / delta_lon_dem * ((double)dem.at<short>(m, n_p1) - (double)dem.at<short>(m, n));
			lower = (double)dem.at<short>(m_p1, n) + delta_x / delta_lon_dem * ((double)dem.at<short>(m_p1, n_p1) - (double)dem.at<short>(m_p1, n));
			dem_interp = upper + delta_y / delta_lat_dem * (lower - upper);

			tmp.at<double>(0, 0) = lat.at<double>(i, j);
			tmp.at<double>(0, 1) = lon.at<double>(i, j);
			tmp.at<double>(0, 2) = dem_interp;
			DEM.at<double>(i, j) = dem_interp;
			util.ell2xyz(tmp, xyz);
			r = xyz - sate1(Range(i, i + 1), Range(0, 3));
			r1 = sqrt(sum(r.mul(r))[0]);//主星斜距
			lat.at<double>(i, j) =  4 * PI * dem_interp * B_effect / wavelength / r1 / sin(inc / 180.0 * PI) / ((double)mode);
		}
	}
	lat = phase + lat;
	util.wrap(lat, lat);
	lat.copyTo(phase_detopo);
	return 0;
}

int Deflat::topography_simulation(
	Mat& topography_phase, 
	Mat& statevector1,
	Mat& statevector2, 
	Mat& lon_cofficient,
	Mat& lat_cofficient,
	Mat& inc_cofficient, 
	double prf1, double prf2, 
	int sceneHeight, int sceneWidth,
	int offset_row, int offset_col,
	double nearRangeTime, double rangeSpacing, double wavelength,
	double acquisition_start_time, double acquisition_stop_time, 
	const char* DEMpath,
	int interp_times
)
{
	if (wavelength < 0.0 ||
		inc_cofficient.type() != CV_64F ||
		inc_cofficient.rows != 1 ||
		inc_cofficient.cols != 11 ||
		nearRangeTime <= 0.0 ||
		prf1 <= 0.0 ||
		prf2 <= 0.0 ||
		wavelength < 0.0 ||
		rangeSpacing < 0.0 ||
		statevector1.type() != CV_64F ||
		statevector1.rows < 5 ||
		statevector1.cols != 7 ||
		statevector2.type() != CV_64F ||
		statevector2.rows < 5 ||
		statevector2.cols != 7 ||
		acquisition_start_time <= 0 ||
		acquisition_stop_time <= 0 ||
		sceneHeight < 1 ||
		sceneWidth < 1 ||
		!DEMpath
		)
	{
		fprintf(stderr, "topography_simulation(): input check failed!\n");
		return -1;
	}
	int ret;
	Mat dem, dem_out;
	Utils util;
	double lon_upperleft, lat_upperleft, lonMax, lonMin, latMax, latMin, B_effect, B_para;
	ret = computeImageGeoBoundry(lat_cofficient, lon_cofficient, sceneHeight, sceneWidth, offset_row, offset_col,
		&lonMax, &latMax, &lonMin, &latMin);
	if (return_check(ret, "computeImageGeoBoundry()", error_head)) return -1;
	ret = getSRTMDEM(DEMpath, dem, &lon_upperleft, &lat_upperleft, lonMin, lonMax, latMin, latMax);
	if (return_check(ret, "getSRTMDEM()", error_head)) return -1;
	ret = demMapping(dem, dem_out, lon_upperleft, lat_upperleft, offset_row, offset_col,
		sceneHeight, sceneWidth, prf1, rangeSpacing, wavelength,
		nearRangeTime, acquisition_start_time, acquisition_stop_time, statevector1, interp_times);
	if (return_check(ret, "demMapping()", error_head)) return -1;
	//Mat out; dem_out.convertTo(out, CV_64F);
	//util.cvmat2bin("E:\\zgb1\\functions\\out.bin", out);
	ret = util.baseline_estimation(statevector1, statevector2, lon_cofficient, lat_cofficient, offset_row,
		offset_col, dem_out.rows, dem_out.cols,
		1 / prf1, 1 / prf2, &B_effect, &B_para);
	if (return_check(ret, "baseline_estimation()", error_head)) return -1;
	ret = topography_phase_simulation(dem_out, topography_phase, inc_cofficient, B_effect, nearRangeTime,
		offset_row, offset_col, wavelength, rangeSpacing);
	if (return_check(ret, "topography_phase_simulation()", error_head)) return -1;
	return 0;
}

int Deflat::demMapping(
	Mat& DEM84,
	Mat& mappedDEM,
	double lon_upperleft,
	double lat_upperleft,
	int offset_row,
	int offset_col,
	int sceneHeight,
	int sceneWidth,
	double prf, 
	double rangeSpacing,
	double wavelength,
	double nearRangeTime,
	double acquisitionStartTime,
	double acquisitionStopTime,
	Mat& stateVector,
	int interp_times,
	double lon_spacing,
	double lat_spacing
)
{
	if (DEM84.empty() ||
		DEM84.type() != CV_16S ||
		sceneHeight < 10 ||
		sceneWidth < 10 ||
		prf <= 0 ||
		wavelength <= 0 ||
		rangeSpacing <= 0 ||
		nearRangeTime <= 0 ||
		acquisitionStartTime <= 0 ||
		acquisitionStopTime <= 0 ||
		lon_spacing <= 0.0 ||
		lat_spacing <= 0.0 ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		stateVector.type() != CV_64F ||
		stateVector.rows < 5 ||
		stateVector.cols != 7
		)
	{
		fprintf(stderr, "demMapping(): input check failed!\n");
		return -1;
	}

	//84坐标系DEM插值
	Mat DEM, stateVector_interp;
	interp_times = interp_times < 1 ? 1 : interp_times;
	cv::resize(DEM84, DEM, cv::Size(DEM84.cols * interp_times, DEM84.rows * interp_times));
	Mat DEM_out = Mat::zeros(sceneHeight, sceneWidth, CV_16S);
	short invalid = -999;
	DEM_out = DEM_out + invalid;
	lon_spacing = lon_spacing / (double)interp_times;
	lat_spacing = lat_spacing / (double)interp_times;
	//初始化轨道类
	double delta_t = stateVector.at<double>(1, 0) - stateVector.at<double>(0, 0);
	orbitStateVectors stateVectors(stateVector, acquisitionStartTime, acquisitionStopTime, delta_t);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;

	int DEM_rows = DEM.rows; int DEM_cols = DEM.cols;
	double dopplerFrequency = 0.0;
	//采用迭代计算每个DEM点在SAR图像中的坐标，以减小计算量
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < DEM_rows; i++)
	{
		for (int j = 0; j < DEM_cols; j++)
		{
			Position groundPosition;
			double lat, lon, height;
			lat = lat_upperleft - (double)i * lat_spacing;
			lon = lon_upperleft + (double)j * lon_spacing;
			lon = lon > 180.0 ? (lon - 360.0) : lon;
			height = DEM.at<short>(i, j);
			Utils::ell2xyz(lon, lat, height, groundPosition);
			int numOrbitVec = stateVectors.newStateVectors.rows;
			double firstVecTime = 0.0;
			double secondVecTime = 0.0;
			double firstVecFreq = 0.0;
			double secondVecFreq = 0.0;
			double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
			for (int ii = 0; ii < numOrbitVec; ii++) {
				Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
					stateVectors.newStateVectors.at<double>(ii, 3));
				Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
					stateVectors.newStateVectors.at<double>(ii, 6));
				currentFreq = 0;
				xdiff = groundPosition.x - orb_pos.x;
				ydiff = groundPosition.y - orb_pos.y;
				zdiff = groundPosition.z - orb_pos.z;
				distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
				currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
				if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
					firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
					firstVecFreq = currentFreq;
				}
				else {
					secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
					secondVecFreq = currentFreq;
					break;
				}
			}

			if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
				continue;
			}

			double lowerBoundTime = firstVecTime;
			double upperBoundTime = secondVecTime;
			double lowerBoundFreq = firstVecFreq;
			double upperBoundFreq = secondVecFreq;
			double midTime, midFreq;
			double diffTime = fabs(upperBoundTime - lowerBoundTime);
			double absLineTimeInterval = time_interval;

			int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
			int numIterations = 0; Position pos; Velocity vel;
			while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

				midTime = (upperBoundTime + lowerBoundTime) / 2.0;
				stateVectors.getPosition(midTime, pos);
				stateVectors.getVelocity(midTime, vel);
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
					zeroDopplerTime =  midTime;
					break;
				}

				diffTime = fabs(upperBoundTime - lowerBoundTime);
				numIterations++;
			}


			zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
			int azimuthIndex = (zeroDopplerTime - acquisitionStartTime) / time_interval;
			int rangeIndex = (distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing;
			azimuthIndex = azimuthIndex - offset_row;
			rangeIndex = rangeIndex - offset_col;
			if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
			{
				
			}
			else
			{
				DEM_out.at<short>(azimuthIndex, rangeIndex) = DEM.at<short>(i, j);
			}
		}
	}
	
	//投影DEM插值
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			if (DEM_out.at<short>(i, j) != invalid) continue;
			int up, down, left, right, up_count, down_count, left_count, right_count;
			double value1, value2, ratio1, ratio2;
			//寻找上面有值的点
			up = i;
			while (true)
			{
				up--;
				if (up < 0) break;
				if (DEM_out.at<short>(up, j) != invalid) break;
			}
			//寻找下面有值的点
			down = i;
			while (true)
			{
				down++;
				if (down > sceneHeight - 1) break;
				if (DEM_out.at<short>(down, j) != invalid) break;
			}
			//寻找左边有值的点
			left = j;
			while (true)
			{
				left--;
				if (left < 0) break;
				if (DEM_out.at<short>(i, left) != invalid) break;
			}
			//寻找右边有值的点
			right = j;
			while (true)
			{
				right++;
				if (right > sceneWidth - 1) break;
				if (DEM_out.at<short>(i, right) != invalid) break;
			}

			//上下左右都有值
			if (left >= 0 && right <= sceneWidth - 1 && up >= 0 && down <= sceneHeight - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(DEM_out.at<short>(i, left)) + 
					double(DEM_out.at<short>(i, right) - DEM_out.at<short>(i, left)) * ratio1;
				ratio2 = double(i - up) / double(down - up);
				value2 = double(DEM_out.at<short>(up, j)) +
					double(DEM_out.at<short>(down, j) - DEM_out.at<short>(up, j)) * ratio2;
				DEM_out.at<short>(i, j) = (value1 + value2) / 2.0;
				continue;
			}
			//上下有值
			if (up >= 0 && down <= sceneHeight - 1)
			{
				ratio2 = double(i - up) / double(down - up);
				value2 = double(DEM_out.at<short>(up, j)) +
					double(DEM_out.at<short>(down, j) - DEM_out.at<short>(up, j)) * ratio2;
				DEM_out.at<short>(i, j) = value2;
				continue;
			}
			//左右有值
			if (left >= 0 && right <= sceneWidth - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(DEM_out.at<short>(i, left)) +
					double(DEM_out.at<short>(i, right) - DEM_out.at<short>(i, left)) * ratio1;
				DEM_out.at<short>(i, j) = value1;
				continue;
			}
			//上边有值
			if (up >= 0)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(up, j);
				continue;
			}
			//下边有值
			if (down <= sceneHeight - 1)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(down, j);
				continue;
			}
			//左边有值
			if (left >= 0)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(i, left);
				continue;
			}
			//右边有值
			if (right <= sceneWidth - 1)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(i, right);
				continue;
			}
			//上下左右都没有值
			DEM_out.at<short>(i, j) = 0;
			
		}
	}
	cv::GaussianBlur(DEM_out, mappedDEM, cv::Size(5, 5), 1, 1);
	return 0;
}

int Deflat::demMapping(
	Mat& DEM84,
	Mat& mappedDEM,
	Mat& mappedLat,
	Mat& mappedLon,
	double lon_upperleft, 
	double lat_upperleft, 
	int offset_row,
	int offset_col,
	int sceneHeight,
	int sceneWidth,
	double prf,
	double rangeSpacing, 
	double wavelength,
	double nearRangeTime,
	double acquisitionStartTime,
	double acquisitionStopTime, 
	Mat& stateVector, 
	int interp_times, 
	double lon_spacing,
	double lat_spacing
)
{
	if (DEM84.empty() ||
		DEM84.type() != CV_16S ||
		sceneHeight < 10 ||
		sceneWidth < 10 ||
		prf <= 0 ||
		wavelength <= 0 ||
		rangeSpacing <= 0 ||
		nearRangeTime <= 0 ||
		acquisitionStartTime <= 0 ||
		acquisitionStopTime <= 0 ||
		lon_spacing <= 0.0 ||
		lat_spacing <= 0.0 ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		stateVector.type() != CV_64F ||
		stateVector.rows < 5 ||
		stateVector.cols != 7
		)
	{
		fprintf(stderr, "demMapping(): input check failed!\n");
		return -1;
	}

	//84坐标系DEM插值
	Mat DEM, stateVector_interp;
	interp_times = interp_times < 1 ? 1 : interp_times;
	cv::resize(DEM84, DEM, cv::Size(DEM84.cols * interp_times, DEM84.rows * interp_times));
	Mat DEM_out = Mat::zeros(sceneHeight, sceneWidth, CV_16S);
	//Mat temp_lonlat = Mat::zeros(sceneHeight, sceneWidth, CV_32F);
	mappedLat.create(sceneHeight, sceneWidth, CV_64F); mappedLon.create(sceneHeight, sceneWidth, CV_64F);
	mappedLat = -999.0; mappedLon = -999.0;
	//temp_lonlat = temp_lonlat - 999.0;
	//temp_lonlat.copyTo(mappedLat); temp_lonlat.copyTo(mappedLon);
	short invalid = -999;
	DEM_out = DEM_out + invalid;
	lon_spacing = lon_spacing / (double)interp_times;
	lat_spacing = lat_spacing / (double)interp_times;
	//初始化轨道类
	double delta_t = stateVector.at<double>(1, 0) - stateVector.at<double>(0, 0);
	orbitStateVectors stateVectors(stateVector, acquisitionStartTime, acquisitionStopTime, delta_t);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;

	int DEM_rows = DEM.rows; int DEM_cols = DEM.cols;
	double dopplerFrequency = 0.0;
	//采用迭代计算每个DEM点在SAR图像中的坐标，以减小计算量
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < DEM_rows; i++)
	{
		for (int j = 0; j < DEM_cols; j++)
		{
			Position groundPosition;
			double lat, lon, height;
			lat = lat_upperleft - (double)i * lat_spacing;
			lon = lon_upperleft + (double)j * lon_spacing;
			lon = lon > 180.0 ? (lon - 360.0) : lon;
			height = DEM.at<short>(i, j);
			Utils::ell2xyz(lon, lat, height, groundPosition);
			int numOrbitVec = stateVectors.newStateVectors.rows;
			double firstVecTime = 0.0;
			double secondVecTime = 0.0;
			double firstVecFreq = 0.0;
			double secondVecFreq = 0.0;
			double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
			for (int ii = 0; ii < numOrbitVec; ii++) {
				Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
					stateVectors.newStateVectors.at<double>(ii, 3));
				Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
					stateVectors.newStateVectors.at<double>(ii, 6));
				currentFreq = 0;
				xdiff = groundPosition.x - orb_pos.x;
				ydiff = groundPosition.y - orb_pos.y;
				zdiff = groundPosition.z - orb_pos.z;
				distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
				currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
				if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
					firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
					firstVecFreq = currentFreq;
				}
				else {
					secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
					secondVecFreq = currentFreq;
					break;
				}
			}

			if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
				continue;
			}

			double lowerBoundTime = firstVecTime;
			double upperBoundTime = secondVecTime;
			double lowerBoundFreq = firstVecFreq;
			double upperBoundFreq = secondVecFreq;
			double midTime, midFreq;
			double diffTime = fabs(upperBoundTime - lowerBoundTime);
			double absLineTimeInterval = time_interval;

			int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
			int numIterations = 0; Position pos; Velocity vel;
			while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

				midTime = (upperBoundTime + lowerBoundTime) / 2.0;
				stateVectors.getPosition(midTime, pos);
				stateVectors.getVelocity(midTime, vel);
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
			int azimuthIndex = (zeroDopplerTime - acquisitionStartTime) / time_interval;
			int rangeIndex = (distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing;
			azimuthIndex = azimuthIndex - offset_row;
			rangeIndex = rangeIndex - offset_col;
			if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
			{

			}
			else
			{
				DEM_out.at<short>(azimuthIndex, rangeIndex) = DEM.at<short>(i, j);
				mappedLon.at<double>(azimuthIndex, rangeIndex) = lon;
				mappedLat.at<double>(azimuthIndex, rangeIndex) = lat;
			}
		}
	}
	//投影DEM插值
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			if (DEM_out.at<short>(i, j) != invalid) continue;
			int up, down, left, right, up_count, down_count, left_count, right_count;
			double value1, value2, ratio1, ratio2;
			//寻找上面有值的点
			up = i;
			while (true)
			{
				up--;
				if (up < 0) break;
				if (DEM_out.at<short>(up, j) != invalid) break;
			}
			//寻找下面有值的点
			down = i;
			while (true)
			{
				down++;
				if (down > sceneHeight - 1) break;
				if (DEM_out.at<short>(down, j) != invalid) break;
			}
			//寻找左边有值的点
			left = j;
			while (true)
			{
				left--;
				if (left < 0) break;
				if (DEM_out.at<short>(i, left) != invalid) break;
			}
			//寻找右边有值的点
			right = j;
			while (true)
			{
				right++;
				if (right > sceneWidth - 1) break;
				if (DEM_out.at<short>(i, right) != invalid) break;
			}

			//上下左右都有值
			if (left >= 0 && right <= sceneWidth - 1 && up >= 0 && down <= sceneHeight - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(DEM_out.at<short>(i, left)) +
					double(DEM_out.at<short>(i, right) - DEM_out.at<short>(i, left)) * ratio1;
				ratio2 = double(i - up) / double(down - up);
				value2 = double(DEM_out.at<short>(up, j)) +
					double(DEM_out.at<short>(down, j) - DEM_out.at<short>(up, j)) * ratio2;
				DEM_out.at<short>(i, j) = (value1 + value2) / 2.0;
				continue;
			}
			//上下有值
			if (up >= 0 && down <= sceneHeight - 1)
			{
				ratio2 = double(i - up) / double(down - up);
				value2 = double(DEM_out.at<short>(up, j)) +
					double(DEM_out.at<short>(down, j) - DEM_out.at<short>(up, j)) * ratio2;
				DEM_out.at<short>(i, j) = value2;
				continue;
			}
			//左右有值
			if (left >= 0 && right <= sceneWidth - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(DEM_out.at<short>(i, left)) +
					double(DEM_out.at<short>(i, right) - DEM_out.at<short>(i, left)) * ratio1;
				DEM_out.at<short>(i, j) = value1;
				continue;
			}
			//上边有值
			if (up >= 0)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(up, j);
				continue;
			}
			//下边有值
			if (down <= sceneHeight - 1)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(down, j);
				continue;
			}
			//左边有值
			if (left >= 0)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(i, left);
				continue;
			}
			//右边有值
			if (right <= sceneWidth - 1)
			{
				DEM_out.at<short>(i, j) = DEM_out.at<short>(i, right);
				continue;
			}
			//上下左右都没有值
			DEM_out.at<short>(i, j) = 0;

		}
	}

	//投影经度插值
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			if (mappedLon.at<double>(i, j) > -998.0) continue;
			int up, down, left, right, up_count, down_count, left_count, right_count;
			double value1, value2, ratio1, ratio2;
			//寻找上面有值的点
			up = i;
			while (true)
			{
				up--;
				if (up < 0) break;
				if (mappedLat.at<double>(up, j) > -998.0) break;
			}
			//寻找下面有值的点
			down = i;
			while (true)
			{
				down++;
				if (down > sceneHeight - 1) break;
				if (mappedLat.at<double>(down, j) > -998.0) break;
			}
			//寻找左边有值的点
			left = j;
			while (true)
			{
				left--;
				if (left < 0) break;
				if (mappedLat.at<double>(i, left) > -998.0) break;
			}
			//寻找右边有值的点
			right = j;
			while (true)
			{
				right++;
				if (right > sceneWidth - 1) break;
				if (mappedLat.at<double>(i, right) > -998.0) break;
			}

			//上下左右都有值
			if (left >= 0 && right <= sceneWidth - 1 && up >= 0 && down <= sceneHeight - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(mappedLat.at<double>(i, left)) +
					double(mappedLat.at<double>(i, right) - mappedLat.at<double>(i, left)) * ratio1;
				ratio2 = double(i - up) / double(down - up);
				value2 = double(mappedLat.at<double>(up, j)) +
					double(mappedLat.at<double>(down, j) - mappedLat.at<double>(up, j)) * ratio2;
				mappedLat.at<double>(i, j) = (value1 + value2) / 2.0;


				ratio1 = double(j - left) / double(right - left);
				value1 = double(mappedLon.at<double>(i, left)) +
					double(mappedLon.at<double>(i, right) - mappedLon.at<double>(i, left)) * ratio1;
				ratio2 = double(i - up) / double(down - up);
				value2 = double(mappedLon.at<double>(up, j)) +
					double(mappedLon.at<double>(down, j) - mappedLon.at<double>(up, j)) * ratio2;
				mappedLon.at<double>(i, j) = (value1 + value2) / 2.0;

				continue;
			}
			//上下有值
			if (up >= 0 && down <= sceneHeight - 1)
			{
				ratio2 = double(i - up) / double(down - up);
				value2 = double(mappedLat.at<double>(up, j)) +
					double(mappedLat.at<double>(down, j) - mappedLat.at<double>(up, j)) * ratio2;
				mappedLat.at<double>(i, j) = value2;


				ratio2 = double(i - up) / double(down - up);
				value2 = double(mappedLon.at<double>(up, j)) +
					double(mappedLon.at<double>(down, j) - mappedLon.at<double>(up, j)) * ratio2;
				mappedLon.at<double>(i, j) = value2;
				continue;
			}
			//左右有值
			if (left >= 0 && right <= sceneWidth - 1)
			{
				ratio1 = double(j - left) / double(right - left);
				value1 = double(mappedLat.at<double>(i, left)) +
					double(mappedLat.at<double>(i, right) - mappedLat.at<double>(i, left)) * ratio1;
				mappedLat.at<double>(i, j) = value1;


				ratio1 = double(j - left) / double(right - left);
				value1 = double(mappedLon.at<double>(i, left)) +
					double(mappedLon.at<double>(i, right) - mappedLon.at<double>(i, left)) * ratio1;
				mappedLon.at<double>(i, j) = value1;
				continue;
			}
			//上边有值
			if (up >= 0)
			{
				mappedLat.at<double>(i, j) = mappedLat.at<double>(up, j);

				mappedLon.at<double>(i, j) = mappedLon.at<double>(up, j);
				continue;
			}
			//下边有值
			if (down <= sceneHeight - 1)
			{
				mappedLat.at<double>(i, j) = mappedLat.at<double>(down, j);

				mappedLon.at<double>(i, j) = mappedLon.at<double>(down, j);
				continue;
			}
			//左边有值
			if (left >= 0)
			{
				mappedLat.at<double>(i, j) = mappedLat.at<double>(i, left);

				mappedLon.at<double>(i, j) = mappedLon.at<double>(i, left);
				continue;
			}
			//右边有值
			if (right <= sceneWidth - 1)
			{
				mappedLat.at<double>(i, j) = mappedLat.at<double>(i, right);

				mappedLon.at<double>(i, j) = mappedLon.at<double>(i, right);
				continue;
			}
			//上下左右都没有值
			mappedLat.at<double>(i, j) = 0;
			mappedLon.at<double>(i, j) = 0;
		}
	}
	//投影纬度插值
	cv::GaussianBlur(mappedLat, mappedLat, cv::Size(5, 5), 1, 1);
	cv::GaussianBlur(mappedLon, mappedLon, cv::Size(5, 5), 1, 1);
	cv::GaussianBlur(DEM_out, mappedDEM, cv::Size(5, 5), 1, 1);
	return 0;
}

int Deflat::SLC_deramp(ComplexMat& slc, Mat& mappedDEM, Mat& mappedLat, Mat& mappedLon, const char* slcH5File, int mode)
{
	if (mappedDEM.rows != mappedLat.rows ||
		mappedDEM.rows != mappedLon.rows ||
		mappedDEM.cols != mappedLat.cols ||
		mappedDEM.cols != mappedLon.cols ||
		mappedDEM.type() != CV_16S ||
		mappedLat.type() != CV_64F ||
		mappedLon.type() != CV_64F ||
		mappedDEM.empty() ||
		!slcH5File
		)
	{
		fprintf(stderr, "SLC_deramp(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	int ret;
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, prf, start, end;
	int sceneHeight, sceneWidth, offset_row, offset_col;
	Mat lon_coef, lat_coef, statevec;
	string start_time, end_time;
	ret = conversion.read_int_from_h5(slcH5File, "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File, "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	//ret = conversion.read_int_from_h5(slcH5File, "offset_row", &offset_row);
	//if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	//ret = conversion.read_int_from_h5(slcH5File, "offset_col", &offset_col);
	//if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	//ret = conversion.read_array_from_h5(slcH5File, "lon_coefficient", lon_coef);
	//if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	//ret = conversion.read_array_from_h5(slcH5File, "lat_coefficient", lat_coef);
	//if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File, "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	//ret = conversion.read_double_from_h5(slcH5File, "range_spacing", &rangeSpacing);
	//if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	//ret = conversion.read_double_from_h5(slcH5File, "slant_range_first_pixel", &nearRangeTime);
	//if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	//nearRangeTime = 2.0 * nearRangeTime / VEL_C;
	ret = conversion.read_str_from_h5(slcH5File, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(slcH5File, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(slcH5File, "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_slc_from_h5(slcH5File, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
	Mat sate1 = Mat::zeros(sceneHeight, 3, CV_64F);
	double delta_t = statevec.at<double>(1, 0) - statevec.at<double>(0, 0);
	orbitStateVectors stateVectors(statevec, start, end, delta_t);
	stateVectors.applyOrbit();

	double dopplerFrequency = 0.0;

	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<double>(0, 0);
	lon = mappedLon.at<double>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

	int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	int numIterations = 0; Position pos; Velocity vel;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors.getPosition(midTime, pos);
		stateVectors.getVelocity(midTime, vel);
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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
	}
	double constant = (mode == 1 ? 4.0 * PI : 2.0 * PI);

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<double>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<double>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			r = -r / wavelength * constant;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	return 0;
}

int Deflat::slantrange_compute(Mat& slant_range, Mat& sate_pos,
	Mat& sate_vel, Mat& mappedDEM, Mat& mappedLat, Mat& mappedLon, const char* slcH5File)
{
	if (mappedDEM.rows != mappedLat.rows ||
		mappedDEM.rows != mappedLon.rows ||
		mappedDEM.cols != mappedLat.cols ||
		mappedDEM.cols != mappedLon.cols ||
		mappedDEM.type() != CV_16S ||
		mappedLat.type() != CV_64F ||
		mappedLon.type() != CV_64F ||
		mappedDEM.empty() ||
		!slcH5File
		)
	{
		fprintf(stderr, "SLC_deramp(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	int ret;
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, prf, start, end;
	int sceneHeight, sceneWidth, offset_row, offset_col;
	Mat lon_coef, lat_coef, statevec;
	string start_time, end_time;
	ret = conversion.read_int_from_h5(slcH5File, "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File, "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File, "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	ret = conversion.read_str_from_h5(slcH5File, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(slcH5File, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(slcH5File, "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	slant_range.create(sceneHeight, sceneWidth, CV_64F);
	sate_pos.create(sceneHeight, 3, CV_64F);
	sate_vel.create(sceneHeight, 3, CV_64F);
	double delta_t = statevec.at<double>(1, 0) - statevec.at<double>(0, 0);
	orbitStateVectors stateVectors(statevec, start, end, delta_t);
	stateVectors.applyOrbit();

	double dopplerFrequency = 0.0;

	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<double>(0, 0);
	lon = mappedLon.at<double>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

	int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	int numIterations = 0; Position pos; Velocity vel;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors.getPosition(midTime, pos);
		stateVectors.getVelocity(midTime, vel);
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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		stateVectors.getVelocity(time, vel);
		sate_pos.at<double>(i, 0) = pos.x;
		sate_pos.at<double>(i, 1) = pos.y;
		sate_pos.at<double>(i, 2) = pos.z;
		sate_vel.at<double>(i, 0) = vel.vx;
		sate_vel.at<double>(i, 1) = vel.vy;
		sate_vel.at<double>(i, 2) = vel.vz;
	}
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<double>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<double>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate_pos(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			slant_range.at<double>(i, j) = r;
		}
	}
	return 0;
}

int Deflat::SLCs_deramp(
	vector<string>& SLCH5Files,
	int reference,
	const char* demPath,
	vector<string>& outSLCH5Files
)
{
	if (SLCH5Files.size() < 1 ||
		reference < 1 || reference > SLCH5Files.size() ||
		!demPath ||
		outSLCH5Files.size() != SLCH5Files.size()
		)
	{
		fprintf(stderr, "SLCs_deramp(): input check failed!\n");
		return -1;
	}
	/*
	* 准备输出h5文件
	*/
	int ret;
	FormatConversion conversion;
	int images_num = outSLCH5Files.size();
	for (int i = 0; i < images_num; i++)
	{
		ret = conversion.creat_new_h5(outSLCH5Files[i].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	}
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, prf, start, end;
	int sceneHeight, sceneWidth, offset_row, offset_col;
	Mat lon_coef, lat_coef, dem, mappedDem, statevec, rangePos, azimuthPos;
	ComplexMat slc;
	string start_time, end_time, master_file;
	master_file = SLCH5Files[reference - 1];
	ret = conversion.read_int_from_h5(master_file.c_str(), "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(master_file.c_str(), "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(master_file.c_str(), "offset_row", &offset_row);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(master_file.c_str(), "offset_col", &offset_col);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(master_file.c_str(), "lon_coefficient", lon_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(master_file.c_str(), "lat_coefficient", lat_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(master_file.c_str(), "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(master_file.c_str(), "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	ret = conversion.read_double_from_h5(master_file.c_str(), "range_spacing", &rangeSpacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(master_file.c_str(), "slant_range_first_pixel", &nearRangeTime);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	nearRangeTime = 2.0 * nearRangeTime / VEL_C;
	ret = conversion.read_str_from_h5(master_file.c_str(), "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(master_file.c_str(), "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(master_file.c_str(), "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = Utils::computeImageGeoBoundry(lat_coef, lon_coef, sceneHeight, sceneWidth, offset_row, offset_col,
		&lonMax, &latMax, &lonMin, &latMin);
	if (return_check(ret, "computeImageGeoBoundry()", error_head)) return -1;
	ret = Utils::getSRTMDEM(demPath, dem, &lon_upperleft, &lat_upperleft, lonMin, lonMax, latMin, latMax);
	if (return_check(ret, "getSRTMDEM()", error_head)) return -1;
	Mat mappedLon, mappedLat;
	ret = demMapping(dem, mappedDem, mappedLat, mappedLon, lon_upperleft, lat_upperleft, offset_row, offset_col, sceneHeight, sceneWidth,
		prf, rangeSpacing, wavelength, nearRangeTime, start, end, statevec, 20);
	if (return_check(ret, "demMapping()", error_head)) return -1;
	ret = conversion.write_array_to_h5(outSLCH5Files[reference - 1].c_str(), "mapped_lat", mappedLat);
	ret = conversion.write_array_to_h5(outSLCH5Files[reference - 1].c_str(), "mapped_lon", mappedLon);
	for (int i = 0; i < images_num; i++)
	{
		ret = SLC_deramp(slc, mappedDem, mappedLat, mappedLon, SLCH5Files[i].c_str());
		if (return_check(ret, "SLC_deramp()", error_head)) return -1;
		ret = conversion.write_slc_to_h5(outSLCH5Files[i].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(SLCH5Files[i].c_str(), outSLCH5Files[i].c_str());
		if (return_check(ret, "Copy_para_from_h5_2_h5()", error_head)) return -1;
		ret = conversion.read_int_from_h5(SLCH5Files[i].c_str(), "offset_row", &offset_row);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(outSLCH5Files[i].c_str(), "offset_row", offset_row);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.read_int_from_h5(SLCH5Files[i].c_str(), "offset_col", &offset_col);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(outSLCH5Files[i].c_str(), "offset_col", offset_col);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(outSLCH5Files[i].c_str(), "range_len", sceneWidth);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(outSLCH5Files[i].c_str(), "azimuth_len", sceneHeight);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
	}
	
	return 0;
}

int Deflat::topography_phase_simulation(
	Mat& mappedDEM, 
	Mat& topography_phase,
	Mat& inc_coef, 
	double B_effect, 
	double nearRangeTime,
	int offset_row,
	int offset_col,
	double wavelength,
	double rangeSpacing
)
{
	if (mappedDEM.empty() ||
		mappedDEM.type() != CV_16S ||
		wavelength < 0.0 ||
		inc_coef.type() != CV_64F ||
		inc_coef.rows != 1 ||
		inc_coef.cols != 11 ||
		nearRangeTime <= 0.0 ||
		wavelength < 0.0 ||
		rangeSpacing < 0.0
		)
	{
		fprintf(stderr, "topography_phase_simulation(): input check failed!\n");
		return -1;
	}
	int rows = mappedDEM.rows;
	int cols = mappedDEM.cols;
	topography_phase.create(rows, cols, CV_64F);
    double offset_inc, scale_inc, offset_y, scale_y, a0, a1, a2, a3, a4, a5;
	offset_inc = inc_coef.at<double>(0, 0);
	scale_inc = inc_coef.at<double>(0, 1);
	offset_y = inc_coef.at<double>(0, 2);
	scale_y = inc_coef.at<double>(0, 3);
	a0 = inc_coef.at<double>(0, 4);
	a1 = inc_coef.at<double>(0, 5);
	a2 = inc_coef.at<double>(0, 6);
	a3 = inc_coef.at<double>(0, 7);
	a4 = inc_coef.at<double>(0, 8);
	a5 = inc_coef.at<double>(0, 9);

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		double r1, inc, jj;
		for (int j = 0; j < cols; j++)
		{
			jj = (double)j + (double)offset_col;
			jj = (jj - offset_y) / scale_y;
			inc = a0 + a1 * jj + a2 * jj * jj + a3 * jj * jj * jj + a4 * jj * jj * jj * jj + a5 * jj * jj * jj * jj * jj;
			inc = inc * scale_inc + offset_inc;

			r1 = nearRangeTime * VEL_C / 2.0 + rangeSpacing * (double)j;//主星斜距
			topography_phase.at<double>(i, j) = - 4 * PI * mappedDEM.at<short>(i, j) * B_effect / wavelength / r1 / sin(inc / 180.0 * PI);
		}
	}
	return 0;
}

int Deflat::computeImageGeoBoundry(
	Mat& lat_coefficient,
	Mat& lon_coefficient,
	int sceneHeight,
	int sceneWidth,
	int offset_row,
	int offset_col,
	double* lonMax,
	double* latMax,
	double* lonMin, 
	double* latMin
)
{
	if ( lon_coefficient.rows != 1 || 
		lon_coefficient.cols != 32 || 
		lon_coefficient.type() != CV_64F ||
		lat_coefficient.rows != 1 ||
		lat_coefficient.cols != 32 || 
		lat_coefficient.type() != CV_64F ||
		sceneHeight < 1 ||
		sceneWidth < 1 ||
		!lonMax || !lonMin || !latMax || !latMin
		)
	{
		fprintf(stderr, "computeImageGeoBoundry(): input check failed! \n");
		return -1;
	}
	int ret;
	Utils util;
	/*
	* 图像坐标转经纬坐标
	*/
	Mat row, col;
	row.create(4, 1, CV_64F); col.create(4, 1, CV_64F);
	row.at<double>(0, 0) = offset_row;//左上角
	col.at<double>(0, 0) = offset_col;

	row.at<double>(1, 0) = offset_row;//右上角
	col.at<double>(1, 0) = offset_col + sceneWidth;

	row.at<double>(2, 0) = offset_row + sceneHeight;//左下角
	col.at<double>(2, 0) = offset_col;

	row.at<double>(3, 0) = offset_row + sceneHeight;//右下角
	col.at<double>(3, 0) = offset_col + sceneWidth;
	Mat lon, lat;
	ret = util.coord_conversion(lon_coefficient, row, col, lon);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	ret = util.coord_conversion(lat_coefficient, row, col, lat);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	cv::minMaxLoc(lon, lonMin, lonMax);
	cv::minMaxLoc(lat, latMin, latMax);
	double extra = 5.0 / 6000;
	*lonMin = *lonMin - extra * 20;
	*lonMax = *lonMax + extra * 20;
	*latMin = *latMin - extra * 20;
	*latMax = *latMax + extra * 20;

	return 0;
}

int Deflat::getSRTMFileName(double lonMin, double lonMax, double latMin, double latMax, vector<string>& name)
{
	if (fabs(lonMin) > 180.0 ||
		fabs(lonMax) > 180.0 ||
		fabs(latMin) >= 60.0 ||
		fabs(latMax) >= 60.0
		)
	{
		fprintf(stderr, "getSRTMFileName(): input check failed!\n");
		return -1;
	}
	name.clear();
	char tmp[512];
	int maxRows = 24; int maxCols = 72; int startRow, endRow, startCol, endCol;
	double spacing = 5.0;
	startRow = (int)((60.0 - latMax) / spacing) + 1;
	endRow = (int)((60.0 - latMin) / spacing) + 1;
	startCol = (int)((lonMin + 180.0) / spacing) + 1;
	endCol = (int)((lonMax + 180.0) / spacing) + 1;
	if (startRow == endRow)
	{
		if (startCol == endCol)
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));
		}
		else
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, startRow);
			}
			name.push_back(string(tmp));
		}
	}
	else
	{
		if (startCol == endCol)
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));

			memset(tmp, 0, 512);
			if (startCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, endRow);
			}
			name.push_back(string(tmp));
		}
		else
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, endRow);
			}
			else if (endCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, endRow);
			}
			else if (endCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, endRow);
			}
			name.push_back(string(tmp));

			memset(tmp, 0, 512);
			if (startCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, endRow);
			}
			name.push_back(string(tmp));
		}
	}

	return 0;
}

int Deflat::downloadSRTM(const char* name)
{
	bool isConnect;
	DWORD dw;
	isConnect = IsNetworkAlive(&dw);
	if (!isConnect)
	{
		fprintf(stderr, "downloadSRTM(): network is not connected!\n");
		return -1;
	}
	int ret;
	string url = this->SRTMURL + name;
	string savefile = this->DEMPath + string("\\") + name;
	std::replace(savefile.begin(), savefile.end(), '/', '\\');
	HRESULT Result = URLDownloadToFileA(NULL, url.c_str(), savefile.c_str(), 0, NULL);
	if (Result != S_OK)
	{
		fprintf(stderr, "downloadSRTM(): download failded!\n");
		return -1;
	}
	return 0;
}

int Deflat::getSRTMDEM(
	const char* filepath,
	Mat& DEM_out,
	double* lonUL,
	double* latUL,
	double lonMin, 
	double lonMax, 
	double latMin,
	double latMax
)
{
	if (!filepath || !lonUL || !latUL) return -1;
	this->DEMPath = filepath;
	if (GetFileAttributesA(filepath) == -1)
	{
		if (_mkdir(filepath) != 0) return -1;
	}
	vector<string> srtmFileName;
	vector<bool> bAlreadyExist;
	int ret = getSRTMFileName(lonMin, lonMax, latMin, latMax, srtmFileName);
	if (return_check(ret, "getSRTMFileName()", error_head)) return -1;
	//判断文件是否已经存在
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string tmp = this->DEMPath + "\\" + srtmFileName[i];
		std::replace(tmp.begin(), tmp.end(), '/', '\\');
		if (-1 != GetFileAttributesA(tmp.c_str()))bAlreadyExist.push_back(true);
		else bAlreadyExist.push_back(false);
	}
	//不存在则下载
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		if (!bAlreadyExist[i])
		{
			ret = downloadSRTM(srtmFileName[i].c_str());
			if (return_check(ret, "downloadSRTM()", error_head)) return -1;
		}
	}
	//解压文件
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string folderName = srtmFileName[i];
		folderName = folderName.substr(0, folderName.length() - 4);
		string path = this->DEMPath + string("\\") + folderName;
		std::replace(path.begin(), path.end(), '/', '\\');
		if (-1 != GetFileAttributesA(path.c_str())) continue;
		string srcFile = this->DEMPath + "\\" + srtmFileName[i];
		std::replace(srcFile.begin(), srcFile.end(), '/', '\\');
		if (GetFileAttributesA(srcFile.c_str()) == -1) continue;
		ret = DigitalElevationModel::unzip(srcFile.c_str(), path.c_str());
		if (return_check(ret, "unzip()", error_head)) return -1;
	}


	int startRow, startCol, endRow, endCol;
	double lonUpperLeft, lonLowerRight, latUpperLeft, latLowerRight;
	int total_rows, total_cols;

	//DEM在一个SRTM方格内
	if (srtmFileName.size() == 1)
	{
		total_rows = 6000, total_cols = 6000;
		int xx, yy;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		latUpperLeft = 60.0 - (yy - 1) * 5.0;
		latLowerRight = latUpperLeft - 5.0;
		lonUpperLeft = -180.0 + (xx - 1) * 5.0;
		lonLowerRight = lonUpperLeft + 5.0;

		startRow = (latUpperLeft - latMax) / this->latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / this->latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
		startCol = startCol < 1 ? 1 : startCol;
		startCol = startCol > total_cols ? total_cols : startCol;
		endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
		endCol = endCol < 1 ? 1 : endCol;
		endCol = endCol > total_cols ? total_cols : endCol;

		string folderName = srtmFileName[0];
		folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
		string path = this->DEMPath + string("\\") + folderName;
		path = path + string("\\") + folderName + string(".tif");
		Mat outDEM = Mat::zeros(6000, 6000, CV_16S);
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
		//if (return_check(ret, "geotiffread()", error_head)) return -1;
		outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
		*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
		*latUL = latUpperLeft - (startRow - 1) * latSpacing;
	}
	//DEM在2个方格内
	else if (srtmFileName.size() == 2)
	{
		int xx, yy, xx2, yy2;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		sscanf(srtmFileName[1].c_str(), "srtm_%d_%d.zip", &xx2, &yy2);
		//同一列
		if (xx == xx2)
		{
			total_rows = 6000 * 2; total_cols = 6000;
			latUpperLeft = 60.0 - ((yy < yy2 ? yy : yy2) - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = -180.0 + (xx - 1) * 5.0;
			lonLowerRight = lonUpperLeft + 5.0;

			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;


			Mat outDEM, outDEM2;

			if (yy < yy2)
			{
				string folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				string path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[1];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM2 = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;
				cv::vconcat(outDEM, outDEM2, outDEM);
			}
			else
			{
				string folderName = srtmFileName[1];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				string path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM2 = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;
				cv::vconcat(outDEM, outDEM2, outDEM);
			}

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		}
		//同一行
		else if (yy == yy2)
		{
			total_cols = 6000 * 2; total_rows = 6000;
			//跨越-180.0/180.0线
			if ((xx == 1 && xx2 == 72) || (xx == 72 && xx2 == 1))
			{
				latUpperLeft = 60.0 - (yy - 1) * 5.0;
				latLowerRight = latUpperLeft - 5.0;
				lonUpperLeft = 175.0;
				lonLowerRight = -175.0;
				startRow = (latUpperLeft - latMax) / this->latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / this->latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMax - lonUpperLeft) / this->lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMin - lonUpperLeft + 360.0) / this->lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;

				Mat outDEM, outDEM2;

				if (xx > xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}
				else
				{
					string folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}

				outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
				*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
				*latUL = latUpperLeft - (startRow - 1) * latSpacing;
			}
			else
			{
				latUpperLeft = 60.0 - (yy - 1) * 5.0;
				latLowerRight = latUpperLeft - 5.0;
				lonUpperLeft = -180.0 + ((xx < xx2 ? xx : xx2) - 1) * 5.0;
				lonLowerRight = lonUpperLeft + 10.0;

				startRow = (latUpperLeft - latMax) / this->latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / this->latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;


				Mat outDEM, outDEM2;

				if (xx < xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}
				else
				{
					string folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}

				outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
				*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
				*latUL = latUpperLeft - (startRow - 1) * latSpacing;
			}
		}
		else
		{
			return -1;
		}



	}
	//DEM在4个方格内
	else if (srtmFileName.size() == 4)
	{
		int xx, yy, xx2, yy2, xx3, yy3, xx4, yy4, temp;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		sscanf(srtmFileName[1].c_str(), "srtm_%d_%d.zip", &xx2, &yy2);
		sscanf(srtmFileName[2].c_str(), "srtm_%d_%d.zip", &xx3, &yy3);
		sscanf(srtmFileName[3].c_str(), "srtm_%d_%d.zip", &xx4, &yy4);
		total_rows = 6000 * 2; total_cols = 6000 * 2;
		//跨越-180.0/180.0线
		if (lonMax * lonMin < 0 && (fabs(lonMin) + fabs(lonMax)) > 180.0)
		{
			startRow = (int)((60.0 - latMax) / 5.0) + 1;
			endRow = (int)((60.0 - latMin) / 5.0) + 1;
			endCol = (int)((lonMin + 180.0) / 5.0) + 1;
			startCol = (int)((lonMax + 180.0) / 5.0) + 1;
			latUpperLeft = 60.0 - (startRow - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = 175.0;
			lonLowerRight = -175.0;



			Mat outDEM, outDEM2, outDEM3;

			char tmpstr[512];
			const char* format = NULL;
			if (startCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, startRow);
			string folderName(tmpstr);
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			string path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, startRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM, outDEM2, outDEM);

			if (startCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);


			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin - lonUpperLeft + 360.0) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		}
		else
		{
			startRow = (int)((60.0 - latMax) / 5.0) + 1;
			endRow = (int)((60.0 - latMin) / 5.0) + 1;
			startCol = (int)((lonMin + 180.0) / 5.0) + 1;
			endCol = (int)((lonMax + 180.0) / 5.0) + 1;
			latUpperLeft = 60.0 - (startRow - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = -180.0 + (startCol - 1) * 5.0;
			lonLowerRight = lonUpperLeft + 10.0;



			Mat outDEM, outDEM2, outDEM3;

			char tmpstr[512];
			const char* format = NULL;
			if (startCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, startRow);
			string folderName(tmpstr);
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			string path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, startRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM, outDEM2, outDEM);

			if (startCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);

			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		}

	}
	else return -1;
	this->rows = DEM_out.rows;
	this->cols = DEM_out.cols;
	return 0;
}
