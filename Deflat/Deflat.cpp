// Deflat.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include"..\include\Deflat.h"
#include"..\include\FormatConversion.h"
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

int Deflat::deflat(const Mat& stateVec1, const Mat& stateVec2, const Mat& lon_coef, const Mat& lat_coef, const Mat& phase, int offset_row, int offset_col, double height, double time_interval, double time_interval2, int mode, double wave_length, Mat& phase_deflated, Mat& flat_phase)
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
		offset_row < 0 ||
		offset_col < 0 ||
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
	lat.copyTo(flat_phase);
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
		offset_row < 0 || offset_col < 0 || interp_interval1 < 0.0 || interp_interval2 < 0.0 ||
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
