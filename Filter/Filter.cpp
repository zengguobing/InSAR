// Filter.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include"..\include\Filter.h"
#include<tchar.h>
#include <atlconv.h>
#ifdef _DEBUG
#pragma comment(lib,"ComplexMat_d.lib")
#pragma comment(lib, "Utils_d.lib")
#else
#pragma comment(lib,"ComplexMat.lib")
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

Filter::Filter()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "FILTER_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "FILTER_DLL_ERROR: error happens when using parallel computing in function: ");
}

Filter::~Filter()
{
}

int Filter::czt2(Mat& src, Mat& dst, int M, int N, double theta0, double phi0)
{
	if (src.cols < 1 ||
		src.rows < 1 ||
		M < 1 ||
		N < 1 ||
		src.type() != CV_64FC2)
	{
		fprintf(stderr, "czt2(): input check failed!\n\n");
		return -1;
	}
	theta0 = -theta0;
	int nc = src.cols;
	int L, i;
	Mat g, h, y;

	i = 1;
	do
	{
		i = i * 2;
	} while (i < M + N);
	L = i;

	Mat planes[] = { Mat::zeros(L, nc, CV_64F), Mat::zeros(L, nc, CV_64F) };
	merge(planes, 2, h);
	merge(planes, 2, g);
	merge(planes, 2, y);


	int j;
	//////////////////////////////////////////h(n)赋值
	for (j = 0; j < nc; j++)
	{
		for (i = 0; i < M; i++)
		{
			h.at<Vec2d>(i, j)[0] = cos(phi0 * 0.5 * i * i);//实部
			h.at<Vec2d>(i, j)[1] = sin(phi0 * 0.5 * i * i);//虚部
		}
		for (i = M; i < L - N + 1; i++)
		{
			h.at<Vec2d>(i, j)[0] = 0;//实部
			h.at<Vec2d>(i, j)[1] = 0;//虚部
		}
		for (i = L - N + 1; i < L; i++)
		{
			h.at<Vec2d>(i, j)[0] = cos(phi0 * 0.5 * ((double)(L - i)) * ((double)(L - i)));//实部
			h.at<Vec2d>(i, j)[1] = sin(phi0 * 0.5 * ((double)(L - i)) * ((double)(L - i)));//虚部
		}
	}


	////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////x(n)加权并补零得到g(n)
	Mat W(N, nc, CV_64FC2, Scalar::all(0));
	//for (j = 0; j < nc; j++)
	//{
	//	for (i = 0; i < N; i++)
	//	{

	//		W.at<Vec2d>(i, j)[0] = cos(theta0 * i - phi0 * 0.5 * i * i);/*实部*/
	//		W.at<Vec2d>(i, j)[1] = sin(theta0 * i - phi0 * 0.5 * i * i);/*虚部*/

	//	}
	//}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < nc; j++)
		{
			W.at<Vec2d>(i, j)[0] = cos(theta0 * i - phi0 * 0.5 * i * i);/*实部*/
			W.at<Vec2d>(i, j)[1] = sin(theta0 * i - phi0 * 0.5 * i * i);/*虚部*/
		}
	}
	Mat tmp;
	mulSpectrums(W, src, tmp, 0, false);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < nc; j++)
		{
			g.at<Vec2d>(i, j)[0] = tmp.at<Vec2d>(i, j)[0];
			g.at<Vec2d>(i, j)[1] = tmp.at<Vec2d>(i, j)[1];
		}
	}
	//mulSpectrums(W, src, g(Range(0, N), Range(0, nc)), 0, false);
	for (j = 0; j < nc; j++)
	{
		for (i = N; i < L; i++)
		{
			g.at<Vec2d>(i, j)[0] = 0;
			g.at<Vec2d>(i, j)[1] = 0;
		}
	}

	/////////////////////////////////////////////////////

	/*h(n)和g(n)进行傅里叶变换
	得到H(k)和G(k)
	*/
	transpose(h, h);
	transpose(g, g);
	dft(h, h, DFT_ROWS);
	dft(g, g, DFT_ROWS);
	transpose(h, h);
	transpose(g, g);
	mulSpectrums(g, h, y, 0, false);
	transpose(y, y);
	idft(y, y, DFT_ROWS);
	transpose(y, y);
	Mat result;
	Mat planes1[] = { Mat::zeros(M, nc, CV_64F), Mat::zeros(M, nc, CV_64F) };
	merge(planes1, 2, result);
	for (j = 0; j < nc; j++)
	{
		for (i = 0; i < M; i++)
		{
			result.at<Vec2d>(i, j)[0] = cos(phi0 * 0.5 * i * i);
			result.at<Vec2d>(i, j)[1] = -sin(phi0 * 0.5 * i * i);
		}
	}

	mulSpectrums(result, y(Range(0, M), Range(0, nc)), result, 0, false);
	dst =  result / double(L);
	return 0;
}

int Filter::meanfilter(Mat& Src, int WndSize)
{
	if (WndSize % 2 == 0 ||
		WndSize < 1 ||
		Src.rows < 1 ||
		Src.cols < 1 ||
		Src.type() != CV_64FC2)
	{
		fprintf(stderr, "meanfilter(): input check failed!\n\n");
		return -1;
	}
	Mat planes[] = { Mat::zeros(Src.rows, Src.cols, CV_64F), Mat::zeros(Src.rows, Src.cols, CV_64F) };
	split(Src, planes);
	blur(planes[0], planes[0], Size(WndSize, WndSize));
	blur(planes[1], planes[1], Size(WndSize, WndSize));
	merge(planes, 2, Src);
	return 0;
}

int Filter::fftshift2(Mat& matrix)
{
	if (matrix.rows < 2 ||
		matrix.cols < 2 ||
		matrix.channels() != 1)
	{
		fprintf(stderr, "fftshift2(): input check failed!\n\n");
		return -1;
	}
	matrix = matrix(Rect(0, 0, matrix.cols & -2, matrix.rows & -2));
	int cx = matrix.cols / 2;
	int cy = matrix.rows / 2;
	Mat tmp;
	Mat q0(matrix, Rect(0, 0, cx, cy));
	Mat q1(matrix, Rect(cx, 0, cx, cy));
	Mat q2(matrix, Rect(0, cy, cx, cy));
	Mat q3(matrix, Rect(cx, cy, cx, cy));

	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);
	return 0;
}

int Filter::slope_adaptive_filter(Mat& phase, Mat& phase_filter, int wndsize_filter, int wndsize_prefilter)
{
	if (wndsize_filter % 2 == 0 ||
		wndsize_prefilter % 2 == 0 ||
		wndsize_filter < 3 ||
		wndsize_prefilter < 3 ||
		wndsize_filter > int(phase.rows / 2) ||
		wndsize_filter > int(phase.cols / 2) ||
		wndsize_prefilter > int(phase.rows / 2) ||
		wndsize_prefilter > int(phase.cols / 2) ||
		phase.type() != CV_64F ||
		phase.channels() > 1)
	{
		fprintf(stderr, "slope_adaptive_filter(): input check failed!\n\n");
		return -1;
	}
	int nn, mm;
	nn = getOptimalDFTSize(2 * wndsize_filter);
	mm = nn;
	double pi = 3.1415926535;
	int Radius = (wndsize_filter - 1) / 2; /*窗半径*/
	int nr_orig = phase.rows;/*原始尺寸rows*/
	int nc_orig = phase.cols;/*原始尺寸cols*/
	Mat phase_enlarged(nr_orig + wndsize_filter - 1, nc_orig + wndsize_filter - 1, CV_64F, Scalar::all(0));/*扩充矩阵*/
	int nr_new = phase_enlarged.rows;
	int nc_new = phase_enlarged.cols;

	phase.copyTo(phase_enlarged(Range(Radius, nr_new - Radius), Range(Radius, nc_new - Radius)));/*原始矩阵复制到中心*/

	/*边缘赋值*/
	Mat temp;
	flip(phase_enlarged(Range(Radius, 2 * Radius), Range(0, nc_new)), temp, 0);/*up - down翻转*/

	temp.copyTo(phase_enlarged(Range(0, Radius), Range(0, nc_new)));

	flip(phase_enlarged(Range(0, nr_new), Range(Radius, 2 * Radius)), temp, 1);/*左右反转*/

	temp.copyTo(phase_enlarged(Range(0, nr_new), Range(0, Radius)));

	flip(phase_enlarged(Range(nr_new - 2 * Radius, nr_new - Radius), Range(0, nc_new)), temp, 0);/*上下翻转*/

	temp.copyTo(phase_enlarged(Range(nr_new - Radius, nr_new), Range(0, nc_new)));

	flip(phase_enlarged(Range(0, nr_new), Range(nc_new - 2 * Radius, nc_new - Radius)), temp, 1);/*左右翻转*/

	temp.copyTo(phase_enlarged(Range(0, nr_new), Range(nc_new - Radius, nc_new)));

	Mat phase_update(nr_new, nc_new, CV_64FC2, Scalar::all(0));/*相位矩阵转换为复数*/

	int i, j;
	for (i = 0; i < nr_new; i++)
	{
		for (j = 0; j < nc_new; j++)
		{
			phase_update.at<Vec2d>(i, j)[0] = cos(phase_enlarged.at<double>(i, j));/*实部*/
			phase_update.at<Vec2d>(i, j)[1] = sin(phase_enlarged.at<double>(i, j));/*虚部*/

		}
	}


	//int nn = 32;/*2D FFT 点数*/
	//int mm = 32;/*CZT 点数*/
	int mn = mm * nn;/*频谱总分辨率*/

	Mat phase_filtered = Mat::zeros(nr_new, nc_new, CV_64F);


	Mat tempi(wndsize_filter, 1, CV_64F, Scalar::all(0));/*行线性相位网格*/
	Mat tempj(1, wndsize_filter, CV_64F, Scalar::all(0));/*列线性相位网格*/

	for (i = 0; i < wndsize_filter; i++)
	{
		tempi.at<double>(i, 0) = (double)i;
		tempj.at<double>(0, i) = (double)i;
	}

	double phi0 = 2.0 * pi / ((double)mn); /*CZT变换参数*/
	int ret;
	volatile bool parallel_flag = true;
//#pragma omp parallel for schedule(guided) \
//	private(ret)
	for (int i = Radius; i < nr_new - Radius; i++)
	{
#pragma omp parallel for schedule(guided) \
		private(ret)
		//if (!parallel_flag) continue;
		for (int j = Radius; j < nc_new - Radius; j++)
		{
			if (!parallel_flag) continue;
			int k, kk;
			Mat phase_estimation(wndsize_filter, wndsize_filter, CV_64FC2, Scalar::all(0));/*滤波窗口内相位*/
			Mat window_mean;
			Mat planes[] = { Mat::zeros(wndsize_filter, wndsize_filter, CV_64F),
				Mat::zeros(wndsize_filter, wndsize_filter, CV_64F) };
			Point peak_loc;
			double fi, fj, fii, fjj;
			double theta0_i, theta0_j;
			Mat phase_czt1, phase_czt2, AA, phase_0_matrix;
			Mat aa(wndsize_filter, wndsize_filter, CV_64FC2, Scalar::all(0));
			Mat phase_0(1, 1, CV_64FC2, Scalar::all(0));
			phase_update(Range(i - Radius, i + Radius + 1), Range(j - Radius, j + Radius + 1)).copyTo(phase_estimation);
			phase_estimation.copyTo(window_mean);
			/*频率估计前预滤波*/
			ret = meanfilter(window_mean, wndsize_prefilter);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			
			//if (parallel_flag_change(parallel_flag, ret)) continue;

			/*nn点傅里叶变换*/
			copyMakeBorder(window_mean, window_mean, 0, nn - wndsize_filter, 0, nn - wndsize_filter, BORDER_CONSTANT, Scalar::all(0));
			dft(window_mean, window_mean);
			split(window_mean, planes);
			magnitude(planes[0], planes[1], planes[0]);
			ret = fftshift2(planes[0]);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			if (ret < 0) continue;
			if (parallel_flag_change(parallel_flag, ret)) continue;
			minMaxLoc(planes[0], NULL, NULL, NULL, &peak_loc);
			fi = ((double)(peak_loc.y - nn / 2 - 1)) / ((double)nn);/*此处减2是为了扩大CZT变换的搜索范围*/
			fj = ((double)(peak_loc.x - nn / 2 - 1)) / ((double)nn);

			theta0_i = 2.0 * pi * fi;
			theta0_j = 2.0 * pi * fj;

			ret = czt2(phase_estimation, phase_czt1, 3 * mm, phase_estimation.rows, theta0_i, phi0);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			//if (parallel_flag_change(parallel_flag, ret)) continue;
			transpose(phase_czt1, phase_czt1);
			ret = czt2(phase_czt1, phase_czt2, 3 * mm, phase_czt1.rows, theta0_j, phi0);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			//if (parallel_flag_change(parallel_flag, ret)) continue;
			transpose(phase_czt2, phase_czt2);
			split(phase_czt2, planes);
			magnitude(planes[0], planes[1], planes[0]);
			minMaxLoc(planes[0], NULL, NULL, NULL, &peak_loc);

			fii = fi + ((double)(peak_loc.y) / (double)mn);
			fjj = fj + ((double)(peak_loc.x) / (double)mn);/*频谱细化后峰值位置*/

			theta0_i = 2.0 * pi * fii;
			theta0_j = 2.0 * pi * fjj;/*更新细化值*/

			/*tempi = tempi * theta0_i;
			tempj = tempj * theta0_j;*/
			Mat one(tempi.rows, tempi.cols, CV_64F, Scalar::all(1));
			Mat one_t(tempj.rows, tempj.cols, CV_64F, Scalar::all(1));

			AA = ((tempi * theta0_i) * one_t) + (one * (tempj * theta0_j));
			for (k = 0; k < wndsize_filter; k++)
			{
				for (kk = 0; kk < wndsize_filter; kk++)
				{
					aa.at<Vec2d>(k, kk)[0] = cos(AA.at<double>(k, kk));/*实部*/
					aa.at<Vec2d>(k, kk)[1] = sin(AA.at<double>(k, kk));/*虚部*/
				}
			}
			//#pragma omp critical 
			//{
			mulSpectrums(phase_update(Range(i - Radius, i + Radius + 1), Range(j - Radius, j + Radius + 1)), aa, phase_0_matrix, 0, true);/*异常*/
			//}

			phase_0.at<Vec2d>(0, 0)[0] = mean(phase_0_matrix)[0];
			phase_0.at<Vec2d>(0, 0)[1] = mean(phase_0_matrix)[1];

			mulSpectrums(phase_0, aa(Range((wndsize_filter - 3) / 2, (wndsize_filter - 1) / 2),
				Range((wndsize_filter - 3) / 2, (wndsize_filter - 1) / 2)), phase_0, 0, false);

			/*转换为相位*/
			phase_filtered.at<double>(i, j) = atan2(phase_0.at<Vec2d>(0, 0)[1], phase_0.at<Vec2d>(0, 0)[0]);
		}

		fprintf(stdout, "process: %lf %\n", double(i - Radius) / double(nr_new - 2 * Radius + 1) * 100);
	}
	//if (parallel_check(parallel_flag, "slope_adaptive_filter()", parallel_error_head)) return -1;
	phase_filtered(Range(Radius, nr_new - Radius), Range(Radius, nc_new - Radius)).copyTo(phase_filter);
	return 0;
}

int Filter::filter_dl(const char* filter_dl_path, const char* tmp_path, const char* dl_model_file, Mat& phase, Mat& phase_filtered)
{
	if (filter_dl_path == NULL ||
		dl_model_file == NULL ||
		tmp_path == NULL||
		phase.rows < 1 ||
		phase.cols < 1 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F)
	{
		fprintf(stderr, "filter_dl(): input check failed!\n\n");
		return -1;
	}
	int nr = phase.rows;
	int nc = phase.cols;
	int ret;
	Utils util;
	Mat cos = Mat::zeros(nr, nc, CV_64F);
	Mat sin = Mat::zeros(nr, nc, CV_64F);
	phase_filtered = Mat::zeros(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			cos.at<double>(i, j) = std::cos(phase.at<double>(i, j));
			sin.at<double>(i, j) = std::sin(phase.at<double>(i, j));
		}
	}
	std::string cos_file(tmp_path);
	std::replace(cos_file.begin(), cos_file.end(), '/', '\\');
	std::string sin_file(tmp_path);
	std::replace(sin_file.begin(), sin_file.end(), '/', '\\');
	cos_file.append("\\cos.dat");
	sin_file.append("\\sin.dat");
	ret = util.cvmat2bin(cos_file.c_str(), cos);
	if (return_check(ret, "util.cvmat2bin(*, *)", error_head)) return -1;
	ret = util.cvmat2bin(sin_file.c_str(), sin);
	if (return_check(ret, "util.cvmat2bin(*, *)", error_head)) return -1;
	///////////////////////////创建并调用深度学习滤波进程//////////////////////
	USES_CONVERSION;
	LPWSTR szCommandLine = new TCHAR[512];
	string Filter_dl_path(filter_dl_path);
	std::replace(Filter_dl_path.begin(), Filter_dl_path.end(), '/', '\\');
	wcscpy(szCommandLine, A2W(Filter_dl_path.c_str()));
	wcscat(szCommandLine, L"\\filter_dl.exe ");
	wcscat(szCommandLine, A2W(dl_model_file));
	wcscat(szCommandLine, L" ");
	wcscat(szCommandLine, A2W(cos_file.c_str()));
	wcscat(szCommandLine, L" ");
	wcscat(szCommandLine, A2W(sin_file.c_str()));
	STARTUPINFO si;
	PROCESS_INFORMATION p_i;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&p_i, sizeof(p_i));
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = FALSE;
	BOOL bRet = ::CreateProcess(
		NULL,           // 不在此指定可执行文件的文件名
		szCommandLine,      // 命令行参数
		NULL,           // 默认进程安全性
		NULL,           // 默认线程安全性
		FALSE,          // 指定当前进程内的句柄不可以被子进程继承
		CREATE_NEW_CONSOLE, // 为新进程创建一个新的控制台窗口
		NULL,           // 使用本进程的环境变量
		NULL,           // 使用本进程的驱动器和目录
		&si,
		&p_i);
	if (bRet)
	{
		char filter_job_name[512]; filter_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(filter_job_name, "FILTER_%lld", tt);
		string filter_job_name_string(filter_job_name);
		HANDLE hd = CreateJobObjectA(NULL, filter_job_name_string.c_str());
		if (hd)
		{
			JOBOBJECT_EXTENDED_LIMIT_INFORMATION extLimitInfo;
			extLimitInfo.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
			BOOL retval = SetInformationJobObject(hd, JobObjectExtendedLimitInformation, &extLimitInfo, sizeof(extLimitInfo));
			if (retval)
			{
				if (p_i.hProcess)
				{
					retval = AssignProcessToJobObject(hd, p_i.hProcess);
				}
			}
		}
		WaitForSingleObject(p_i.hProcess, INFINITE);
		if (szCommandLine != NULL) delete[] szCommandLine;
		::CloseHandle(p_i.hThread);
		::CloseHandle(p_i.hProcess);
	}
	else
	{
		fprintf(stderr, "filter_dl(): create filter_dl.exe process failed!\n\n");
		if (szCommandLine != NULL) delete[] szCommandLine;
		return -1;
	}
	///////////////////////////创建并调用深度学习滤波进程//////////////////////
	cos_file.append(".out");
	sin_file.append(".out");
	Mat cos1, sin1;
	ret = util.bin2cvmat(cos_file.c_str(), cos1);
	if (return_check(ret, "util.bin2cvmat(*, *)", error_head)) return -1;
	ret = util.bin2cvmat(sin_file.c_str(), sin1);
	if (return_check(ret, "util.bin2cvmat(*, *)", error_head)) return -1;
	cos = cos - cos1;
	sin = sin - sin1;
	ComplexMat tmp(cos, sin);
	tmp.GetPhase().copyTo(phase_filtered);
	return 0;
}

int Filter::Goldstein_filter(Mat& phase, Mat& phase_filter, double alpha, int n_win, int n_pad)
{
	if (phase.cols < 3 ||
		phase.rows < 3 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F ||
		alpha <= 0 ||
		n_win < 5 ||
		n_pad < 0 
		)
	{
		fprintf(stderr, "Goldstein_filter(): input check failed!\n\n");
		return -1;
	}
	int n_i = phase.rows;
	int n_j = phase.cols;
	ComplexMat ph;
	Mat cos, sin;
	int ret;
	Utils util;
	ret = util.phase2cos(phase, ph.re, ph.im);
	//if (return_check(ret, "util.phase2cos(*, *, *)", error_head)) return -1;
	//ph.SetIm(sin);
	//ph.SetRe(cos);

	ComplexMat ph_out(n_i, n_j);
	int n_inc = floor(n_win / 4);
	int n_win_i = ceil(n_i / n_inc) - 1;
	int n_win_j = ceil(n_j / n_inc) - 1;
	int x = floor(n_win / 2 - 1);
	Mat qua_wnd = Mat::zeros(x + 1, x + 1, CV_64F);
	for (int i = 0; i <= x; i++)
	{
		for (int j = 0; j <= x; j++)
		{
			qua_wnd.at<double>(i, j) = double(i + j);
		}
	}
	qua_wnd.at<double>(0, 0) = 1e-6;
	Mat fliped_wnd;
	flip(qua_wnd, fliped_wnd, 1);
	hconcat(qua_wnd, fliped_wnd, qua_wnd);
	flip(qua_wnd, fliped_wnd, 0);
	vconcat(qua_wnd, fliped_wnd, qua_wnd);
	Mat guasswin = Mat::zeros(7, 1, CV_64F);
	double val[] = { 0.0439369336234074, 0.249352208777296, 0.706648277857716, 1, 0.706648277857716, 0.249352208777296, 0.0439369336234074 };
	memcpy(guasswin.data, val, sizeof(double) * 7);
	Mat guasswin_t;
	transpose(guasswin, guasswin_t);
	guasswin = guasswin * guasswin_t;
	GenerateGaussMask(guasswin, 7, 7, 1.0);
	int n_win_ex = n_win + n_pad;
	ComplexMat ph_bit(n_win_ex, n_win_ex);
	ComplexMat temp, temp1, temp2, fft_out, ph_filt;
	Mat wf, wf2, wind_func, tmp, tmp1, tmp2, H;
	qua_wnd.copyTo(wind_func);
	double median = 1.0;
	int i1, i2, i_shift, j_shift, ix1, ix2, j1, j2;
	for (ix1 = 1; ix1 <= n_win_i; ix1++)
	{
		wind_func.copyTo(wf);
		i1 = (ix1 - 1) * n_inc + 1;
		i2 = i1 + n_win - 1;
		if (i2 > n_i)
		{
			i_shift = i2 - n_i;
			i2 = n_i;
			i1 = n_i - n_win + 1;
			tmp1 = Mat::zeros(i_shift, n_win, CV_64F);
			wf(Range(0, n_win - i_shift), Range(0, wf.cols)).copyTo(tmp2);
			vconcat(tmp1, tmp2, wf);
		}
		for (ix2 = 1; ix2 <= n_win_j; ix2++)
		{
			wf.copyTo(wf2);
			j1 = (ix2 - 1) * n_inc + 1;
			j2 = j1 + n_win - 1;
			if (j2 > n_j)
			{
				j_shift = j2 - n_j;
				j2 = n_j;
				j1 = n_j - n_win + 1;
				tmp1 = Mat::zeros(n_win, j_shift, CV_64F);
				wf2(Range(0, wf2.rows), Range(0, n_win - j_shift)).copyTo(tmp2);
				hconcat(tmp1, tmp2, wf2);
			}
			if (wf2.cols != n_win || wf2.rows != n_win)
			{
				fprintf(stderr, "Goldstein_filter(): wf2.size and n_win mismatch, please check to make sure n_win is even!\n\n");
				return -1;
			}
			temp = ph(cv::Range(i1 - 1, i2), cv::Range(j1 - 1, j2));
			ret = ph_bit.SetValue(cv::Range(0, n_win), cv::Range(0, n_win), temp);
			if (return_check(ret, "ComplexMat::SetValue(*, *, *)", error_head)) return -1;
			ret = util.fft2(ph_bit, fft_out);
			if (return_check(ret, "util.fft2(*, *)", error_head)) return -1;
			H = fft_out.GetMod();
			ret = fftshift2(H);
			if (return_check(ret, "fftshift2(*, *)", error_head)) return -1;
			GaussianFilter(H, H, guasswin);
			//filter2D(H, H, -1, guasswin, Point(-1, -1), 0.0, BORDER_CONSTANT);
			ret = util.ifftshift(H);
			if (return_check(ret, "ifftshift(*, *)", error_head)) return -1;
			H.copyTo(tmp);
			tmp = tmp.reshape(0, 1);
			cv::sort(tmp, tmp, SORT_ASCENDING + SORT_EVERY_ROW);
			if (tmp.cols % 2 == 1)
			{
				median = tmp.at<double>(0, int((tmp.cols + 1) / 2) - 1);
			}
			else
			{
				median = tmp.at<double>(0, int(tmp.cols / 2) - 1) + tmp.at<double>(0, int(tmp.cols / 2));
				median = median / 2.0;
			}
			if (fabs(median) > 1e-8)
			{
				H = H / median;
			}
			cv::pow(H, alpha, H);
			fft_out = fft_out * H;
			ret = util.ifft2(fft_out, temp);
			if (return_check(ret, "util.ifft2(*, *)", error_head)) return -1;
			Mat re, im;
			re = temp.GetRe();
			im = temp.GetIm();
			temp1 = temp(Range(0, n_win), Range(0, n_win));
			ph_filt = temp1 * wf2;
			temp = ph_out(Range(i1 - 1, i2), Range(j1 - 1, j2)) + ph_filt;
			ret = ph_out.SetValue(Range(i1 - 1, i2), Range(j1 - 1, j2), temp);
			if (return_check(ret, "ph_out.SetValue(*, *, *)", error_head)) return -1;
		}
		fprintf(stdout, "Goldstein filtering process: %d / %d\n", ix1, n_win_i);
	}
	ph_out.GetPhase().copyTo(phase_filter);
	return 0;
}

int Filter::Goldstein_filter_parallel(Mat& phase, Mat& phase_filter, double alpha, int n_win, int n_pad)
{
	if (phase.cols < 3 ||
		phase.rows < 3 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F ||
		alpha <= 0 ||
		n_win < 5 ||
		n_pad < 0
		)
	{
		fprintf(stderr, "Goldstein_filter_parallel(): input check failed!\n\n");
		return -1;
	}
	int n_i = phase.rows;
	int n_j = phase.cols;
	ComplexMat ph;
	Mat cos, sin;
	
	Utils util;
	util.phase2cos(phase, ph.re, ph.im);

	ComplexMat ph_out(n_i, n_j);
	int n_inc = floor(n_win / 4);
	int n_win_i = ceil(n_i / n_inc) - 1;
	int n_win_j = ceil(n_j / n_inc) - 1;
	int x = floor(n_win / 2 - 1);
	Mat qua_wnd = Mat::zeros(x + 1, x + 1, CV_64F);
	for (int i = 0; i <= x; i++)
	{
		for (int j = 0; j <= x; j++)
		{
			qua_wnd.at<double>(i, j) = double(i + j);
		}
	}
	qua_wnd.at<double>(0, 0) = 1e-6;
	Mat fliped_wnd, guasswin;
	flip(qua_wnd, fliped_wnd, 1);
	hconcat(qua_wnd, fliped_wnd, qua_wnd);
	flip(qua_wnd, fliped_wnd, 0);
	vconcat(qua_wnd, fliped_wnd, qua_wnd);
	GenerateGaussMask(guasswin, 7, 7, 1.0);
	int n_win_ex = n_win + n_pad;
	
	for (int ix1 = 1; ix1 <= n_win_i; ix1++)
	{
		Mat wf, wind_func, tmp1, tmp2;
		qua_wnd.copyTo(wind_func);
		int i1, i2, i_shift;
		wind_func.copyTo(wf);
		i1 = (ix1 - 1) * n_inc + 1;
		i2 = i1 + n_win - 1;
		if (i2 > n_i)
		{
			i_shift = i2 - n_i;
			i2 = n_i;
			i1 = n_i - n_win + 1;
			tmp1 = Mat::zeros(i_shift, n_win, CV_64F);
			wf(Range(0, n_win - i_shift), Range(0, wf.cols)).copyTo(tmp2);
			vconcat(tmp1, tmp2, wf);
		}
#pragma omp parallel for schedule(guided)
		for (int ix2 = 1; ix2 <= n_win_j; ix2++)
		{
			ComplexMat ph_bit(n_win_ex, n_win_ex);
			ComplexMat temp, temp1, temp2, fft_out, ph_filt;
			Mat wf2, tmp, tmp11, tmp22, H;
			int j_shift, j1, j2; int ret;
			double median = 1.0;
			wf.copyTo(wf2);
			j1 = (ix2 - 1) * n_inc + 1;
			j2 = j1 + n_win - 1;
			if (j2 > n_j)
			{
				j_shift = j2 - n_j;
				j2 = n_j;
				j1 = n_j - n_win + 1;
				tmp11 = Mat::zeros(n_win, j_shift, CV_64F);
				wf2(Range(0, wf2.rows), Range(0, n_win - j_shift)).copyTo(tmp22);
				hconcat(tmp11, tmp22, wf2);
			}
			if (wf2.cols != n_win || wf2.rows != n_win)
			{
				fprintf(stderr, "Goldstein_filter_parallel(): wf2.size and n_win mismatch, please check to make sure n_win is even!\n\n");
				//return -1;
			}
			temp = ph(cv::Range(i1 - 1, i2), cv::Range(j1 - 1, j2));
			ret = ph_bit.SetValue(cv::Range(0, n_win), cv::Range(0, n_win), temp);
			ret = util.fft2(ph_bit, fft_out);
			H = fft_out.GetMod();
			ret = fftshift2(H);
			GaussianFilter(H, H, guasswin);
			//filter2D(H, H, -1, guasswin, Point(-1, -1), 0.0, BORDER_CONSTANT);
			ret = util.ifftshift(H);
			H.copyTo(tmp);
			tmp = tmp.reshape(0, 1);
			cv::sort(tmp, tmp, SORT_ASCENDING + SORT_EVERY_ROW);
			if (tmp.cols % 2 == 1)
			{
				median = tmp.at<double>(0, int((tmp.cols + 1) / 2) - 1);
			}
			else
			{
				median = tmp.at<double>(0, int(tmp.cols / 2) - 1) + tmp.at<double>(0, int(tmp.cols / 2));
				median = median / 2.0;
			}
			if (fabs(median) > 1e-8)
			{
				H = H / median;
			}
			cv::pow(H, alpha, H);
			fft_out = fft_out * H;
			ret = util.ifft2(fft_out, temp);
			Mat re, im;
			re = temp.GetRe();
			im = temp.GetIm();
			temp1 = temp(Range(0, n_win), Range(0, n_win));
			ph_filt = temp1 * wf2;
			temp = ph_out(Range(i1 - 1, i2), Range(j1 - 1, j2)) + ph_filt;
			ret = ph_out.SetValue(Range(i1 - 1, i2), Range(j1 - 1, j2), temp);
		}
		fprintf(stdout, "Goldstein filtering process: %d / %d\n", ix1, n_win_i);
	}
	ph_out.GetPhase().copyTo(phase_filter);
	return 0;
}

// 按二维高斯函数实现高斯滤波
int Filter::GaussianFilter(cv::Mat& src, cv::Mat& dst, cv::Mat window) 
{
	int hh = (window.rows - 1) / 2;
	int hw = (window.cols - 1) / 2;
	Mat Dst = cv::Mat::zeros(src.size(), src.type());
	//边界填充
	cv::Mat Newsrc;
	cv::copyMakeBorder(src, Newsrc, hh, hh, hw, hw, cv::BORDER_REPLICATE);//边界复制
	Dst.zeros(src.size(), src.type());
	//高斯滤波
	for (int i = hh; i < src.rows + hh; ++i) {
		for (int j = hw; j < src.cols + hw; ++j) {
			double sum = 0.0;

			for (int r = -hh; r <= hh; ++r) {
				for (int c = -hw; c <= hw; ++c) {
					sum = sum + Newsrc.ptr<double>(i + r)[j + c] * window.ptr<double>(r + hh)[c + hw];
				}
			}
			Dst.ptr<double>(i - hh)[j - hw] = sum;
		}
	}
	Dst.copyTo(dst);
	return 0;
}

int Filter::GenerateGaussMask(Mat& Mask, int window_height, int window_width, double sigma)
{
	Mask.create(window_height, window_width, CV_64F);
	int h = window_height;
	int w = window_width;
	int center_h = (h - 1) / 2;
	int center_w = (w - 1) / 2;
	double sum = 0.0;
	double x, y;
	for (int i = 0; i < h; ++i) {
		y = pow(i - center_h, 2);
		for (int j = 0; j < w; ++j) {
			x = pow(j - center_w, 2);
			//因为最后都要归一化的，常数部分可以不计算，也减少了运算量
			double g = exp(-(x + y) / (2 * sigma * sigma));
			Mask.ptr<double>(i)[j] = g;
			sum += g;
		}
	}
	Mask = Mask / sum;
	return 0;
}
