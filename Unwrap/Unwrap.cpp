// Unwrap.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include"..\include\Unwrap.h"
#include"..\include\FormatConversion.h"
#include<tchar.h>
#include <atlconv.h>
#include<queue>
#define CHECK_RETURN(detail_info) \
{ \
	if (ret < 0) \
	   strncpy( &message[40], detail_info, 150); \
	   fprintf(stderr, "%s\n\n", message); \
	   return -1; \
}
#ifdef _DEBUG
#pragma comment(lib, "ComplexMat_d.lib")
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "FormatConversion_d.lib")
#else
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "FormatConversion.lib")
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
Unwrap::Unwrap()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "UNWRAP_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "UNWRAP_DLL_ERROR: error happens when using parallel computing in function: ");
}

Unwrap::~Unwrap()
{
}

int Unwrap::MCF(
	Mat& wrapped_phase,
	Mat& unwrapped_phase,
	Mat& coherence,
	Mat& residue,
	const char* MCF_problem_file,
	const char* MCF_EXE_PATH
)

{
	if (wrapped_phase.rows < 2 ||
		wrapped_phase.cols < 2 ||
		(wrapped_phase.rows - residue.rows) != 1 ||
		((wrapped_phase.cols - residue.cols)) != 1 ||
		wrapped_phase.rows != coherence.rows ||
		wrapped_phase.cols != coherence.cols ||
		wrapped_phase.type() != CV_64F||
		coherence.type() != CV_64F||
		residue.type() != CV_64F)
	{
		fprintf(stderr, "MCF(): input check failed!\n\n");
		return -1;
	}
	USES_CONVERSION;
	Utils util;
	int ret;
	ret = util.write_DIMACS(MCF_problem_file, residue, coherence, 0.7);
	if (return_check(ret, "write_DIMACS(*, *, *)", error_head)) return -1;
	//////////////////////////创建并调用最小费用流法进程///////////////////////////////
	LPWSTR szCommandLine = new TCHAR[256];
	wcscpy(szCommandLine, A2W(MCF_EXE_PATH));
	wcscat(szCommandLine, L"\\mcf.exe ");
	wcscat(szCommandLine, A2W(MCF_problem_file));
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
		char mcf_job_name[512]; mcf_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(mcf_job_name, "MCF_%lld", tt);
		string mcf_job_name_string(mcf_job_name);
		HANDLE hd = CreateJobObjectA(NULL, mcf_job_name_string.c_str());
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
		fprintf(stderr, "MCF(): create mcf.exe process failed!\n\n");
		if (szCommandLine != NULL) delete[] szCommandLine;
		return -1;
	}
	Mat k1, k2;
	string solution(MCF_problem_file);
	solution.append(".sol");
	ret = util.read_DIMACS(solution.c_str(), k1, k2, wrapped_phase.rows, wrapped_phase.cols);
	if (return_check(ret, "read_DIMACS(*, *, *)", error_head)) return -1;
	Mat diff_1, diff_2;
	ret = util.diff(wrapped_phase, diff_1, diff_2, false);
	if (return_check(ret, "diff(*, *, *, *)", error_head)) return -1;
	ret = util.wrap(diff_1, diff_1);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;
	ret = util.wrap(diff_2, diff_2);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;
	double pi = 3.1415926535;
	diff_1 = diff_1 / (2 * pi);
	diff_2 = diff_2 / (2 * pi);
	diff_1 = diff_1 - k1;
	diff_2 = diff_2 - k2;
	Mat tmp = diff_2(Range(0, 1), Range(0, diff_2.cols));
	copyMakeBorder(tmp, tmp, 0, 0, 1, 0, BORDER_CONSTANT, Scalar(0.0));
	ret = util.cumsum(tmp, 2);
	if (return_check(ret, "cumsum(*, *)", error_head)) return -1;
	copyMakeBorder(diff_1, diff_1, 1, 0, 0, 0, BORDER_CONSTANT, Scalar(0));
	for (int i = 0; i < diff_1.cols; i++)
	{
		diff_1.at<double>(0, i) = tmp.at<double>(0, i);
	}
	ret = util.cumsum(diff_1, 1);
	/*Mat tmp = diff_1(Range(0, diff_1.rows), Range(0, 1));
	copyMakeBorder(tmp, tmp, 1, 0, 0, 0, BORDER_CONSTANT, Scalar(0.0));
	ret = util.cumsum(tmp, 1);
	if (return_check(ret, "cumsum(*, *)", error_head)) return -1;
	copyMakeBorder(diff_2, diff_2, 0, 0, 1, 0, BORDER_CONSTANT, Scalar(0));
	for (int i = 0; i < diff_2.rows; i++)
	{
		diff_2.at<double>(i, 0) = tmp.at<double>(i, 0);
	}
	ret = util.cumsum(diff_2, 2);*/
	if (return_check(ret, "cumsum(*, *)", error_head)) return -1;
	unwrapped_phase = (diff_1)* 2 * pi;
	unwrapped_phase = unwrapped_phase + wrapped_phase.at<double>(0, 0);
	return 0;
}


int Unwrap::MCF(
	Mat& wrapped_phase,
	Mat& unwrapped_phase,
	Mat& mask,
	vector<tri_node>& nodes,
	tri_edge* edges,
	int num_edges, 
	int start,
	bool pass,
	double thresh
)
{
	if (wrapped_phase.rows < 2 ||
		wrapped_phase.cols < 2 ||
		wrapped_phase.type() != CV_64F ||
		wrapped_phase.channels() != 1||
		mask.rows != wrapped_phase.rows||
		mask.cols != wrapped_phase.cols||
		mask.type() != CV_32S||
		mask.channels() != 1||
		nodes.size() < 3||
		edges == NULL||
		num_edges < 1||
		start < 1||
		start > nodes.size()
		)
	{
		fprintf(stderr, "MCF(): input check failed!\n\n");
		return -1;
	}
	wrapped_phase.copyTo(unwrapped_phase);
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2;
	double distance, grad, phi1, phi2, gain, tt, min, max;
	min = 1000000000.0;
	max = -1000000000.0;
	if (pass) tt = 0.5;
	else
	{
		tt = 100000.0;
	}
	long* ptr_neigh = NULL;
	queue<int> que;
	//int start = 1;//起始点默认为第一个点，后续可以自己设定
	ret = nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
	nodes[start - 1].set_status(true);
	for (int i = 0; i < num_neigh; i++)
	{
		if (*(ptr_neigh + i) < 1 || *(ptr_neigh + i) > num_edges)
		{
			fprintf(stderr, "MCF(): edge index exceed legal range!\n");
			return -1;
		}
		end2 = (edges + *(ptr_neigh + i) - 1)->end1 == start ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
		if (end2 < 1 || end2 > num_nodes)
		{
			fprintf(stderr, "MCF(): node index exceed legal range!\n");
			return -1;
		}
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= thresh &&
			fabs((edges + *(ptr_neigh + i) - 1)->gain) < tt &&
			nodes[end2 - 1].get_balance() &&
			!((edges + *(ptr_neigh + i) - 1)->isBoundry && fabs((edges + *(ptr_neigh + i) - 1)->gain) > 0.5)
			)
		{
			que.push(end2);
			//解缠
			nodes[start - 1].get_phase(&phi1);
			nodes[end2 - 1].get_phase(&phi2);
			grad = phi2 - phi1;
			grad = atan2(sin(grad), cos(grad));
			gain = start > end2 ? 2 * PI * (edges + *(ptr_neigh + i) - 1)->gain : -2 * PI * (edges + *(ptr_neigh + i) - 1)->gain;
			nodes[end2 - 1].set_phase(grad + phi1 + gain);
			min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
			max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
			nodes[end2 - 1].set_status(true);
		}
	}
	while (que.size() != 0)
	{
		number = que.front();
		que.pop();
		if (number < 1 || number > num_nodes)
		{
			fprintf(stderr, "MCF(): node index exceed legal range!\n");
			return -1;
		}
		ret = nodes[number - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		for (int i = 0; i < num_neigh; i++)
		{
			if (*(ptr_neigh + i) < 1 || *(ptr_neigh + i) > num_edges)
			{
				fprintf(stderr, "MCF(): edge index exceed legal range!\n");
				return -1;
			}
			end2 = (edges + *(ptr_neigh + i) - 1)->end1 == number ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
			if (end2 < 1 || end2 > num_nodes)
			{
				fprintf(stderr, "MCF(): node index exceed legal range!\n");
				return -1;
			}
			nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
			if (!nodes[end2 - 1].get_status() &&
				distance <= thresh &&
				fabs((edges + *(ptr_neigh + i) - 1)->gain) < tt &&
				nodes[end2 - 1].get_balance()&&
				!((edges + *(ptr_neigh + i) - 1)->isBoundry && fabs((edges + *(ptr_neigh + i) - 1)->gain) > 0.5)
				)
			{
				que.push(end2);
				//解缠
				nodes[number - 1].get_phase(&phi1);
				nodes[end2 - 1].get_phase(&phi2);
				grad = phi2 - phi1;
				grad = atan2(sin(grad), cos(grad));
				gain = number > end2 ? 2 * PI * (edges + *(ptr_neigh + i) - 1)->gain : -2 * PI * (edges + *(ptr_neigh + i) - 1)->gain;
				min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
				max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
				nodes[end2 - 1].set_phase(grad + phi1 + gain);
				nodes[end2 - 1].set_status(true);
			}
		}
	}
	int rows, cols;
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;
	double phi;
	Mat _mask = Mat::zeros(nr, nc, CV_32S);
	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			ret = nodes[i].get_pos(&rows, &cols);
			if (rows > nr - 1 || cols > nc - 1 || rows < 0 || cols < 0)
			{
				fprintf(stderr, "MCF(): node posistion exceed legal range!\n");
				return -1;
			}
			ret = nodes[i].get_phase(&phi);
			unwrapped_phase.at<double>(rows, cols) = phi;
			_mask.at<int>(rows, cols) = 1;
		}
	}
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (_mask.at<int>(i, j) < 1)
			{
				unwrapped_phase.at<double>(i, j) = min - 0.1*(max - min);
			}
		}
	}
	_mask.copyTo(mask);
	return 0;
}

int Unwrap::MCF(
	const Mat& wrapped_phase,
	Mat& unwrapped_phase,
	Mat& out_mask,
	const Mat& mask,
	vector<tri_node>& nodes,
	vector<tri_edge>& edges,
	int start,
	bool pass,
	double thresh
)
{
	if (wrapped_phase.rows < 2 ||
		wrapped_phase.cols < 2 ||
		wrapped_phase.type() != CV_64F ||
		wrapped_phase.channels() != 1 ||
		mask.rows != wrapped_phase.rows ||
		mask.cols != wrapped_phase.cols ||
		mask.type() != CV_32S ||
		mask.channels() != 1 ||
		nodes.size() < 3 ||
		edges.size() < 3 ||
		start < 1 ||
		start > nodes.size()
		)
	{
		fprintf(stderr, "MCF(): input check failed!\n\n");
		return -1;
	}
	if (unwrapped_phase.rows != wrapped_phase.rows ||
		unwrapped_phase.cols != wrapped_phase.cols ||
		unwrapped_phase.type() != CV_64F)
	{
		wrapped_phase.copyTo(unwrapped_phase);
	}
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2;
	double distance, grad, phi1, phi2, gain, tt, min, max;
	min = 1000000000.0;
	max = -1000000000.0;
	if (pass) tt = 0.5;
	else
	{
		tt = 100000.0;
	}
	int num_edges = edges.size();
	long* ptr_neigh = NULL;
	queue<int> que;
	ret = nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
	nodes[start - 1].set_status(true);
	for (int i = 0; i < num_neigh; i++)
	{
		if (*(ptr_neigh + i) < 1 || *(ptr_neigh + i) > num_edges)
		{
			fprintf(stderr, "MCF(): edge index exceed legal range!\n");
			return -1;
		}
		end2 = edges[*(ptr_neigh + i) - 1].end1 == start ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
		if (end2 < 1 || end2 > num_nodes)
		{
			fprintf(stderr, "MCF(): node index exceed legal range!\n");
			return -1;
		}
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= thresh &&
			/*!edges[*(ptr_neigh + i) - 1].isBoundry &&*/
			!(edges[*(ptr_neigh + i) - 1].isBoundry && fabs(edges[*(ptr_neigh + i) - 1].gain) > 0.5) &&
			fabs(edges[*(ptr_neigh + i) - 1].gain) < tt /*&&
			nodes[end2 - 1].get_balance()*/
			)
		{
			que.push(end2);
			//解缠
			nodes[start - 1].get_phase(&phi1);
			nodes[end2 - 1].get_phase(&phi2);
			grad = phi2 - phi1;
			grad = atan2(sin(grad), cos(grad));
			gain = start < end2 ? 2 * PI * edges[*(ptr_neigh + i) - 1].gain : -2 * PI * edges[*(ptr_neigh + i) - 1].gain;
			nodes[end2 - 1].set_phase(grad + phi1 + gain);
			min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
			max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
			nodes[end2 - 1].set_status(true);
		}
	}
	while (que.size() != 0)
	{
		number = que.front();
		que.pop();
		if (number < 1 || number > num_nodes)
		{
			fprintf(stderr, "MCF(): node index exceed legal range!\n");
			return -1;
		}
		ret = nodes[number - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		for (int i = 0; i < num_neigh; i++)
		{
			int end1_row, end2_row, end1_col, end2_col;
			if (*(ptr_neigh + i) < 1 || *(ptr_neigh + i) > num_edges)
			{
				fprintf(stderr, "MCF(): edge index exceed legal range!\n");
				return -1;
			}
			end2 = edges[*(ptr_neigh + i) - 1].end1 == number ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
			if (end2 < 1 || end2 > num_nodes)
			{
				fprintf(stderr, "MCF(): node index exceed legal range!\n");
				return -1;
			}
			nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
			if (!nodes[end2 - 1].get_status() &&
				distance <= thresh &&
				/*!edges[*(ptr_neigh + i) - 1].isBoundry &&*/
				!(edges[*(ptr_neigh + i) - 1].isBoundry && fabs(edges[*(ptr_neigh + i) - 1].gain) > 0.5) &&
				fabs(edges[*(ptr_neigh + i) - 1].gain) < tt/* &&
				nodes[end2 - 1].get_balance() */
				)
			{
				que.push(end2);
				//解缠
				nodes[number - 1].get_phase(&phi1);
				nodes[end2 - 1].get_phase(&phi2);
				grad = phi2 - phi1;
				grad = atan2(sin(grad), cos(grad));
				gain = number < end2 ? 2 * PI * edges[*(ptr_neigh + i) - 1].gain : -2 * PI * edges[*(ptr_neigh + i) - 1].gain;
				min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
				max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
				nodes[end2 - 1].set_phase(grad + phi1 + gain);
				nodes[end2 - 1].set_status(true);
			}
		}
	}
	int rows, cols;
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;
	double phi;
	Mat _mask = Mat::zeros(nr, nc, CV_32S);
	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			ret = nodes[i].get_pos(&rows, &cols);
			ret = nodes[i].get_phase(&phi);
			unwrapped_phase.at<double>(rows, cols) = phi;
			_mask.at<int>(rows, cols) = 1;
		}
	}
//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < nr; i++)
//	{
//		for (int j = 0; j < nc; j++)
//		{
//			if (_mask.at<int>(i, j) < 1)
//			{
//				unwrapped_phase.at<double>(i, j) = min - 0.01 * (max - min);
//			}
//		}
//	}
	_mask.copyTo(out_mask);
	return 0;
}

int Unwrap::MCF_second(Mat& unwrapped_phase, vector<tri_node>& nodes, tri_edge* edges, int num_edges, bool pass, double thresh)
{
	if (unwrapped_phase.rows < 2 ||
		unwrapped_phase.cols < 2 ||
		unwrapped_phase.type() != CV_64F ||
		unwrapped_phase.channels() != 1 ||
		nodes.size() < 3 ||
		edges == NULL ||
		num_edges < 3
		)
	{
		fprintf(stderr, "MCF_second(): input check failed!\n\n");
		return -1;
	}
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2, row_start, col_start;
	double distance, grad, phi1, phi2, gain, tt;
	if (pass) tt = 100000.0;
	else
	{
		tt = 100000.0;
	}
	long* ptr_neigh = NULL;
	queue<int> que;
	queue<int> start_que;
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;
	//未解缠序列号矩阵

	//Mat wrapped_num = Mat::zeros(nr, nc, CV_32S);
	//Mat wrapped_mask = Mat::zeros(nr, nc, CV_32S);
	//Mat grad_updown, grad_leftright;
	//int wrapped_row, wrapped_col;
	//int row_unwrap_start, col_unwrap_start;
	//for (int i = 0; i < num_nodes; i++)
	//{
	//	nodes[i].get_pos(&wrapped_row, &wrapped_col);
	//	wrapped_num.at<int>(wrapped_row, wrapped_col) = i + 1;
	//	wrapped_mask.at<int>(wrapped_row, wrapped_col) = 1;
	//}
	//grad_updown = wrapped_mask(cv::Range(1, nr), cv::Range(0, nc)) - wrapped_mask(cv::Range(0, nr - 1), cv::Range(0, nc));
	//grad_leftright = wrapped_mask(cv::Range(0, nr), cv::Range(1, nc)) - wrapped_mask(cv::Range(0, nr), cv::Range(0, nc - 1));
	//for (int i = 0; i < grad_updown.rows; i++)
	//{
	//	for (int j = 0; j < grad_updown.cols; j++)
	//	{
	//		if (grad_updown.at<int>(i, j) > 0)
	//		{
	//			start_que.push(wrapped_num.at<int>(i + 1, j));
	//		}
	//		if (grad_updown.at<int>(i, j) < 0)
	//		{
	//			start_que.push(wrapped_num.at<int>(i, j));
	//		}
	//	}
	//}
	//for (int i = 0; i < grad_leftright.rows; i++)
	//{
	//	for (int j = 0; j < grad_leftright.cols; j++)
	//	{
	//		if (grad_leftright.at<int>(i, j) > 0)
	//		{
	//			start_que.push(wrapped_num.at<int>(i, j + 1));
	//		}
	//		if (grad_leftright.at<int>(i, j) < 0)
	//		{
	//			start_que.push(wrapped_num.at<int>(i, j));
	//		}
	//	}
	//}
	//int row_search_start, col_search_start, row_search_end, col_search_end;
	//bool b_continue = true;
	//while (start_que.size() != 0)
	//{
	//	number = start_que.front();
	//	start_que.pop();
	//	if (!nodes[number - 1].get_status())
	//	{
	//		nodes[number - 1].get_pos(&row_unwrap_start, &col_unwrap_start);
	//		nodes[number - 1].get_phase(&phi1);
	//		row_search_start = row_unwrap_start - 1 >= 0 ? row_unwrap_start - 1 : row_unwrap_start;
	//		col_search_start = col_unwrap_start - 1 >= 0 ? col_unwrap_start - 1 : col_unwrap_start;
	//		row_search_end = row_unwrap_start + 1 < nr ? row_unwrap_start - 1 : row_unwrap_start;
	//		col_search_end = col_unwrap_start + 1 < nc ? col_unwrap_start - 1 : col_unwrap_start;
	//		b_continue = true;
	//		for (int i = row_search_start; i <= row_search_end; i++)
	//		{
	//			if (b_continue)
	//			{
	//				for (int j = col_search_start; j <= col_search_end; j++)
	//				{
	//					if (wrapped_num.at<int>(i, j) == 0)//找到附近已经解缠的点
	//					{
	//						grad = atan2(sin(phi1 - unwrapped_phase.at<double>(i, j)), cos(phi1 - unwrapped_phase.at<double>(i, j)));
	//						phi1 = unwrapped_phase.at<double>(i, j) + grad;
	//						nodes[number - 1].set_phase(phi1);
	//						nodes[number - 1].set_status(true);
	//						b_continue = false;
	//						que.push(number);
	//						break;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				break;
	//			}
	//		}
	//	}
	//	
	//}

	//找到已解缠邻接节点，并以已解缠邻接节点为起始点开始解缠

	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			que.push(i + 1);
		}
	}
	while (que.size() != 0)
	{
		number = que.front();
		que.pop();
		nodes[number - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			end2 = (edges + *(ptr_neigh + i) - 1)->end1 == number ?
				(edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
			nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
			if (!nodes[end2 - 1].get_status() &&
				distance <= thresh &&
				fabs((edges + *(ptr_neigh + i) - 1)->gain) < tt &&
				!(edges + *(ptr_neigh + i) - 1)->isResidueEdge /*非超过阈值的残差边*/
				)
			{
				que.push(end2);
				nodes[number - 1].get_phase(&phi1);
				nodes[end2 - 1].get_phase(&phi2);
				grad = phi2 - phi1;
				grad = atan2(sin(grad), cos(grad));
				gain = number > end2 ? 2 * PI * (edges + *(ptr_neigh + i) - 1)->gain : -2 * PI * (edges + *(ptr_neigh + i) - 1)->gain;
				nodes[end2 - 1].set_phase(grad + phi1 + gain);
				nodes[end2 - 1].set_status(true);
			}
		}
	}
	
	int row, col;
	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			ret = nodes[i].get_pos(&row, &col);
			ret = nodes[i].get_phase(&phi1);
			unwrapped_phase.at<double>(row, col) = phi1;
		}
	}

	return 0;
}

int Unwrap::mcf_delaunay(const char* MCF_problem_file, const char* MCF_EXE_PATH)
{
	if (MCF_problem_file == NULL ||
		MCF_EXE_PATH == NULL
		)
	{
		fprintf(stderr, "mcf_delaunay(): input check failed!\n\n");
		return -1;
	}
	USES_CONVERSION;
	Utils util;
	//////////////////////////创建并调用最小费用流法进程///////////////////////////////
	LPWSTR szCommandLine = new TCHAR[256];
	wcscpy(szCommandLine, A2W(MCF_EXE_PATH));
	wcscat(szCommandLine, L"\\mcf.exe ");
	wcscat(szCommandLine, A2W(MCF_problem_file));
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
		char mcf_job_name[512]; mcf_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(mcf_job_name, "MCF_%lld", tt);
		string mcf_job_name_string(mcf_job_name);
		HANDLE hd = CreateJobObjectA(NULL, mcf_job_name_string.c_str());
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
		fprintf(stderr, "MCF(): create mcf.exe process failed!\n\n");
		if (szCommandLine != NULL) delete[] szCommandLine;
		return -1;
	}
	return 0;
}

int Unwrap::QualityMap_MCF(Mat& wrapped_phase, Mat& unwrapped_phase, Mat& mask, vector<tri_node>& nodes, tri_edge* edges, int num_edges, int start, bool pass, double thresh)
{
	if (wrapped_phase.rows < 2 ||
		wrapped_phase.cols < 2 ||
		wrapped_phase.type() != CV_64F ||
		wrapped_phase.channels() != 1 ||
		mask.rows != wrapped_phase.rows ||
		mask.cols != wrapped_phase.cols ||
		mask.type() != CV_32S ||
		mask.channels() != 1 ||
		nodes.size() < 3 ||
		edges == NULL ||
		num_edges < 1 ||
		start < 1 ||
		start > nodes.size()
		)
	{
		fprintf(stderr, "MCF(): input check failed!\n\n");
		return -1;
	}
	wrapped_phase.copyTo(unwrapped_phase);
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2;
	double distance, grad, phi1, phi2, gain, tt, min, max;
	min = 1000000000.0;
	max = -1000000000.0;
	if (pass) tt = 0.5;
	else
	{
		tt = 100000.0;
	}
	long* ptr_neigh = NULL;
	//queue<int> que;
	priority_queue<edge_index> neighbour_que;
	edge_index tmp_edge_index;
	bool early_break = false;
	//int start = 1;//起始点默认为第一个点，后续可以自己设定
	if (start > num_nodes) start = 1;
	//////////寻找相关系数最大的边为起始边////////////
	int ix = 0;
	double qua = 100000.0;
	for (int i = 0; i < num_edges; i++)
	{
		if ((edges + i)->quality < qua)
		{
			ix = i;
			qua = (edges + i)->quality;
		}
	}
	start = (edges + ix)->end1;
	ret = nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
	nodes[start - 1].set_status(true);
	for (int i = 0; i < num_neigh; i++)
	{
		end2 = (edges + *(ptr_neigh + i) - 1)->end1 == start ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= thresh &&
			fabs((edges + *(ptr_neigh + i) - 1)->gain) < tt &&
			nodes[end2 - 1].get_balance() /*&&
			!nodes[end2 - 1].is_residue_node()*/
			)
		{
			tmp_edge_index.num = *(ptr_neigh + i);
			tmp_edge_index.quality = (edges + *(ptr_neigh + i) - 1)->quality;
			neighbour_que.push(tmp_edge_index);
		}
	}


	while (neighbour_que.size() != 0)
	{
		tmp_edge_index = neighbour_que.top();
		neighbour_que.pop();
		if (nodes[(edges + tmp_edge_index.num - 1)->end1 - 1].get_status())
		{
			number = (edges + tmp_edge_index.num - 1)->end1;
			end2 = (edges + tmp_edge_index.num - 1)->end2;
		}
		else
		{
			number = (edges + tmp_edge_index.num - 1)->end2;
			end2 = (edges + tmp_edge_index.num - 1)->end1;
		}
		//number = nodes[(edges + tmp_edge_index.num - 1)->end1 - 1].get_status() ? (edges + tmp_edge_index.num - 1)->end2 : (edges + tmp_edge_index.num - 1)->end1;
		//end2 = (edges + tmp_edge_index.num - 1)->end1 == number ? (edges + tmp_edge_index.num - 1)->end2 : (edges + tmp_edge_index.num - 1)->end1;
		nodes[number - 1].get_phase(&phi1);
		nodes[end2 - 1].get_phase(&phi2);
		grad = phi2 - phi1;
		if (!nodes[end2 - 1].get_status() &&
			fabs((edges + tmp_edge_index.num - 1)->gain) < tt &&
			nodes[end2 - 1].get_balance()/*&&
			!nodes[end2 - 1].is_residue_node()*/
			)
		{
			grad = atan2(sin(grad), cos(grad));
			gain = 0.0;
			//gain = number > end2 ? 2 * PI * (edges + tmp_edge_index.num - 1)->gain : -2 * PI * (edges + tmp_edge_index.num - 1)->gain;
			min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
			max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
			nodes[end2 - 1].set_phase(grad + phi1 + gain);
			nodes[end2 - 1].set_status(true);


			ret = nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
			number = end2;
			for (int i = 0; i < num_neigh; i++)
			{

				end2 = (edges + *(ptr_neigh + i) - 1)->end1 == number ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
				nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
				if (!nodes[end2 - 1].get_status() &&
					distance <= thresh &&
					fabs((edges + *(ptr_neigh + i) - 1)->gain) < tt &&
					nodes[end2 - 1].get_balance()/*&&
					!nodes[end2 - 1].is_residue_node()*/
					)
				{
					tmp_edge_index.num = *(ptr_neigh + i);
					tmp_edge_index.quality = (edges + *(ptr_neigh + i) - 1)->quality;
					neighbour_que.push(tmp_edge_index);
				}
			}
		}
		/*else
		{
			nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
			if (nodes[end2 - 1].get_status() && fabs(grad) >= 2 * PI && distance < 1.5)
			{
				early_break = true;
				break;
			}
		}*/
		
	}

	int rows, cols;
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;
	double phi;
	Mat _mask = Mat::zeros(nr, nc, CV_32S);
	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			ret = nodes[i].get_pos(&rows, &cols);
			if (rows > nr - 1 || cols > nc - 1 || rows < 0 || cols < 0)
			{
				fprintf(stderr, "MCF(): node posistion exceed legal range!\n");
				return -1;
			}
			ret = nodes[i].get_phase(&phi);
			unwrapped_phase.at<double>(rows, cols) = phi;
			_mask.at<int>(rows, cols) = 1;
		}
	}
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (_mask.at<int>(i, j) < 1)
			{
				unwrapped_phase.at<double>(i, j) = min - 0.1 * (max - min);
			}
		}
	}
	_mask.copyTo(mask);
	return 0;
}

int Unwrap::_QualityGuided_MCF_1(
	const Mat& wrapped_phase, 
	Mat& unwrapped_phase,
	Mat& out_mask,
	vector<tri_node>& nodes,
	vector<tri_edge>& edges,
	double distance_thresh, 
	bool pass
)
{
	if (wrapped_phase.empty() ||
		wrapped_phase.type() != CV_64F ||
		nodes.size() < 3 ||
		edges.size() < 3 ||
		distance_thresh < 1.0
		)
	{
		fprintf(stderr, "_Quality_MCF_1(): input check failed!\n");
		return -1;
	}

	wrapped_phase.copyTo(unwrapped_phase);
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2, start;
	double distance, grad, phi1, phi2, gain, tt, min, max;
	min = 1000000000.0;
	max = -1000000000.0;
	long* ptr_neigh = NULL;
	priority_queue<edge_index> neighbour_que;
	edge_index tmp_edge_index;
	bool early_break = false;
	size_t num_edges = edges.size();
	//////////寻找相关系数最大的边为起始边(对应quality最小的边)////////////
	int ix = 0;
	double qua = 100000.0;
	for (int i = 0; i < num_edges; i++)
	{
		if (edges[i].quality < qua)
		{
			ix = i;
			qua = edges[i].quality;
		}
	}
	start = edges[ix].end1;

	ret = nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
	nodes[start - 1].set_status(true);
	for (int i = 0; i < num_neigh; i++)
	{
		end2 = edges[*(ptr_neigh + i) - 1].end1 == start ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= distance_thresh &&
			!edges[*(ptr_neigh + i) - 1].isBoundry &&
			nodes[end2 - 1].get_balance() &&
			!edges[*(ptr_neigh + i) - 1].isResidueEdge
			)
		{
			tmp_edge_index.num = *(ptr_neigh + i);
			tmp_edge_index.quality = edges[*(ptr_neigh + i) - 1].quality;
			neighbour_que.push(tmp_edge_index);
		}
	}


	while (neighbour_que.size() != 0)
	{
		tmp_edge_index = neighbour_que.top();
		neighbour_que.pop();
		if (nodes[edges[tmp_edge_index.num - 1].end1 - 1].get_status())
		{
			number = edges[tmp_edge_index.num - 1].end1;
			end2 = edges[tmp_edge_index.num - 1].end2;
		}
		else
		{
			number = edges[tmp_edge_index.num - 1].end2;
			end2 = edges[tmp_edge_index.num - 1].end1;
		}
		if (!nodes[end2 - 1].get_status() &&
			!edges[tmp_edge_index.num - 1].isBoundry &&
			nodes[end2 - 1].get_balance() && 
			!edges[tmp_edge_index.num - 1].isResidueEdge
			)
		{
			nodes[number - 1].get_phase(&phi1);
			nodes[end2 - 1].get_phase(&phi2);
			grad = phi2 - phi1;
			grad = atan2(sin(grad), cos(grad));
			gain = 0.0;
			//min = min > (grad + phi1 + gain) ? (grad + phi1 + gain) : min;
			//max = max < (grad + phi1 + gain) ? (grad + phi1 + gain) : max;
			nodes[end2 - 1].set_phase(grad + phi1 + gain);
			nodes[end2 - 1].set_status(true);


			ret = nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
			number = end2;
			for (int i = 0; i < num_neigh; i++)
			{

				end2 = edges[*(ptr_neigh + i) - 1].end1 == number ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
				nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
				if (!nodes[end2 - 1].get_status() &&
					distance <= distance_thresh &&
					!edges[*(ptr_neigh + i) - 1].isBoundry &&
					nodes[end2 - 1].get_balance()&&
					!edges[*(ptr_neigh + i) - 1].isResidueEdge
					)
				{
					tmp_edge_index.num = *(ptr_neigh + i);
					tmp_edge_index.quality = edges[*(ptr_neigh + i) - 1].quality;
					neighbour_que.push(tmp_edge_index);
				}
			}
		}


	}

	
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;
	double phi;
	Mat _mask = Mat::zeros(nr, nc, CV_32S);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < num_nodes; i++)
	{
		int rows, cols; double phi3;
		if (nodes[i].get_status())
		{
			nodes[i].get_pos(&rows, &cols);
			nodes[i].get_phase(&phi3);
			unwrapped_phase.at<double>(rows, cols) = phi3;
			_mask.at<int>(rows, cols) = 1;
		}
	}
//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < nr; i++)
//	{
//		for (int j = 0; j < nc; j++)
//		{
//			if (_mask.at<int>(i, j) < 1)
//			{
//				unwrapped_phase.at<double>(i, j) = min - 0.1 * (max - min);
//			}
//		}
//	}
	_mask.copyTo(out_mask);

	return 0;
}

int Unwrap::_QualityGuided_MCF_2(
	Mat& unwrapped_phase,
	vector<tri_node>& nodes, 
	vector<tri_edge>& edges,
	double distance_thresh
)
{
	if (unwrapped_phase.rows < 2 ||
		unwrapped_phase.cols < 2 ||
		unwrapped_phase.type() != CV_64F ||
		unwrapped_phase.channels() != 1 ||
		nodes.size() < 3 ||
		edges.size() < 3 
		)
	{
		fprintf(stderr, "_QualityGuided_MCF_2(): input check failed!\n\n");
		return -1;
	}
	int num_nodes = nodes.size();
	int num_neigh, number, ret, end2, row_start, col_start;
	double distance, grad, phi1, phi2, gain, tt;
	long* ptr_neigh = NULL;
	queue<int> que;
	//queue<int> start_que;
	int nr = unwrapped_phase.rows;
	int nc = unwrapped_phase.cols;

	//找到已解缠邻接节点，并以已解缠邻接节点为起始点开始解缠

	for (int i = 0; i < num_nodes; i++)
	{
		if (nodes[i].get_status())
		{
			que.push(i + 1);
		}
	}
	while (que.size() != 0)
	{
		number = que.front();
		que.pop();
		nodes[number - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			end2 = edges[*(ptr_neigh + i) - 1].end1 == number ?
				edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
			nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
			if (!nodes[end2 - 1].get_status() &&
				distance <= distance_thresh &&
				!edges[*(ptr_neigh + i) - 1].isBoundry
				)
			{
				que.push(end2);
				nodes[number - 1].get_phase(&phi1);
				nodes[end2 - 1].get_phase(&phi2);
				grad = phi2 - phi1;
				grad = atan2(sin(grad), cos(grad));
				gain = number < end2 ? 2 * PI * edges[*(ptr_neigh + i) - 1].gain : -2 * PI * edges[*(ptr_neigh + i) - 1].gain;
				nodes[end2 - 1].set_phase(grad + phi1 + gain);
				nodes[end2 - 1].set_status(true);
			}
		}
	}

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < num_nodes; i++)
	{
		int row, col; double phi3;
		if (nodes[i].get_status())
		{
			ret = nodes[i].get_pos(&row, &col);
			ret = nodes[i].get_phase(&phi3);
			unwrapped_phase.at<double>(row, col) = phi3;
		}
	}

	return 0;
}

int Unwrap::QualityGuided_MCF(
	const Mat& wrapped_phase,
	Mat& unwrapped_phase, 
	double coherence_thresh,
	double distance_thresh,
	const char* tmp_path, 
	const char* EXE_path
)
{
	if (wrapped_phase.empty() ||
		tmp_path == NULL ||
		EXE_path == NULL ||
		coherence_thresh < 0.0 ||
		coherence_thresh > 1.0
		)
	{
		fprintf(stderr, "QualityGuided_MCF(): input check failed!\n");
		return -1;
	}
	distance_thresh = distance_thresh < 2.0 ? 2.0 : distance_thresh;
	Mat coherence, phase, mask, quality_index;
	Utils util;
	int ret, nr, nc, count = 0;
	nr = wrapped_phase.rows; nc = wrapped_phase.cols;
	wrapped_phase.copyTo(phase);
	ret = util.phase_coherence(phase, 3, 3, coherence);
	if (return_check(ret, "phase_coherence()", error_head)) return -1;
	quality_index = 1 - coherence;
	//mask = Mat::zeros(nr, nc, CV_32S);

	ret = util.gen_mask(coherence, mask, 7, coherence_thresh);
	if (return_check(ret, "gen_mask()", error_head)) return -1;
	count = cv::countNonZero(mask);
	if (count < 100)//高质量像素小于100，直接使用规则网络的MCF
	{
		Mat residue;
		ret = util.residue(phase, residue);
		if (return_check(ret, "residue()", error_head)) return -1;
		string mcf_problem_file(tmp_path);
		mcf_problem_file.append("\\mcf_problem.net");
		ret = MCF(phase, unwrapped_phase, coherence, residue, mcf_problem_file.c_str(), EXE_path);
		if (return_check(ret, "MCF()", error_head)) return -1;
		return 0;
	}

	string folder(tmp_path);
	string node_file = folder + "\\triangle.node";
	string edge_file = folder + "\\triangle.1.edge";
	string ele_file = folder + "\\triangle.1.ele";
	string neigh_file = folder + "\\triangle.1.neigh";
	string mcf_problem = folder + "\\mcf_delaunay.net";
	string mcf_solution = folder + "\\mcf_delaunay.net.sol";
	vector<tri_node> nodes, nodes_sub; vector<tri_edge> edges, edges_sub; vector<triangle> tri, tri_sub;
	vector<int> node_neighbour, node_neighbour_sub;
	long num_nodes = count;
	ret = util.write_node_file(node_file.c_str(), mask);
	if (return_check(ret, "write_node_file()", error_head)) return -1;
	ret = util.gen_delaunay(node_file.c_str(), EXE_path);
	if (return_check(ret, "gen_delaunay()", error_head)) return -1;
	ret = util.read_edges(edge_file.c_str(), edges, node_neighbour, num_nodes);
	if (return_check(ret, "read_edges()", error_head)) return -1;
	ret = util.init_tri_node(nodes, phase, mask, edges, node_neighbour, num_nodes);
	if (return_check(ret, "init_tri_node()", error_head)) return -1;
	ret = util.init_edges_quality(quality_index, edges, nodes);
	if (return_check(ret, "init_edges_quality()", error_head)) return -1;
	ret = util.read_triangle(ele_file.c_str(), neigh_file.c_str(), tri, nodes, edges);
	if (return_check(ret, "read_triangle()", error_head)) return -1;
	ret = util.residue(tri, nodes, edges, distance_thresh);
	if (return_check(ret, "residue()", error_head)) return -1;

	Mat out_mask;
	ret = _QualityGuided_MCF_1(phase, unwrapped_phase, out_mask, nodes, edges, distance_thresh);
	if (return_check(ret, "residue()", error_head)) return -1;

	out_mask.convertTo(mask, CV_64F);
	util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\out_mask.bin", mask);
	out_mask.convertTo(mask, CV_32S);
	util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\unwrapped_phase1.bin", unwrapped_phase);

	Mat mask_2 = Mat::zeros(nr, nc, CV_32S);
	count = 0;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (out_mask.at<int>(i, j) == 0)
			{
				mask_2.at<int>(i, j) = 1; count++;
			}
		}
	}
	Mat _mask_sentinel; mask_2.copyTo(_mask_sentinel);//质量图解缠未解出的区域掩膜
	if (count == 0) return 0;
	//找出邻接已解缠节点
	Mat grad_updown, grad_leftright;
	grad_updown = mask_2(cv::Range(1, nr), cv::Range(0, nc)) - mask_2(cv::Range(0, nr - 1), cv::Range(0, nc));
	grad_leftright = mask_2(cv::Range(0, nr), cv::Range(1, nc)) - mask_2(cv::Range(0, nr), cv::Range(0, nc - 1));
	for (int i = 0; i < grad_updown.rows; i++)
	{
		for (int j = 0; j < grad_updown.cols; j++)
		{
			if (grad_updown.at<int>(i, j) > 0)
			{
				mask_2.at<int>(i, j) = 2;//2代表邻接已解缠节点
			}
			if (grad_updown.at<int>(i, j) < 0)
			{
				mask_2.at<int>(i + 1, j) = 2;
			}
		}
	}
	for (int i = 0; i < grad_leftright.rows; i++)
	{
		for (int j = 0; j < grad_leftright.cols; j++)
		{
			if (grad_leftright.at<int>(i, j) > 0)
			{
				mask_2.at<int>(i, j) = 2;
			}
			if (grad_leftright.at<int>(i, j) < 0)
			{
				mask_2.at<int>(i, j + 1) = 2;
			}
		}
	}

	/*
	* 找到未解缠的点以及其相邻已解缠的点，形成三角网络
	加入总的队列wrapped_que，已解缠的边缘点加入到unwrapped_neighbour_que
	*/
	queue<int> wrapped_que, unwrapped_neighbour_que, low_quality_que;
	num_nodes = cv::countNonZero(mask_2);
	ret = util.write_node_file(node_file.c_str(), mask_2);
	if (return_check(ret, "write_node_file()", error_head)) return -1;
	ret = util.gen_delaunay(node_file.c_str(), EXE_path);
	if(return_check(ret, "gen_delaunay()", error_head)) return -1;
	ret = util.read_edges(edge_file.c_str(), edges, node_neighbour, num_nodes);
	if (return_check(ret, "read_edges()", error_head)) return -1;
	ret = util.init_tri_node(nodes, unwrapped_phase, mask_2, edges, node_neighbour, num_nodes);
	if (return_check(ret, "init_tri_node()", error_head)) return -1;

	int start, end1, end2, ambig, i;
	long* ptr_neigh = NULL; int num_neigh, row, col, num_triangle, positive, negative;
	double distance, phi, cluster_distance_thresh = 1.2;//低质量聚类距离阈值
	Mat zeros = Mat::zeros(nr, nc, CV_32S);
	Mat ambiguity, new_mask;
	int  mask_sentinel_new = cv::countNonZero(_mask_sentinel);
	while (mask_sentinel_new != 0)//只要还存在未解缠的像素就继续循环
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i].get_pos(&row, &col);
			if (_mask_sentinel.at<int>(row, col) == 1)
			{
				wrapped_que.push(i + 1);
				break;
			}
		}
		//mask_sentinel_old = cv::countNonZero(_mask_sentinel);
		zeros.copyTo(new_mask);//循环前将new_mask清零
		while (!wrapped_que.empty())//寻找低质量点cluster
		{
			
			start = wrapped_que.front();
			wrapped_que.pop();
			if (nodes[start - 1].get_status())
			{
				unwrapped_neighbour_que.push(start);
			}
			nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			nodes[start - 1].get_pos(&row, &col);//设置新的mask
			new_mask.at<int>(row, col) = 1;
			_mask_sentinel.at<int>(row, col) = 0;//未解缠的像素掩膜更新
			nodes[start - 1].set_balance(false);//已加入队列设置为不平衡，避免重复加入队列
			for (int i = 0; i < num_neigh; i++)
			{
				end2 = edges[*(ptr_neigh + i) - 1].end1 == start ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
				nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
				if (
					distance <= cluster_distance_thresh && //小于低质量聚类距离阈值则为同一类
					nodes[end2 - 1].get_balance()
					)
				{
					wrapped_que.push(end2);
					nodes[end2 - 1].set_balance(false);//已加入队列设置为不平衡，避免重复加入队列
				}
			}
		}

		//new_mask.convertTo(new_mask, CV_64F);
		//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\out_mask.bin", new_mask);
		//new_mask.convertTo(new_mask, CV_32S);

		ambiguity = Mat::zeros(1, unwrapped_neighbour_que.size(), CV_32S);
		num_nodes = cv::countNonZero(new_mask);
		ret = util.write_node_file(node_file.c_str(), new_mask);
		ret = util.gen_delaunay(node_file.c_str(), EXE_path);
		ret = util.read_edges(edge_file.c_str(), edges_sub, node_neighbour_sub, num_nodes);
		ret = util.init_tri_node(nodes_sub, wrapped_phase, new_mask, edges_sub, node_neighbour_sub, num_nodes);
		ret = util.read_triangle(ele_file.c_str(), neigh_file.c_str(), tri_sub, nodes_sub, edges_sub);
		ret = util.residue(tri_sub, nodes_sub, edges_sub, 1000.0);
		if (return_check(ret, "residue()", error_head)) return -1;
		/*
		* 检查残差点数，若无残差点则不使用mcf.exe求解
		*/
		num_triangle = tri_sub.size(); positive = 0; negative = 0;
		for (int ii = 0; ii < num_triangle; ii++)
		{
			if (tri_sub[ii].residue > 0.7)
			{
				positive++;
			}
			if (tri_sub[ii].residue < -0.7)
			{
				negative++;
			}
		}
		if (positive == 0 && negative == 0)
		{

			ret = MCF(wrapped_phase, unwrapped_phase, out_mask, new_mask, nodes_sub, edges_sub, 1, false, distance_thresh);
			if (return_check(ret, "MCF()", error_head)) return -1;
		}
		else
		{
			ret = util.write_DIMACS(mcf_problem.c_str(), tri_sub, nodes_sub, edges_sub, coherence);
			if (return_check(ret, "write_DIMACS()", error_head)) return -1;
			ret = mcf_delaunay(mcf_problem.c_str(), EXE_path);
			if (return_check(ret, "mcf_delaunay()", error_head)) return -1;
			ret = util.read_DIMACS(mcf_solution.c_str(), edges_sub, nodes_sub, tri_sub);
			if (return_check(ret, "read_DIMACS()", error_head)) return -1;
			ret = MCF(wrapped_phase, unwrapped_phase, out_mask, new_mask, nodes_sub, edges_sub, 1, false, distance_thresh);
			if (return_check(ret, "MCF()", error_head)) return -1;
		}

		//校正模糊数
		i = 0;
		while (!unwrapped_neighbour_que.empty())
		{
			end2 = unwrapped_neighbour_que.front();
			unwrapped_neighbour_que.pop();
			nodes[end2 - 1].get_pos(&row, &col);
			nodes[end2 - 1].get_phase(&phi);
			ambig = (int)round((phi - unwrapped_phase.at<double>(row, col)) / (2 * 3.141592653589793238));
			ambiguity.at<int>(0, i++) = ambig;
		}
		ret = util.get_mode_index(ambiguity, &ambig);
		//校正相位
		for (int i = 0; i < nodes_sub.size(); i++)
		{
			if (nodes_sub[i].get_status())
			{
				nodes_sub[i].get_phase(&phi);
				nodes_sub[i].get_pos(&row, &col);
				phi += (double)ambig * 2 * 3.141592653589793238;
				nodes_sub[i].set_phase(phi);
				unwrapped_phase.at<double>(row, col) = phi;
			}
		}
		mask_sentinel_new = cv::countNonZero(_mask_sentinel);
	}
	


	return 0;
}

int Unwrap::snaphu(
	const char* wrapped_phase_file,
	Mat& unwrapped_phase,
	const char* project_path,
	const char* tmp_folder,
	const char* exe_path
)
{
	if (wrapped_phase_file == NULL || 
		project_path == NULL ||
		tmp_folder == NULL ||
		exe_path == NULL)
	{
		fprintf(stderr, "snaphu(): input check failed!\n");
		return -1;
	}

	FormatConversion conversion;
	Utils util;
	int ret, nr, nc;
	double B_effect, B_parallel;
	Mat wrapped_phase, coherence, amplitude1, amplitude2, lon_coef, lat_coef, state_vec1, state_vec2, prf1, prf2,
		carrier_frequency, offset_row, offset_col;
	ComplexMat master, slave;
	FILE* fp = NULL;
	string EXE_path(exe_path);
	std::replace(EXE_path.begin(), EXE_path.end(), '/', '\\');
	string project(project_path);
	std::replace(project.begin(), project.end(), '/', '\\');
	string source_1, source_2;
	string folder(tmp_folder);
	std::replace(folder.begin(), folder.end(), '/', '\\');
	string config_file = folder + "\\snaphu.config";
	string ampfile1 = folder + "\\ampfile1.dat";
	string ampfile2 = folder + "\\ampfile2.dat";
	string coherence_file = folder + "\\coherence.dat";
	string IN_file = folder + "\\wrapped_phase.dat";
	string OUT_file = folder + "\\unwrapped_phase.dat";

	ret = conversion.read_array_from_h5(wrapped_phase_file, "phase", wrapped_phase);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nr = wrapped_phase.rows; nc = wrapped_phase.cols;
	//估计基线
	bool b_baseline = true;
	bool b_source = true;
	bool b_amp = true;
	bool b_coh = true;
	if (0 > conversion.read_str_from_h5(wrapped_phase_file, "source_1", source_1)) b_source = false;
	else source_1 = project + source_1;
	if (0 > conversion.read_str_from_h5(wrapped_phase_file, "source_2", source_2)) b_source = false;
	else source_2 = project + source_2;

	if (0 > conversion.read_array_from_h5(wrapped_phase_file, "coherence", coherence)) b_coh = false;
	if (b_source)
	{
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "lon_coefficient", lon_coef)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "lat_coefficient", lat_coef)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "offset_row", offset_row)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "offset_col", offset_col)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "carrier_frequency", carrier_frequency)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "prf", prf1)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_1.c_str(), "state_vec", state_vec1)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_2.c_str(), "prf", prf2)) b_baseline = false;
		if (0 > conversion.read_array_from_h5(source_2.c_str(), "state_vec", state_vec2)) b_baseline = false;

		if (0 > conversion.read_slc_from_h5(source_1.c_str(), master)) b_amp = false;
		if (0 > conversion.read_slc_from_h5(source_2.c_str(), slave)) b_amp = false;
	}
	if (b_source && b_baseline)
	{
		ret = util.baseline_estimation(state_vec1, state_vec2, lon_coef, lat_coef, offset_row.at<int>(0, 0),
			offset_col.at<int>(0, 0), nr, nc, 1.0 / prf1.at<double>(0, 0), 1.0 / prf2.at<double>(0, 0), &B_effect, &B_parallel);
		if (ret < 0) b_baseline = false;
	}
	Mat phase;
	wrapped_phase.convertTo(phase, CV_32F);
	fopen_s(&fp, IN_file.c_str(), "wb");
	if (!fp)
	{
		fprintf(stderr, "snaphu(): can't open %s!\n", IN_file.c_str());
		return -1;
	}
	fwrite(phase.data, sizeof(float), nr * nc, fp);
	fclose(fp);
	fp = NULL;
	if (b_source && b_amp)//有幅度信息
	{
		if (master.type() != CV_64F) master.convertTo(master, CV_64F);
		amplitude1 = master.GetMod();
		amplitude1.convertTo(amplitude1, CV_32F);
		fopen_s(&fp, ampfile1.c_str(), "wb");
		if (!fp)
		{
			fprintf(stderr, "snaphu(): can't open %s!\n", ampfile1.c_str());
			return -1;
		}
		fwrite(amplitude1.data, sizeof(float), nr * nc, fp);
		fclose(fp);
		fp = NULL;

		if (slave.type() != CV_64F) master.convertTo(slave, CV_64F);
		amplitude1 = slave.GetMod();
		amplitude1.convertTo(amplitude1, CV_32F);
		fopen_s(&fp, ampfile2.c_str(), "wb");
		if (!fp)
		{
			fprintf(stderr, "snaphu(): can't open %s!\n", ampfile2.c_str());
			return -1;
		}
		fwrite(amplitude1.data, sizeof(float), nr * nc, fp);
		fclose(fp);
		fp = NULL;
	}
	if (b_coh)//相关系数信息
	{
		coherence.convertTo(coherence, CV_32F);
		fopen_s(&fp, coherence_file.c_str(), "wb");
		if (!fp)
		{
			fprintf(stderr, "snaphu(): can't open %s!\n", coherence_file.c_str());
			return -1;
		}
		fwrite(coherence.data, sizeof(float), nr * nc, fp);
		fclose(fp);
		fp = NULL;
		b_coh = true;
	}
	else
	{
		ret = util.phase_coherence(wrapped_phase, coherence);
		if (ret < 0) b_coh = false;
		else
		{
			coherence.convertTo(coherence, CV_32F);
			fopen_s(&fp, coherence_file.c_str(), "wb");
			if (!fp)
			{
				fprintf(stderr, "snaphu(): can't open %s!\n", coherence_file.c_str());
				return -1;
			}
			fwrite(coherence.data, sizeof(float), nr * nc, fp);
			fclose(fp);
			fp = NULL;
			b_coh = true;
		}
	}



	//写入配置参数
	fopen_s(&fp, config_file.c_str(), "wt");
	if (!fp)
	{
		fprintf(stderr, "snaphu(): can't open %s!\n", coherence_file.c_str());
		return -1;
	}
	fprintf(fp, "INFILEFORMAT FLOAT_DATA\n");
	fprintf(fp, "OUTFILEFORMAT FLOAT_DATA\n");
	fprintf(fp, "CORRFILEFORMAT FLOAT_DATA\n");
	fprintf(fp, "AMPFILEFORMAT FLOAT_DATA\n");
	fprintf(fp, "LINELENGTH %d\n", nc);
	fprintf(fp, "INFILE %s\n", IN_file.c_str());
	fprintf(fp, "OUTFILE %s\n", OUT_file.c_str());
	if (b_coh) fprintf(fp, "CORRFILE %s\n", coherence_file.c_str());
	if (b_source && b_amp)
	{
		fprintf(fp, "AMPFILE1 %s\n", ampfile1.c_str());
		fprintf(fp, "AMPFILE2 %s\n", ampfile2.c_str());
	}
	if (b_source && b_baseline)
	{
		B_parallel = sqrt(B_effect * B_effect + B_parallel * B_parallel);
		fprintf(fp, "BASELINE %lf\n", B_parallel);
		fprintf(fp, "BPERP %lf\n", B_effect);
		fprintf(fp, "LAMBDA %lf\n", 3e8 / carrier_frequency.at<double>(0, 0));
	}
	if (b_source)
	{
		Mat DR, DA;
		if (0 == conversion.read_array_from_h5(source_1.c_str(), "range_spacing", DR))
		{
			fprintf(fp, "DR %lf\n", DR.at<double>(0, 0));
		}
		if (0 == conversion.read_array_from_h5(source_1.c_str(), "azimuth_spacing", DA))
		{
			fprintf(fp, "DA %lf\n", DA.at<double>(0, 0));
		}
		if (0 == conversion.read_array_from_h5(source_1.c_str(), "range_resolution", DR))
		{
			fprintf(fp, "RANGERES %lf\n", DR.at<double>(0, 0));
		}
		if (0 == conversion.read_array_from_h5(source_1.c_str(), "azimuth_resolution", DA))
		{
			fprintf(fp, "AZRES %lf\n", DA.at<double>(0, 0));
		}
	}

	fclose(fp);
	fp = NULL;

	USES_CONVERSION;
	//////////////////////////创建并调用snaphu.exe进程///////////////////////////////
	LPWSTR szCommandLine = new TCHAR[256];
	wcscpy(szCommandLine, A2W(EXE_path.c_str()));
	wcscat(szCommandLine, L"\\snaphu.exe -f ");
	wcscat(szCommandLine, A2W(config_file.c_str()));
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
		char snaphu_job_name[512]; snaphu_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(snaphu_job_name, "SNAPHU_%lld", tt);
		string snaphu_job_name_string(snaphu_job_name);
		HANDLE hd = CreateJobObjectA(NULL, snaphu_job_name_string.c_str());
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
		fprintf(stderr, "snaphu(): create snaphu.exe process failed!\n\n");
		if (szCommandLine != NULL) delete[] szCommandLine;
		return -1;
	}

	//读取结果
	unwrapped_phase.create(nr, nc, CV_32F);
	fopen_s(&fp, OUT_file.c_str(), "rb");
	if (!fp)
	{
		fprintf(stderr, "snaphu(): can't open %s!\n", OUT_file.c_str());
		return -1;
	}
	fread(unwrapped_phase.data, sizeof(float), nr * nc, fp);
	fclose(fp);
	fp = NULL;
	unwrapped_phase.convertTo(unwrapped_phase, CV_64F);
	return 0;
}

int Unwrap::GetSPD(Mat& wrapped_phase, Mat& SPD)
{
	if (wrapped_phase.rows < 1 ||
		wrapped_phase.cols < 1 ||
		wrapped_phase.type() != CV_64F ||
		wrapped_phase.channels() != 1)
	{
		fprintf(stderr, "GetSPD(): input check failed!\n\n");
		return -1;
	}
	int win_w = 3;
	int win_h = 3;
	int armw = win_w / 2;
	int armh = win_h / 2;
	Mat padded;
	copyMakeBorder(wrapped_phase, padded, armh, armh, armw, armw, BORDER_REFLECT_101);//镜像翻转边缘
	int width = padded.cols;
	int height = padded.rows;
	int i, j;
	volatile bool parallel_flag = true;
	int ret = 0;
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (i = armh; i < height - armh; i++)
	{
		if (!parallel_flag) continue;
		for (j = armw; j < width - armw; j++)
		{
			if (!parallel_flag) continue;
			int m, n;
			double sum = 0;
			double delta = 0;
			/*在3*3的窗口内计算与中心像素的梯度绝对值和*/
			for (m = i - 1; m < i + 2; m++)
				for (n = j - 1; n < j + 2; n++)
				{
					delta = padded.ptr<double>(i)[j] - padded.ptr<double>(m)[n];
					/*梯度取主值*/
					if (delta <= -PI)
						sum += abs(delta + 2 * PI);
					else if (delta > -PI && delta < PI)
						sum += abs(delta);
					else
						sum += abs(delta - 2 * PI);
				}
			///*PSD*/
			//double k = mean(padded(Range(i - 1, i + 2), Range(i - 1, i + 2)))[0];
			//for (m = i - 1; m < i + 2; m++)
			//	for (n = j - 1; n < j + 2; n++)
			//	{
			//		sum += pow((padded.ptr<double>(i)[j] - k), 2);
			//		/*梯度取主值*/
			//		/*if (delta <= -PI)
			//			sum += abs(delta + 2 * PI);
			//		else if (delta > -PI && delta < PI)
			//			sum += abs(delta);
			//		else
			//			sum += abs(delta - 2 * PI);*/
			//	}
			SPD.ptr<double>(i - armh)[j - armw] = sqrt(sum / 8);
		}
	}
	if (parallel_check(parallel_flag, "GetSPD()", parallel_error_head)) return -1;
	return 0;
}

int Unwrap::unwrap(Mat& src, Mat& dst, int x0, int y0, int x1, int y1, Mat& flag, Mat& adjoin, Mat& SPD, Heap& Q)
{
	if (src.rows < 1 ||
		src.cols < 1 ||
		src.type() != CV_64F ||
		src.channels() != 1 ||
		x0 < 0 || x0 >= src.cols || y0 < 0 || y0 >= src.rows ||
		x1 < 0 || x1 >= src.cols || y1 < 0 || y1 >= src.rows ||
		flag.size() != src.size() ||
		flag.type() != CV_64F ||
		flag.channels() != 1 ||
		adjoin.size() != src.size() ||
		adjoin.type() != CV_64F ||
		adjoin.channels() != 1 ||
		SPD.size() != src.size() ||
		SPD.type() != CV_64F ||
		SPD.channels() != 1
		)
	{
		fprintf(stderr, "unwrap(): input check failed!\n\n");
		return -1;
	}
	int width = src.cols;
	int height = src.rows;
	double delta = atan2(sin(src.ptr<double>(y0)[x0] - src.ptr<double>(y1)[x1]), cos(src.ptr<double>(y0)[x0] - src.ptr<double>(y1)[x1]));
	dst.ptr<double>(y1)[x1] = dst.ptr<double>(y0)[x0] - delta;
	flag.ptr<double>(y1)[x1] = 0;
	adjoin.ptr<double>(y1)[x1] = 0;
	int ret;
	if (x1 > 0)
		if ((flag.ptr<double>(y1)[x1 - 1] == 1) && (adjoin.ptr<double>(y1)[x1 - 1] == 0))
		{
			adjoin.ptr<double>(y1)[x1 - 1] = 1;
			ret = Q.push(SPD.ptr<double>(y1)[x1 - 1], x1 - 1, y1);
			if (ret < 0)
				return -1;
		}

	if (x1 < width - 1)
		if ((flag.ptr<double>(y1)[x1 + 1] == 1) && (adjoin.ptr<double>(y1)[x1 + 1] == 0))
		{
			adjoin.ptr<double>(y1)[x1 + 1] = 1;
			ret = Q.push(SPD.ptr<double>(y1)[x1 + 1], x1 + 1, y1);
			if (ret < 0)
				return -1;
		}
	if (y1 > 0)
		if ((flag.ptr<double>(y1 - 1)[x1] == 1) && (adjoin.ptr<double>(y1 - 1)[x1] == 0))
		{
			adjoin.ptr<double>(y1 - 1)[x1] = 1;
			ret = Q.push(SPD.ptr<double>(y1 - 1)[x1], x1, y1 - 1);
			if (ret < 0)
				return -1;
		}

	if (y1 < height - 1)
		if ((flag.ptr<double>(y1 + 1)[x1] == 1) && (adjoin.ptr<double>(y1 + 1)[x1] == 0))
		{
			adjoin.ptr<double>(y1 + 1)[x1] = 1;
			ret = Q.push(SPD.ptr<double>(y1 + 1)[x1], x1, y1 + 1);
			if (ret < 0)
				return -1;
		}
	return 0;
}

int Unwrap::SPD_Guided_Unwrap(Mat& wrapped_phase, Mat& unwrapped_phase)
{
	if (wrapped_phase.rows < 1 ||
		wrapped_phase.cols < 1 ||
		wrapped_phase.type() != CV_64F ||
		wrapped_phase.channels() != 1)
	{
		fprintf(stderr, "SPD_Guided_Unwrap(): input check failed!\n\n");
		return -1;
	}
	Heap Heap;
	Mat tmp = Mat::zeros(wrapped_phase.size(), CV_64FC1);
	Mat SPD = Mat::zeros(wrapped_phase.size(), CV_64FC1);
	int ret = GetSPD(wrapped_phase, SPD);
	if (ret < 0)
		return -1;
	int count = 0;
	int width = wrapped_phase.cols;
	int height = wrapped_phase.rows;
	Mat flag = Mat::ones(wrapped_phase.size(), CV_64FC1);
	Mat adjoin = Mat::zeros(wrapped_phase.size(), CV_64FC1);
	Point Min;
	minMaxLoc(SPD, NULL, NULL, &Min, NULL);
	int x = Min.x;
	int y = Min.y;
	int mark = 0;
	tmp.ptr<double>(y)[x] = tmp.ptr<double>(y)[x];
	count++;
	flag.ptr<double>(y)[x] = 0;
	if (x > 0)
	{
		ret = unwrap(wrapped_phase, tmp, x, y, x - 1, y, flag, adjoin, SPD, Heap);
		if (ret < 0)
			return -1;
		count++;
	}

	if (x < width - 1)
	{
		ret = unwrap(wrapped_phase, tmp, x, y, x + 1, y, flag, adjoin, SPD, Heap);
		if (ret < 0)
			return -1;
		count++;
	}
	if (y > 0)
	{
		ret = unwrap(wrapped_phase, tmp, x, y, x, y - 1, flag, adjoin, SPD, Heap);
		if (ret < 0)
			return -1;
		count++;
	}

	if (y < height - 1)
	{
		ret = unwrap(wrapped_phase, tmp, x, y, x, y + 1, flag, adjoin, SPD, Heap);
		if (ret < 0)
			return -1;
		count++;
	}
	int top = height * width;
	while (count < top)
	{
		mark = 0;
		if (Heap.size != 0)
		{
			ret = Heap.top(&x, &y);
			if (ret < 0)
				return -1;
			Heap.pop();
		}
		if (y > 0)
		{
			if (flag.ptr<double>(y - 1)[x] == 0)
			{
				ret = unwrap(wrapped_phase, tmp, x, y - 1, x, y, flag, adjoin, SPD, Heap);
				if (ret < 0)
					return -1;
				mark = 1;
			}
		}
		if (y < height - 1 && mark == 0)
		{
			if (flag.ptr<double>(y + 1)[x] == 0)
			{
				ret = unwrap(wrapped_phase, tmp, x, y + 1, x, y, flag, adjoin, SPD, Heap);
				if (ret < 0)
					return -1;
				mark = 1;
			}
		}
		if (x > 0 && mark == 0)
		{
			if (flag.ptr<double>(y)[x - 1] == 0)
			{
				ret = unwrap(wrapped_phase, tmp, x - 1, y, x, y, flag, adjoin, SPD, Heap);
				if (ret < 0)
					return -1;
				mark = 1;
			}
		}
		if (x < width - 1 && mark == 0)
		{
			if (flag.ptr<double>(y)[x + 1] == 0)
			{
				ret = unwrap(wrapped_phase, tmp, x + 1, y, x, y, flag, adjoin, SPD, Heap);
				if (ret < 0)
					return -1;
				mark = 1;
			}
		}
		if (mark == 1)
			count++;
	}
	tmp.copyTo(unwrapped_phase);
	free(Heap.x);
	free(Heap.y);
	free(Heap.queue);
	return 0;
}