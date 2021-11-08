// Unwrap.cpp : 定义 DLL 应用程序的导出函数。
//
#include "stdafx.h"
#include"..\include\Unwrap.h"
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
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Utils.lib")
#else
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
		HANDLE hd = CreateJobObjectA(NULL, "MCF");
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
			nodes[end2 - 1].get_balance() /*&&
			!nodes[end2 - 1].is_residue_node()*/
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
				nodes[end2 - 1].get_balance()/*&&
				!nodes[end2 - 1].is_residue_node()*/
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
		HANDLE hd = CreateJobObjectA(NULL, "MCF");
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
