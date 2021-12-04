// Utils.cpp : 定义 DLL 应用程序的导出函数。
//
#include<complex.h>
#include "stdafx.h"
#include"..\include\Utils.h"
#include"..\include\Registration.h"
#include"..\include\Deflat.h"
#include"..\include\Unwrap.h"
#include<tchar.h>
#include <atlconv.h>
#include"../include/FormatConversion.h"
#include"Eigen/Dense"

#ifdef _DEBUG
#pragma comment(lib,"ComplexMat_d.lib")
#pragma comment(lib, "Registration_d.lib")
#pragma comment(lib, "FormatConversion_d.lib")
#pragma comment(lib, "Deflat_d.lib")
#pragma comment(lib, "Unwrap_d.lib")
#else
#pragma comment(lib,"ComplexMat.lib")
#pragma comment(lib, "Registration.lib")
#pragma comment(lib, "FormatConversion.lib")
#pragma comment(lib, "Deflat.lib")
#pragma comment(lib, "Unwrap.lib")
#endif // _DEBUG

using namespace cv;
/*宏定义*/
#define RETURN_MSG \
{ \
    if( fp ) fclose( fp ); \
    return( -1 ); \
}

#define GET_NEXT_LINE \
{ \
    if( !fgets( instring, 256, fp ) ) \
        ch = 0; \
	    else \
        ch = *instring; \
}
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

inline bool read_check(long ret, long ret_ref, const char* detail_info, const char* error_head)
{
	if (ret != ret_ref)
	{
		fprintf(stderr, "%s %s\n\n", error_head, detail_info);
		return true;
	}
	return false;
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
Utils::Utils()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "UTILS_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "UTILS_DLL_ERROR: error happens when using parallel computing in function: ");
}

Utils::~Utils()
{
}

int Utils::get_mode_index(const Mat& input, int* out)
{
	if (input.empty() || input.type() != CV_32S || out == NULL)
	{
		fprintf(stderr, "get_mode_index(): input check failed!\n");
		return -1;
	}
	Mat temp; input.copyTo(temp);
	int nr = temp.rows; int nc = temp.cols;
	temp = temp.reshape(0, 1);
	cv::sort(temp, temp, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);
	int total = nr * nc;
	int max_count = 0; int max_count_ix = 0; int count = 0, i = 0, j = 0;
	while (i < total - 1)
	{
		count = 0;
		for (j = i; j < total - 1; j++)
		{
			if (temp.at<int>(0, j) == temp.at<int>(0, j + 1)) count++;
			else break;
		}
		if (max_count < count)
		{
			max_count = count;
			max_count_ix = j;
		}
		j++;
		i = j;
	}
	*out = temp.at<int>(0, max_count_ix);
	return 0;
}

int Utils::diff(Mat& Src, Mat& diff1, Mat& diff2, bool same)
{
	int nr = Src.rows;
	int nc = Src.cols;
	if (nr < 2 || nc < 2 || Src.type() != CV_64F)
	{
		fprintf(stderr, "diff(): input check failed!\n\n");
		return -1;
	}

	diff1 = Src(Range(1, nr), Range(0, nc)) - Src(Range(0, nr - 1), Range(0, nc));
	diff2 = Src(Range(0, nr), Range(1, nc)) - Src(Range(0, nr), Range(0, nc - 1));
	if (same)
	{
		copyMakeBorder(diff1, diff1, 0, 1, 0, 0, BORDER_CONSTANT, Scalar(0.0));
		copyMakeBorder(diff2, diff2, 0, 0, 0, 1, BORDER_CONSTANT, Scalar(0.0));
	}
	return 0;
}



int Utils::generate_phase(const ComplexMat& Master, const ComplexMat& Slave, Mat& phase)
{
	if (Master.GetRows() < 1 ||
		Master.GetCols() < 1 ||
		Slave.GetRows() != Master.GetRows() ||
		Slave.GetCols() != Master.GetCols())
	{
		fprintf(stderr, "generate_phase(): input check failed!\n\n");
		return -1;
	}
	ComplexMat tmp;
	int ret = Master.Mul(Slave, tmp, true);
	phase = tmp.GetPhase();
	return 0;
}

int Utils::write_DIMACS(const char* DIMACS_file_problem, triangle* tri, int num_triangle, vector<tri_node>& nodes, tri_edge* edges, long num_edges, Mat& cost)
{
	if (DIMACS_file_problem == NULL ||
		tri == NULL ||
		num_triangle < 1 ||
		nodes.size() < 3 ||
		edges == NULL ||
		num_edges < 3||
		cost.rows < 2||
		cost.cols < 2||
		cost.channels() != 1||
		cost.type() != CV_64F
		)
	{
		fprintf(stderr, "write_DIMACS(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(DIMACS_file_problem, "wt");
	if (fp == NULL)
	{
		fprintf(stderr, "write_DIMACS(): can't open %s\n", DIMACS_file_problem);
		return -1;
	}
	
	int ret, num_nodes;

	num_nodes = nodes.size();
	long num_arcs = 0;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i) != NULL)
		{
			if ((tri + i)->neigh1 > 0) num_arcs++;
			if ((tri + i)->neigh2 > 0) num_arcs++;
			if ((tri + i)->neigh3 > 0) num_arcs++;
		}
	}

	//统计正负残差点并写入节点信息
	int positive, negative, total;
	positive = 0;
	negative = 0;
	double thresh = 0.7;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i)->residue > thresh)
		{
			positive++;
		}
		if ((tri + i)->residue < -thresh)
		{
			negative++;
		}
	}
	bool b_balanced = (positive == negative);
	if (negative == 0 || positive == 0)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "write_DIMACS(): no residue point!\n\n");
		return -1;
	}
	fprintf(fp, "c This is MCF problem file.\n");
	fprintf(fp, "c Problem line(nodes, links)\n");
	//统计边缘三角形个数
	long boundry_tri = 0;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i) != NULL)
		{
			if ((edges + (tri + i)->edge1 - 1)->isBoundry ||
				(edges + (tri + i)->edge2 - 1)->isBoundry ||
				(edges + (tri + i)->edge3 - 1)->isBoundry)
			{
				boundry_tri++;
			}
		}
	}
	long n;
	if (!b_balanced)
	{
		n = num_triangle + 1;
		fprintf(fp, "p min %ld %ld\n", n, num_arcs + boundry_tri * 2);
	}
	else
	{
		n = num_triangle;
		fprintf(fp, "p min %ld %ld\n", n, num_arcs);
	}
	fprintf(fp, "c Node descriptor lines\n");
	positive = 0;
	negative = 0;
	int count = 0;
	bool b_positive, b_negative, is_residue;
	double sum = 0.0;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i) != NULL && (tri + i)->residue > thresh)
		{
			fprintf(fp, "n %d %lf\n", i + 1, (tri + i)->residue);
			sum += (tri + i)->residue;
		}
		if ((tri + i) != NULL && (tri + i)->residue < -thresh)
		{
			fprintf(fp, "n %d %lf\n", i + 1, (tri + i)->residue);
			sum += (tri + i)->residue;
		}
	}
	//写入大地节点
	if (!b_balanced)
	{
		fprintf(fp, "n %d %lf\n", num_triangle + 1, -sum);
	}

	//写入流费用
	fprintf(fp, "c Arc descriptor lines(from, to, minflow, maxflow, cost)\n");
	int rows, cols;
	int lower_bound = 0;
	int upper_bound = 5;
	double cost_mean;
	int nr = cost.rows;
	int nc = cost.cols;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i) != NULL &&
			(tri + i)->p1 >= 1 &&
			(tri + i)->p1 <= num_nodes &&
			(tri + i)->p2 >= 1 &&
			(tri + i)->p2 <= num_nodes &&
			(tri + i)->p3 >= 1 &&
			(tri + i)->p3 <= num_nodes
			)
		{
			cost_mean = 0.0;
			nodes[(tri + i)->p1 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			nodes[(tri + i)->p2 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			nodes[(tri + i)->p3 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			cost_mean = cost_mean / 3;
			if ((tri + i)->neigh1 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, (tri + i)->neigh1, lower_bound, upper_bound, cost_mean);
			if ((tri + i)->neigh2 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, (tri + i)->neigh2, lower_bound, upper_bound, cost_mean);
			if ((tri + i)->neigh3 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, (tri + i)->neigh3, lower_bound, upper_bound, cost_mean);
		}
	}
	if (!b_balanced)
	{
		//写入边界流费用
		for (int i = 0; i < num_triangle; i++)
		{
			if ((tri + i) != NULL &&
				(((edges + (tri + i)->edge1 - 1)->isBoundry) ||
					((edges + (tri + i)->edge2 - 1)->isBoundry) ||
					((edges + (tri + i)->edge3 - 1)->isBoundry)
					))
			{
				cost_mean = 0.0;
				nodes[(tri + i)->p1 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				nodes[(tri + i)->p2 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				nodes[(tri + i)->p3 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				cost_mean = cost_mean / 3;
				fprintf(fp, "a %d %d %d %d %lf\n", i + 1, num_triangle + 1, lower_bound, upper_bound, cost_mean);
				fprintf(fp, "a %d %d %d %d %lf\n", num_triangle + 1, i + 1, lower_bound, upper_bound, cost_mean);
			}
		}
	}
	if (fp) fclose(fp);
	fp = NULL;
	return 0;
}

int Utils::write_DIMACS(
	const char* DIMACS_file_problem,
	vector<triangle>& triangle,
	vector<tri_node>& nodes,
	vector<tri_edge>& edges,
	const Mat& cost
)
{
	if (DIMACS_file_problem == NULL ||
		triangle.size() < 1||
		nodes.size() < 3 ||
		edges.size() < 3 ||
		cost.rows < 2 ||
		cost.cols < 2 ||
		cost.channels() != 1 ||
		cost.type() != CV_64F
		)
	{
		fprintf(stderr, "write_DIMACS(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(DIMACS_file_problem, "wt");
	if (fp == NULL)
	{
		fprintf(stderr, "write_DIMACS(): can't open %s\n", DIMACS_file_problem);
		return -1;
	}

	int ret, num_nodes;
	int num_triangle = triangle.size();
	num_nodes = nodes.size();
	long num_arcs = 0;
	for (int i = 0; i < num_triangle; i++)
	{
		if (triangle[i].neigh1 > 0) num_arcs++;
		if (triangle[i].neigh2 > 0) num_arcs++;
		if (triangle[i].neigh3 > 0) num_arcs++;
	}

	//统计正负残差点并写入节点信息
	int positive, negative, total;
	positive = 0;
	negative = 0;
	double thresh = 0.7;
	for (int i = 0; i < num_triangle; i++)
	{
		if (triangle[i].residue > thresh)
		{
			positive++;
		}
		if (triangle[i].residue < -thresh)
		{
			negative++;
		}
	}
	bool b_balanced = (positive == negative);
	if (negative == 0 && positive == 0)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "write_DIMACS(): no residue point!\n\n");
		return -1;
	}
	fprintf(fp, "c This is MCF problem file.\n");
	fprintf(fp, "c Problem line(nodes, links)\n");
	//统计边缘三角形个数
	long boundry_tri = 0;
	for (int i = 0; i < num_triangle; i++)
	{
		if (edges[triangle[i].edge1 - 1].isBoundry ||
			edges[triangle[i].edge2 - 1].isBoundry ||
			edges[triangle[i].edge3 - 1].isBoundry)
		{
			boundry_tri++;
		}
	}
	long n;
	if (!b_balanced)
	{
		n = num_triangle + 1;
		fprintf(fp, "p min %ld %ld\n", n, num_arcs + boundry_tri * 2);
	}
	else
	{
		n = num_triangle;
		fprintf(fp, "p min %ld %ld\n", n, num_arcs);
	}
	fprintf(fp, "c Node descriptor lines\n");
	positive = 0;
	negative = 0;
	int count = 0;
	bool b_positive, b_negative, is_residue;
	double sum = 0.0;
	for (int i = 0; i < num_triangle; i++)
	{
		if (triangle[i].residue > thresh)
		{
			fprintf(fp, "n %d %lf\n", i + 1, triangle[i].residue);
			sum += triangle[i].residue;
		}
		if (triangle[i].residue < -thresh)
		{
			fprintf(fp, "n %d %lf\n", i + 1, triangle[i].residue);
			sum += triangle[i].residue;
		}
	}
	//写入大地节点
	if (!b_balanced)
	{
		fprintf(fp, "n %d %lf\n", num_triangle + 1, -sum);
	}

	//写入流费用
	fprintf(fp, "c Arc descriptor lines(from, to, minflow, maxflow, cost)\n");
	int rows, cols;
	int lower_bound = 0;
	int upper_bound = 1;
	double cost_mean;
	int nr = cost.rows;
	int nc = cost.cols;
	for (int i = 0; i < num_triangle; i++)
	{
		if (triangle[i].p1 >= 1 &&
			triangle[i].p1 <= num_nodes &&
			triangle[i].p2 >= 1 &&
			triangle[i].p2 <= num_nodes &&
			triangle[i].p3 >= 1 &&
			triangle[i].p3 <= num_nodes
			)
		{
			cost_mean = 0.0;
			nodes[triangle[i].p1 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			nodes[triangle[i].p2 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			nodes[triangle[i].p3 - 1].get_pos(&rows, &cols);
			if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
			cost_mean = cost_mean / 3;
			if (triangle[i].neigh1 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, triangle[i].neigh1, lower_bound, upper_bound, cost_mean);
			if (triangle[i].neigh2 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, triangle[i].neigh2, lower_bound, upper_bound, cost_mean);
			if (triangle[i].neigh3 > 0) fprintf(fp, "a %d %d %d %d %lf\n", i + 1, triangle[i].neigh3, lower_bound, upper_bound, cost_mean);
		}
	}
	if (!b_balanced)
	{
		//写入边界流费用
		for (int i = 0; i < num_triangle; i++)
		{
			if (edges[triangle[i].edge1 - 1].isBoundry ||
				edges[triangle[i].edge2 - 1].isBoundry ||
				edges[triangle[i].edge3 - 1].isBoundry)
			{
				cost_mean = 0.0;
				nodes[triangle[i].p1 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				nodes[triangle[i].p2 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				nodes[triangle[i].p3 - 1].get_pos(&rows, &cols);
				if (rows >= 0 && rows <= nr - 1) cost_mean += cost.at<double>(rows, cols);
				cost_mean = cost_mean / 3;
				fprintf(fp, "a %d %d %d %d %lf\n", i + 1, num_triangle + 1, lower_bound, upper_bound, cost_mean);
				fprintf(fp, "a %d %d %d %d %lf\n", num_triangle + 1, i + 1, lower_bound, upper_bound, cost_mean);
			}
		}
	}
	if (fp) fclose(fp);
	fp = NULL;
	return 0;
}

int Utils::read_DIMACS(const char* DIMACS_file_solution, Mat& k1, Mat& k2, int rows, int cols)
{
	if (rows < 2 || cols < 2)
	{
		fprintf(stderr, "read_DIMACS(): input check failed!\n\n");
		return -1;
	}
	char instring[256];
	char ch;
	double obj_value = 0;
	int i, row_index, col_index, symbol;
	long from, to;
	double flow = 0;
	k1 = Mat::zeros(rows - 1, cols, CV_64F);
	k2 = Mat::zeros(rows, cols - 1, CV_64F);
	long earth_node_indx = (rows - 1) * (cols - 1) + 1;
	FILE* fp = NULL;
	fopen_s(&fp, DIMACS_file_solution, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_DIMACS(): can't open file %s\n\n", DIMACS_file_solution);
		return -1;
	}
	/////////////////////读取注释///////////////////////////
	GET_NEXT_LINE;
	while (ch != 's' && ch)
	{
		if (ch != 'c')
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}
		GET_NEXT_LINE;
	}
	/////////////////////读取优化目标值/////////////////////
	for (i = 1; i < 81; i++)
	{
		if (isspace((int)instring[i]) > 0)
		{
			i++;
			break;
		}
	}
	if (sscanf(&(instring[i]), "%lf", &obj_value) != 1)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
		return -1;
	}
	if (obj_value < 0.0)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): this problem can't be solved(unbounded or infeasible)!\n\n");
		return -1;
	}
	////////////////////读取MCF结果/////////////////////////
	GET_NEXT_LINE;
	while (ch && ch == 'f')
	{
		if (sscanf(&(instring[2]), "%ld %ld %lf", &from, &to, &flow) != 3 ||
			flow < 0.0 || from < 0 || to < 0)
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}
		//是否为接地弧
		if (from == earth_node_indx || to == earth_node_indx)
		{
			if (from == earth_node_indx)
			{
				//top
				if (to <= cols - 1)
				{
					row_index = 0;
					col_index = to - 1;
					symbol = 1;
					k2.at<double>(row_index, col_index) = k2.at<double>(row_index, col_index) + symbol * flow;
				}
				//bottom
				else if (to > (rows - 2) * (cols - 1))
				{
					symbol = -1;
					row_index = rows - 1;
					col_index = to - (rows - 2) * (cols - 1) - 1;
					k2.at<double>(row_index, col_index) = k2.at<double>(row_index, col_index) + symbol * flow;
				}
				//right
				else if (to % (cols - 1) == 0 && to < (rows - 2) * (cols - 1) && to > (cols - 1))
				{
					symbol = 1;
					row_index = to / (cols - 1) - 1;
					col_index = cols - 1;
					k1.at<double>(row_index, col_index) = k1.at<double>(row_index, col_index) + symbol * flow;
				}
				//left
				else
				{
					symbol = -1;
					row_index = to / (cols - 1);
					col_index = 0;
					k1.at<double>(row_index, col_index) = k1.at<double>(row_index, col_index) + symbol * flow;
				}
			}
			else
			{
				/*long tmp;
				tmp = from;
				from = to;
				to = tmp;*/
				//top
				if (from <= cols - 1)
				{
					symbol = -1;
					row_index = 0;
					col_index = from - 1;
					k2.at<double>(row_index, col_index) = k2.at<double>(row_index, col_index) + symbol * flow;
				}
				//bottom
				else if (from >= (rows - 2) * (cols - 1))
				{
					symbol = 1;
					row_index = rows - 1;
					col_index = from - (rows - 2) * (cols - 1) - 1;
					k2.at<double>(row_index, col_index) = k2.at<double>(row_index, col_index) + symbol * flow;
				}
				//right
				else if (from % (cols - 1) == 0 && from < (rows - 2) * (cols - 1) && from >(cols - 1))
				{
					symbol = -1;
					row_index = from / (cols - 1) - 1;
					col_index = cols - 1;
					k1.at<double>(row_index, col_index) = k1.at<double>(row_index, col_index) + symbol * flow;
				}
				//left
				else
				{
					symbol = 1;
					row_index = from / (cols - 1);
					col_index = 0;
					k1.at<double>(row_index, col_index) = k1.at<double>(row_index, col_index) + symbol * flow;
				}
			}
		}
		else
		{
			if (abs(from - to) > 1)
			{
				symbol = from > to ? 1 : -1;
				row_index = (from > to ? to : from) / long(cols - 1);
				col_index = (from > to ? to : from) % long(cols - 1);
				if (col_index == 0)
				{
					col_index = cols - 1;
					row_index--;
				}
				col_index--;
				k2.at<double>(row_index + 1, col_index) = k2.at<double>(row_index + 1, col_index) + symbol * flow;
			}
			else
			{
				symbol = from > to ? 1 : -1;
				row_index = (from > to ? to : from) / long(cols - 1);
				col_index = (from > to ? to : from) % long(cols - 1);
				k1.at<double>(row_index, col_index) = k1.at<double>(row_index, col_index) + symbol * flow;
			}
		}
		GET_NEXT_LINE;
	}
	if (ch != 'c')
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
		return -1;
	}
	if (fp)
	{
		fclose(fp);
	}
	return 0;
}

int Utils::write_DIMACS(const char* DIMACS_file_problem, Mat& residue, Mat& coherence, double thresh)
{

	if (residue.cols < 2 ||
		residue.rows < 2 ||
		coherence.cols < 2 ||
		coherence.rows < 2||
		residue.type() != CV_64F||
		coherence.type() != CV_64F||
		(coherence.rows - residue.rows) != 1||
		(coherence.cols - residue.cols) != 1||
		thresh < 0.0)
	{
		fprintf(stderr, "write_DIMACS(): input check failed!\n\n");
		return -1;
	}
	long nr = residue.rows;
	long nc = residue.cols;
	long i, j;
	long node_index = 1;
	double sum = 0.0;
	//统计正负残差点数
	long positive, negative, total, Arcs_num, Nodes_num;
	positive = 0;
	negative = 0;
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++)
		{
			if (residue.at<double>(i, j) > thresh)
			{
				positive++;
			}
			if (residue.at<double>(i, j) < -thresh)
			{
				negative++;
			}
		}
	}
	bool b_balanced = (positive == negative);
	if (!b_balanced)
	{
		Nodes_num = residue.rows * residue.cols + 1;
		Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1) +
			2 * 2 * residue.cols + 2 * 2 * (residue.rows - 2);
	}
	else
	{
		Nodes_num = residue.rows * residue.cols;
		Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1);
	}
	ofstream fout;
	FILE* fp = NULL;
	fopen_s(&fp, DIMACS_file_problem, "wt");
	if (!fp)
	{
		fprintf(stderr, "write_DIMACS(): cant't open file %s\n\n", DIMACS_file_problem);
		return -1;
	}
	fprintf(fp, "c This is a DIMACS file, describing Minimum Cost Flow problem.\n");
	fprintf(fp, "c Problem line (nodes, links)\n");
	fprintf(fp, "p min %ld %ld\n", Nodes_num, Arcs_num);
	fprintf(fp, "c Node descriptor lines (supply+ or demand-)\n");

	
	/*
	* 写入节点的度（残差值1，-1）
	*/
	
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++)
		{
			if (residue.at<double>(i, j) > thresh)
			{
				node_index = i * nc + j + 1;
				fprintf(fp, "n %ld %lf\n", node_index, residue.at<double>(i, j));
				sum += residue.at<double>(i, j);
			}
			if (residue.at<double>(i, j) < -thresh)
			{
				node_index = i * nc + j + 1;
				fprintf(fp, "n %ld %lf\n", node_index, residue.at<double>(i, j));
				sum += residue.at<double>(i, j);
			}
		}
	}

	/*写接地节点*/
	
	node_index = nc * nr + 1;
	if (!b_balanced)
	{
		fprintf(fp, "n %ld %lf\n", node_index, -sum);
	}
	

	long earth_node_index = node_index;

	/*
	* 写入每个有向弧的费用（流费用）
	*/
	long lower_bound = 0;
	long upper_bound = 5;
	double mean_coherence1, mean_coherence2, mean_coherence3, mean_coherence4;
	fprintf(fp, "c Arc descriptor lines (from, to, minflow, maxflow, cost)\n");
	/*接地节点的有向弧流费用*/
	if (!b_balanced)
	{
		//top
		for (i = 0; i < nc; i++)
		{
			node_index = i + 1;
			mean_coherence1 = coherence.at<double>(0, i);
			mean_coherence2 = mean_coherence1;
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_coherence1,
				earth_node_index, node_index, lower_bound, upper_bound, mean_coherence2);
		}
		//bottom
		for (i = 0; i < nc; i++)
		{
			node_index = nc * (nr - 1) + i + 1;
			mean_coherence1 = coherence.at<double>(nr - 1, i);
			mean_coherence2 = mean_coherence1;
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_coherence1,
				earth_node_index, node_index, lower_bound, upper_bound, mean_coherence2);
		}
		//left
		for (i = 1; i < nr - 1; i++)
		{
			node_index = nc * i + 1;
			mean_coherence1 = coherence.at<double>(i, 0);
			mean_coherence2 = mean_coherence1;
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_coherence1,
				earth_node_index, node_index, lower_bound, upper_bound, mean_coherence2);
		}
		//right
		for (i = 1; i < nr - 1; i++)
		{
			node_index = nc * (i + 1);
			mean_coherence1 = coherence.at<double>(i, nc - 1);
			mean_coherence2 = mean_coherence1;
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_coherence1,
				earth_node_index, node_index, lower_bound, upper_bound, mean_coherence2);
		}
	}

	/*非接地节点的有向弧流费用*/
	for (i = 0; i < nr - 1; i++)
	{
		for (j = 0; j < nc - 1; j++)
		{
			node_index = i * nc + j + 1;
			/*正向*/
			mean_coherence1 = mean(coherence(Range(i, i + 2), Range(j, j + 3))).val[0];
			/*逆向*/
			mean_coherence2 = mean(coherence(Range(i, i + 2), Range(j, j + 3))).val[0];


			/*正向*/
			mean_coherence3 = mean(coherence(Range(i, i + 3), Range(j, j + 2))).val[0];
			/*逆向*/
			mean_coherence4 = mean(coherence(Range(i, i + 3), Range(j, j + 2))).val[0];
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, node_index + 1, lower_bound, upper_bound, mean_coherence1,
				node_index + 1, node_index, lower_bound, upper_bound, mean_coherence2,
				node_index, node_index + nc, lower_bound, upper_bound, mean_coherence3,
				node_index + nc, node_index, lower_bound, upper_bound, mean_coherence4);

		}
	}



	for (j = 0; j < nc - 1; j++)
	{
		node_index = (nr - 1) * nc + j + 1;
		/*正向*/
		mean_coherence1 = mean(coherence(Range(nr - 1, nr + 1), Range(j, j + 3))).val[0];
		/*逆向*/
		mean_coherence2 = mean(coherence(Range(nr - 1, nr + 1), Range(j, j + 3))).val[0];
		fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
			node_index, node_index + 1, lower_bound, upper_bound, mean_coherence1,
			node_index + 1, node_index, lower_bound, upper_bound, mean_coherence2);
	}

	for (i = 0; i < nr - 1; i++)
	{
		node_index = (i + 1) * nc;
		/*正向*/
		mean_coherence1 = mean(coherence(Range(i, i + 3), Range(nc - 1, nc + 1))).val[0];
		/*逆向*/
		mean_coherence2 = mean(coherence(Range(i, i + 3), Range(nc - 1, nc + 1))).val[0];
		fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
			node_index, node_index + nc, lower_bound, upper_bound, mean_coherence1,
			node_index + nc, node_index, lower_bound, upper_bound, mean_coherence2);
	}
	fprintf(fp, "c ");
	fprintf(fp, "c End of file");
	if (fp) fclose(fp);
	fp = NULL;
	return 0;
}

int Utils::cumsum(Mat& phase, int dim)
{
	/*
	cumulates along the dimension specified by dim
	dim = 1,按列计算
	dim = 2,按行计算
	*/
	int rows = phase.rows;
	int cols = phase.cols;
	int i, j;
	if (rows > 0 && cols > 0)
	{
		switch (dim)
		{
		case 1:
			if (phase.rows < 2 ||
				phase.cols < 1 ||
				phase.type() != CV_64F)
			{
				fprintf(stderr, "cumsum(): input check failed!\n\n");
				return -1;
			}
			for (j = 0; j < cols; j++)
			{
				for (i = 1; i < rows; i++)
				{
					phase.at<double>(i, j) = phase.at<double>(i - 1, j) + phase.at<double>(i, j);
				}
			}
			break;
		case 2:
			if (phase.rows < 1 ||
				phase.cols < 2 ||
				phase.type() != CV_64F)
			{
				fprintf(stderr, "cumsum(): input check failed!\n\n");
				return -1;
			}
			for (i = 0; i < rows; i++)
			{
				for (j = 1; j < cols; j++)
				{
					phase.at<double>(i, j) = phase.at<double>(i, j - 1) + phase.at<double>(i, j);
				}
			}
			break;
		default:
			break;
		}
	}
	return 0;
}

int Utils::cross(Mat& vec1, Mat& vec2, Mat& out)
{
	if (vec1.cols != 3 ||
		vec1.rows < 1 ||
		vec1.type() != CV_64F ||
		vec1.channels() != 1 ||
		vec1.cols != vec2.cols ||
		vec1.rows != vec2.rows ||
		vec2.channels() != 1 ||
		vec2.type() != CV_64F
		)
	{
		fprintf(stderr, "cross(): input check failed!\n\n");
		return -1;
	}
	Mat out_tmp = Mat::zeros(vec1.rows, vec1.cols, CV_64F);
	int rows = vec1.rows;
	for (int i = 0; i < rows; i++)
	{
		out_tmp.at<double>(i, 0) = vec1.at<double>(i, 1) * vec2.at<double>(i, 2) -
			vec1.at<double>(i, 2) * vec2.at<double>(i, 1);//a(2) * b(3) - a(3) * b(2)

		out_tmp.at<double>(i, 1) = vec1.at<double>(i, 2) * vec2.at<double>(i, 0) -
			vec1.at<double>(i, 0) * vec2.at<double>(i, 2);//a(3) * b(1) - a(1) * b(3)

		out_tmp.at<double>(i, 2) = vec1.at<double>(i, 0) * vec2.at<double>(i, 1) -
			vec1.at<double>(i, 1) * vec2.at<double>(i, 0);//a(1) * b(2) - a(2) * b(1)
	}
	out_tmp.copyTo(out);
	return 0;
}

int Utils::gen_mask(Mat& coherence, Mat& phase_derivatives, Mat& mask, int wnd_size, double coh_thresh, double phase_derivative_thresh)
{
	if (coherence.rows < 2 ||
		coherence.cols < 2 ||
		coherence.channels() != 1 ||
		coherence.type() != CV_64F ||
		coherence.rows != phase_derivatives.rows ||
		coherence.cols  != phase_derivatives.cols ||
		phase_derivatives.channels() != 1 ||
		phase_derivatives.type() != CV_64F ||
		wnd_size < 0 ||
		wnd_size > coherence.rows ||
		coh_thresh < 0.0 ||
		coh_thresh > 1.0||
		phase_derivative_thresh < 0.0
		)
	{
		fprintf(stderr, "gen_mask(): input check failed!\n\n");
		return -1;
	}
	int nr = coherence.rows;
	int nc = coherence.cols;
	Mat temp_coh, temp_phase_derivatives;
	coherence.copyTo(temp_coh);
	phase_derivatives.copyTo(temp_phase_derivatives);
	int radius = (wnd_size - 1) / 2;
	cv::copyMakeBorder(temp_coh, temp_coh, radius, radius, radius, radius, cv::BORDER_DEFAULT);
	cv::copyMakeBorder(temp_phase_derivatives, temp_phase_derivatives, radius, radius, radius, radius, cv::BORDER_DEFAULT);
	Mat tmp = Mat::zeros(nr, nc, CV_32S);
	tmp.copyTo(mask);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double mean, mean1;
		for (int j = 0; j < nc; j++)
		{
			mean = cv::mean(temp_coh(cv::Range(i, i + 2 * radius + 1), cv::Range(j, j + 2 * radius + 1)))[0];
			mean1 = cv::mean(temp_phase_derivatives(cv::Range(i, i + 2 * radius + 1), cv::Range(j, j + 2 * radius + 1)))[0];
			if (mean > coh_thresh && 
				coherence.at<double>(i, j) > coh_thresh&&
				mean1 < phase_derivative_thresh &&
				phase_derivatives.at<double>(i, j) < phase_derivative_thresh)
			{
				mask.at<int>(i, j) = 1;
			}
				
		}
	}
	return 0;
}

int Utils::residue_sift(Mat& residue_src, Mat& residue_dst, double thresh, long* num_residue)
{
	int rows = residue_src.rows;
	int cols = residue_src.cols;
	if ( rows < 1 ||
		 cols < 1 ||
		 residue_src.type() != CV_64F||
		 residue_src.channels() != 1||
		 thresh < 0.0 ||
		 num_residue == NULL
		)
	{
		fprintf(stderr, "residue_sift(): input check failed!\n\n");
		return -1;
	}
	Mat tmp;
	residue_src.copyTo(tmp);
	*num_residue = 0;
//#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (std::fabs(residue_src.at<double>(i, j)) > thresh)
			{
				tmp.at<double>(i, j) = 1.0;
				*num_residue = *num_residue + 1;
			}
			else
			{
				tmp.at<double>(i, j) = 0.0;
			}
		}
	}
	tmp.copyTo(residue_dst);
	return 0;
}

int Utils::wrap(Mat& Src, Mat& Dst)
{
	int rows = Src.rows;
	int cols = Src.cols;
	double pi = 3.1415926535;
	if (rows < 1 || cols < 1 || Src.type() != CV_64F)
	{
		fprintf(stderr, "wrap(): input check failed!\n\n");
		return -1;
	}
	Mat tmp = Mat::zeros(rows, cols, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			tmp.at<double>(i, j) = atan2(sin(Src.at<double>(i, j)), cos(Src.at<double>(i, j)));
		}
	}
	//Dst = tmp;
	tmp.copyTo(Dst);
	return 0;
}

int Utils::residue(Mat& phase, Mat& residuemat)
{
	int rows = phase.rows;
	int cols = phase.cols;
	int ret;
	if (rows < 2 || cols < 2 || phase.type() != CV_64F)
	{
		fprintf(stderr, "residue(): input check failed!\n\n");
		return -1;
	}
	Mat Diff_1 = phase(Range(1, rows), Range(0, cols)) - phase(Range(0, rows - 1), Range(0, cols));
	Mat Diff_2 = phase(Range(0, rows), Range(1, cols)) - phase(Range(0, rows), Range(0, cols - 1));
	ret = this->wrap(Diff_1, Diff_1);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;
	ret = this->wrap(Diff_2, Diff_2);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;

	Diff_1 = Diff_1(Range(0, Diff_1.rows), Range(1, Diff_1.cols)) - 
		Diff_1(Range(0, Diff_1.rows), Range(0, Diff_1.cols - 1));

	Diff_2 = Diff_2(Range(1, Diff_2.rows), Range(0, Diff_2.cols)) - 
		Diff_2(Range(0, Diff_2.rows - 1), Range(0, Diff_2.cols));

	Diff_1 = Diff_2 - Diff_1;
	double pi = 3.1415926535;
	Diff_1 = Diff_1 / (2 * pi);
	//Diff_1.copyTo(residuemat);
	residuemat = Diff_1;
	return 0;
}

int Utils::residue(triangle* tri, int num_triangle, vector<tri_node>& nodes, tri_edge* edges, int num_edges)
{
	if (tri == NULL ||
		num_triangle < 1 ||
		nodes.size() < 1||
		edges == NULL||
		num_edges < 3
		)
	{
		fprintf(stderr, "residue(): input check failed!\n\n");
		return -1;
	}
	double thresh = 50.0;
	int num_nodes = nodes.size();
	int end1, end2, end3, tmp;
	double x1, y1, x2, y2, x3, y3, direction, delta12, delta23, delta31, residue, phi1, phi2, phi3, distance1,
		distance2, distance3;
	int row1, col1, row2, col2, row3, col3;
	bool b_res = false;
	for (int i = 0; i < num_triangle; i++)
	{
		if ((tri + i) != NULL)
		{
			end1 = (tri + i)->p1;
			end2 = (tri + i)->p2;
			end3 = (tri + i)->p3;
		}
		if (end1 > end2)
		{
			tmp = end1;
			end1 = end2;
			end2 = tmp;
		}
		
		nodes[end1 - 1].get_pos(&row1, &col1);
		nodes[end2 - 1].get_pos(&row2, &col2);
		nodes[end3 - 1].get_pos(&row3, &col3);

		nodes[end1 - 1].get_distance(nodes[end2 - 1], &distance1);
		nodes[end2 - 1].get_distance(nodes[end3 - 1], &distance2);
		nodes[end3 - 1].get_distance(nodes[end1 - 1], &distance3);
		b_res = true;
		if ((distance1 > thresh) || (distance2 > thresh) || (distance3 > thresh)) b_res = false;

		nodes[end1 - 1].get_phase(&phi1);
		nodes[end2 - 1].get_phase(&phi2);
		nodes[end3 - 1].get_phase(&phi3);

		x2 = double(col2 - col1);
		y2 = double(row1 - row2);
		x1 = double(col1 - col3);
		y1 = double(row3 - row1);
		direction = x1 * y2 - x2 * y1;

		delta12 = atan2(sin(phi2 - phi1), cos(phi2 - phi1));
		delta23 = atan2(sin(phi3 - phi2), cos(phi3 - phi2));
		delta31 = atan2(sin(phi1 - phi3), cos(phi1 - phi3));

		double res = (delta12 + delta23 + delta31) / 2.0 / PI;
		if (fabs(res) > 0.7 && !b_res)//标注边长超过阈值的残差边和残差节点
		{
			(edges + (tri + i)->edge1 - 1)->isResidueEdge = true;
			(edges + (tri + i)->edge2 - 1)->isResidueEdge = true;
			(edges + (tri + i)->edge3 - 1)->isResidueEdge = true;
			nodes[(tri + i)->p1 - 1].set_residue(true);
			nodes[(tri + i)->p2 - 1].set_residue(true);
			nodes[(tri + i)->p3 - 1].set_residue(true);
		}
		res = b_res ? res : 0.0;
		if (direction > 0.0)//在目标三角形中顺残差方向(残差方向定义为逆时针方向)
		{
			(tri + i)->residue = res;
		}
		else
		{
			(tri + i)->residue = -res;
		}
	}
	return 0;
}

int Utils::residue(vector<triangle>& triangle, vector<tri_node>& nodes, vector<tri_edge>& edges, double distance_thresh)
{
	if (triangle.size() < 1 ||
		nodes.size() < 1 ||
		edges.size() < 3
		)
	{
		fprintf(stderr, "residue(): input check failed!\n\n");
		return -1;
	}
	int num_triangle = triangle.size();
	double thresh;
	thresh = distance_thresh < 2.0 ? 2.0 : distance_thresh;
	int num_nodes = nodes.size();
	int end1, end2, end3, tmp;
	double x1, y1, x2, y2, x3, y3, direction, delta12, delta23, delta31, residue, phi1, phi2, phi3, distance1,
		distance2, distance3;
	int row1, col1, row2, col2, row3, col3;
	bool b_res = false;
	for (int i = 0; i < num_triangle; i++)
	{
		end1 = triangle[i].p1;
		end2 = triangle[i].p2;
		end3 = triangle[i].p3;
		if (end1 > end2)
		{
			tmp = end1;
			end1 = end2;
			end2 = tmp;
		}

		nodes[end1 - 1].get_pos(&row1, &col1);
		nodes[end2 - 1].get_pos(&row2, &col2);
		nodes[end3 - 1].get_pos(&row3, &col3);

		nodes[end1 - 1].get_distance(nodes[end2 - 1], &distance1);
		nodes[end2 - 1].get_distance(nodes[end3 - 1], &distance2);
		nodes[end3 - 1].get_distance(nodes[end1 - 1], &distance3);
		b_res = true;
		if ((distance1 > thresh) || (distance2 > thresh) || (distance3 > thresh)) b_res = false;

		nodes[end1 - 1].get_phase(&phi1);
		nodes[end2 - 1].get_phase(&phi2);
		nodes[end3 - 1].get_phase(&phi3);

		x2 = double(col2 - col1);
		y2 = -double(row1 - row2);
		x1 = double(col1 - col3);
		y1 = -double(row3 - row1);
		direction = x1 * y2 - x2 * y1;

		delta12 = atan2(sin(phi2 - phi1), cos(phi2 - phi1));
		delta23 = atan2(sin(phi3 - phi2), cos(phi3 - phi2));
		delta31 = atan2(sin(phi1 - phi3), cos(phi1 - phi3));

		double res = (delta12 + delta23 + delta31) / 2.0 / PI;
		if (fabs(res) > 0.7 && !b_res)//标注边长超过阈值的残差边和残差节点
		{
			edges[triangle[i].edge1 - 1].isResidueEdge = true;
			edges[triangle[i].edge2 - 1].isResidueEdge = true;
			edges[triangle[i].edge3 - 1].isResidueEdge = true;

			nodes[triangle[i].p1 - 1].set_residue(true);
			nodes[triangle[i].p2 - 1].set_residue(true);
			nodes[triangle[i].p3 - 1].set_residue(true);
		}
		res = b_res ? res : 0.0;
		if (direction < 0.0)//在目标三角形中顺残差方向(残差方向定义为逆时针方向)
		{
			triangle[i].residue = res;
		}
		else
		{
			triangle[i].residue = -res;
		}
	}
	return 0;
}

int Utils::gen_mask(Mat& coherence, Mat& mask, int wnd_size, double thresh)
{
	if (coherence.rows < 2 ||
		coherence.cols < 2 ||
		coherence.channels() != 1 ||
		coherence.type() != CV_64F ||
		wnd_size < 0 ||
		wnd_size > coherence.rows ||
		thresh < 0.0 ||
		thresh > 1.0
		)
	{
		fprintf(stderr, "gen_mask(): input check failed!\n\n");
		return -1;
	}
	int nr = coherence.rows;
	int nc = coherence.cols;
	Mat temp_coh;
	coherence.copyTo(temp_coh);
	int radius = (wnd_size + 1) / 2;
	cv::copyMakeBorder(temp_coh, temp_coh, radius, radius, radius, radius, cv::BORDER_REFLECT);
	Mat tmp = Mat::zeros(nr, nc, CV_32S);
	tmp.copyTo(mask);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double mean;
		for (int j = 0; j < nc; j++)
		{
			mean = cv::mean(temp_coh(cv::Range(i, i + 2 * radius + 1), cv::Range(j, j + 2 * radius + 1)))[0];
			if (mean > thresh && coherence.at<double>(i, j) > thresh) mask.at<int>(i, j) = 1;
		}
	}
	return 0;
}

int Utils::real_coherence(ComplexMat& Mast, ComplexMat& Slave, Mat& coherence)
{
	int wa = 3;  //窗口方位向尺寸
	int wr = 3;  //窗口距离向尺寸

	int na = Mast.GetRows();
	int nr = Mast.GetCols();
	if ((na < 3) ||
		(nr < 3) ||
		Mast.re.type() != CV_64F ||
		Slave.re.type() != CV_64F ||
		Mast.GetCols() != Slave.GetCols()||
		Mast.GetRows() != Slave.GetRows())
	{
		fprintf(stderr, "real_coherence(): input check failed!\n\n");
		return -1;
	}

	int win_a = (wa - 1) / 2; //方位窗半径
	int win_r = (wr - 1) / 2; //距离窗半径

	int na_new = na - 2 * win_a;
	int nr_new = nr - 2 * win_r;

	Mat Coherence(na_new, nr_new, CV_64F, Scalar::all(0));


#pragma omp parallel for schedule(guided)
	for (int i = win_a + 1; i <= na - win_a; i++)
	{
		for (int j = win_r + 1; j <= nr - win_r; j++)
		{
			Mat s1, s2, sum1, sum2;

			double up, down;
			magnitude(Mast.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), Mast.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), s1);
			magnitude(Slave.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), Slave.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), s2);
			up = sum((s1.mul(s1)).mul(s2.mul(s2)))[0];
			pow(s1, 4, s1);
			pow(s2, 4, s2);
			down = sqrt(sum(s1)[0] * sum(s2)[0]);
			if (up / (down + 1e-12) > 1.0)
			{
				Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = 1;
			}
			else
			{
				Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = up / (down + 1e-12);
			}

		}
	}
	copyMakeBorder(Coherence, Coherence, 1, 1, 1, 1, BORDER_REFLECT);
	coherence = Coherence;
	return 0;
}

int Utils::real_coherence(const ComplexMat& master_image, const ComplexMat& slave_image, int est_wndsize_rg, int est_wndsize_az, Mat& coherence)
{


	int na = master_image.GetRows();
	int nr = master_image.GetCols();
	if ((na < est_wndsize_az) ||
		(nr < est_wndsize_rg) ||
		master_image.type() != CV_64F ||
		slave_image.type() != CV_64F ||
		master_image.GetCols() != slave_image.GetCols() ||
		master_image.GetRows() != slave_image.GetRows() ||
		est_wndsize_rg % 2 == 0||
		est_wndsize_az % 2 == 0 ||
		est_wndsize_rg < 3 ||
		est_wndsize_az < 3
		)
	{
		fprintf(stderr, "real_coherence(): input check failed!\n\n");
		return -1;
	}

	int win_a = (est_wndsize_az - 1) / 2; //方位窗半径
	int win_r = (est_wndsize_rg - 1) / 2; //距离窗半径

	int na_new = na - 2 * win_a;
	int nr_new = nr - 2 * win_r;

	Mat Coherence(na_new, nr_new, CV_64F, Scalar::all(0));


#pragma omp parallel for schedule(guided)
	for (int i = win_a + 1; i <= na - win_a; i++)
	{
		for (int j = win_r + 1; j <= nr - win_r; j++)
		{
			Mat s1, s2, sum1, sum2;

			double up, down;
			magnitude(master_image.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), master_image.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), s1);
			magnitude(slave_image.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), slave_image.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)), s2);
			up = sum((s1.mul(s1)).mul(s2.mul(s2)))[0];
			pow(s1, 4, s1);
			pow(s2, 4, s2);
			down = sqrt(sum(s1)[0] * sum(s2)[0]);
			if (up / (down + 1e-12) > 1.0)
			{
				Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = 1;
			}
			else
			{
				Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = up / (down + 1e-12);
			}

		}
	}
	copyMakeBorder(Coherence, Coherence, win_a, win_a, win_r, win_r, BORDER_REFLECT);
	Coherence.copyTo(coherence);
	return 0;
}

int Utils::complex_coherence(ComplexMat& Mast, ComplexMat& Slave, Mat& coherence)
{
	int wa = 3;  //窗口方位向尺寸
	int wr = 3;  //窗口距离向尺寸

	int na = Mast.GetRows();
	int nr = Mast.GetCols();

	if ((na < 3) ||
		(nr < 3) ||
		Mast.re.type() != CV_64F ||
		Slave.re.type() != CV_64F ||
		Mast.GetCols() != Slave.GetCols() ||
		Mast.GetRows() != Slave.GetRows())
	{
		fprintf(stderr, "complex_coherence(): input check failed!\n\n");
		return -1;
	}

	int win_a = (wa - 1) / 2; //方位窗半径
	int win_r = (wr - 1) / 2; //距离窗半径

	int na_new = na - 2 * win_a;
	int nr_new = nr - 2 * win_r;

	Mat Coherence(na_new, nr_new, CV_64F, Scalar::all(0));
#pragma omp parallel for schedule(guided)
	for (int i = win_a + 1; i <= na - win_a; i++)
	{
		for (int j = win_r + 1; j <= nr - win_r; j++)
		{
			Mat planes_master[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat planes_slave[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat planes[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat s1, s2;
			double up, down, sum1, sum2;
			Mast.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_master[0]);
			Mast.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_master[1]);

			Slave.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_slave[0]);
			Slave.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_slave[1]);

			merge(planes_master, 2, s1);
			merge(planes_slave, 2, s2);
			mulSpectrums(s1, s2, s1, 0, true);
			split(s1, planes);
			sum1 = sum(planes[0])[0];
			sum2 = sum(planes[1])[0];
			up = sqrt(sum1 * sum1 + sum2 * sum2);
			magnitude(planes_master[0], planes_master[1], planes_master[0]);
			magnitude(planes_slave[0], planes_slave[1], planes_slave[0]);
			sum1 = sum(planes_master[0].mul(planes_master[0]))[0];
			sum2 = sum(planes_slave[0].mul(planes_slave[0]))[0];
			down = sqrt(sum1 * sum2);
			Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = up / (down + 0.0000001);
		}
	}
	copyMakeBorder(Coherence, Coherence, 1, 1, 1, 1, BORDER_REFLECT);
	coherence = Coherence;
	return 0;
}

int Utils::complex_coherence(
	const ComplexMat& master_image, 
	const ComplexMat& slave_image,
	int est_wndsize_rg, 
	int est_wndsize_az,
	Mat& coherence
)
{
	int na = master_image.GetRows();
	int nr = master_image.GetCols();

	if ((na < est_wndsize_az) ||
		(nr < est_wndsize_rg) ||
		master_image.type() != CV_64F ||
		slave_image.type() != CV_64F ||
		master_image.GetCols() != slave_image.GetCols() ||
		master_image.GetRows() != slave_image.GetRows() ||
		est_wndsize_az % 2 == 0||
		est_wndsize_rg % 2 == 0||
		est_wndsize_rg < 3||
		est_wndsize_az < 3
		)
	{
		fprintf(stderr, "complex_coherence(): input check failed!\n\n");
		return -1;
	}

	int win_a = (est_wndsize_az - 1) / 2; //方位窗半径
	int win_r = (est_wndsize_rg - 1) / 2; //距离窗半径

	int na_new = na - 2 * win_a;
	int nr_new = nr - 2 * win_r;

	Mat Coherence(na_new, nr_new, CV_64F, Scalar::all(0));
#pragma omp parallel for schedule(guided)
	for (int i = win_a + 1; i <= na - win_a; i++)
	{
		for (int j = win_r + 1; j <= nr - win_r; j++)
		{
			Mat planes_master[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat planes_slave[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat planes[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_64F) };
			Mat s1, s2;
			double up, down, sum1, sum2;
			master_image.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_master[0]);
			master_image.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_master[1]);

			slave_image.re(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_slave[0]);
			slave_image.im(Range(i - 1 - win_a, i + win_a), Range(j - 1 - win_r, j + win_r)).copyTo(planes_slave[1]);

			merge(planes_master, 2, s1);
			merge(planes_slave, 2, s2);
			mulSpectrums(s1, s2, s1, 0, true);
			split(s1, planes);
			sum1 = sum(planes[0])[0];
			sum2 = sum(planes[1])[0];
			up = sqrt(sum1 * sum1 + sum2 * sum2);
			magnitude(planes_master[0], planes_master[1], planes_master[0]);
			magnitude(planes_slave[0], planes_slave[1], planes_slave[0]);
			sum1 = sum(planes_master[0].mul(planes_master[0]))[0];
			sum2 = sum(planes_slave[0].mul(planes_slave[0]))[0];
			down = sqrt(sum1 * sum2);
			Coherence.at<double>(i - 1 - win_a, j - 1 - win_r) = up / (down + 0.0000001);
		}
	}
	copyMakeBorder(Coherence, Coherence, win_a, win_a, win_r, win_r, BORDER_REFLECT);
	Coherence.copyTo(coherence);
	return 0;
}

int Utils::phase_coherence(Mat& phase, Mat& coherence)
{
	if (phase.rows < 3 ||
		phase.cols < 3 ||
		phase.type() != CV_64F ||
		phase.channels() != 1)
	{
		fprintf(stderr, "phase_coherence(): input check failed!\n\n");
		return -1;
	}
	ComplexMat master, slave;
	Mat cos, sin;
	int ret;
	ret = this->phase2cos(phase, cos, sin);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	master.SetRe(cos);
	slave.SetRe(cos);
	master.SetIm(sin);
	sin = -sin;
	slave.SetIm(sin);
	ret = this->complex_coherence(master, slave, coherence);
	if (return_check(ret, "complex_coherence(*, *, *)", error_head)) return -1;
	return 0;
}

int Utils::phase_coherence(const Mat& phase, int est_wndsize_rg, int est_wndsize_az, Mat& coherence)
{
	if (phase.rows < 3 ||
		phase.cols < 3 ||
		phase.type() != CV_64F ||
		phase.channels() != 1 ||
		est_wndsize_rg % 2 == 0||
		est_wndsize_az % 2 == 0
		)
	{
		fprintf(stderr, "phase_coherence(): input check failed!\n\n");
		return -1;
	}
	ComplexMat master, slave;
	Mat cos, sin;
	int ret;
	ret = phase2cos(phase, cos, sin);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	master.SetRe(cos);
	slave.SetRe(cos);
	master.SetIm(sin);
	sin = -sin;
	slave.SetIm(sin);
	ret = complex_coherence(master, slave, est_wndsize_rg, est_wndsize_az, coherence);
	if (return_check(ret, "complex_coherence(*, *, *)", error_head)) return -1;
	return 0;
}

int Utils::phase_derivatives_variance(Mat& phase, Mat& phase_derivatives_variance, int wndsize)
{
	if (phase.cols < 2 ||
		phase.rows < 2 ||
		phase.type() != CV_64F ||
		phase.channels() != 1||
		wndsize < 3||
		wndsize % 2 == 0
		)
	{
		fprintf(stderr, "phase_derivatives_variance(): input check failed!\n\n");
		return -1;
	}
	
	int nr = phase.rows;
	int nc = phase.cols;
	if (wndsize > int(nr / 10) || wndsize > int(nc / 10)) wndsize = 3;
	wndsize = (wndsize - 1) / 2;
	int ret;
	phase.copyTo(phase_derivatives_variance);
	Mat derivative_row, derivative_col;
	derivative_row = phase(Range(1, nr), Range(0, nc)) - phase(Range(0, nr - 1), Range(0, nc));
	derivative_col = phase(Range(0, nr), Range(1, nc)) - phase(Range(0, nr), Range(0, nc - 1));
	ret = wrap(derivative_row, derivative_row);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;
	ret = wrap(derivative_col, derivative_col);
	if (return_check(ret, "wrap(*, *)", error_head)) return -1;
	copyMakeBorder(derivative_col, derivative_col, 0, 1, 0, 1, BORDER_DEFAULT);
	copyMakeBorder(derivative_row, derivative_row, 0, 1, 0, 1, BORDER_DEFAULT);
	copyMakeBorder(derivative_col, derivative_col, wndsize, wndsize, wndsize, wndsize, BORDER_DEFAULT);
	copyMakeBorder(derivative_row, derivative_row, wndsize, wndsize, wndsize, wndsize, BORDER_DEFAULT);
	//cv::meanStdDev()
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double std1, std2, mean1, mean2;
		Mat tmp1, tmp2;
		for (int j = 0; j < nc; j++)
		{
			mean1 = cv::mean(derivative_row(Range(i, i + 2*wndsize + 1), Range(j, j + 2 * wndsize + 1)))[0];
			mean2 = cv::mean(derivative_col(Range(i, i + 2 * wndsize + 1), Range(j, j + 2 * wndsize + 1)))[0];
			tmp1 = derivative_row(Range(i, i + 2 * wndsize + 1), Range(j, j + 2 * wndsize + 1)) - mean1;
			tmp2 = derivative_col(Range(i, i + 2 * wndsize + 1), Range(j, j + 2 * wndsize + 1)) - mean2;
			tmp1 = tmp1.mul(tmp1);
			tmp2 = tmp2.mul(tmp2);
			std1 = sum(tmp1)[0];
			std2 = sum(tmp2)[0];
			phase_derivatives_variance.at<double>(i, j) = (sqrt(std1) + sqrt(std2)) / ((2 * wndsize + 1) * (2 * wndsize + 1));
		}
	}
	return 0;
}

int Utils::max_integrable_distance(Mat& phase, Mat& max_integrable_distance, double conservative_thresh)
{
	if (phase.rows < 5 ||
		phase.cols < 5 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F ||
		conservative_thresh < 1.41
		)
	{
		fprintf(stderr, "max_integrable_distance(): input check failed!\n\n");
		return -1;
	}
	int nr = phase.rows;
	int nc = phase.cols;
	max_integrable_distance = Mat::ones(nr, nc, CV_64F);
	max_integrable_distance = max_integrable_distance * 1.5;
#pragma omp parallel for schedule(guided)
	for (int i = 1; i < nr - 1; i++)
	{
		for (int j = 1; j < nc - 1; j++)
		{
			double max = -1.0;
			double delta, ph0;
			ph0 = phase.at<double>(i, j);

			delta = phase.at<double>(i - 1, j - 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.414;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i - 1, j) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.0;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i - 1, j + 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.414;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i, j - 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.0;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i, j + 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.0;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i + 1, j - 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.414;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i + 1, j) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.0;
			max = max > delta ? max : delta;

			delta = phase.at<double>(i + 1, j + 1) - ph0;
			delta = fabs(atan2(sin(delta), cos(delta))) / 1.414;
			max = max > delta ? max : delta;

			if (max < 1.0 / conservative_thresh * PI)
			{
				max = 1.0 / conservative_thresh * PI;
			}
			max_integrable_distance.at<double>(i, j) = PI / max < 1.5 ? 1.5 : PI / max;
		}
	}

	return 0;
}

int Utils::fftshift(Mat& mag)
{
	// rearrange the quadrants of Fourier image
	// so that the origin is at the image center
	if (mag.rows < 2 ||
		mag.cols < 2 ||
		mag.channels() != 1)
	{
		fprintf(stderr, "fftshift(): input check failed!\n\n");
		return -1;
	}
	mag = mag(Rect(0, 0, mag.cols & -2, mag.rows & -2));
	int cx = mag.cols / 2;
	int cy = mag.rows / 2;
	Mat tmp;
	Mat q0(mag, Rect(0, 0, cx, cy));
	Mat q1(mag, Rect(cx, 0, cx, cy));
	Mat q2(mag, Rect(0, cy, cx, cy));
	Mat q3(mag, Rect(cx, cy, cx, cy));

	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);
	return 0;
}

int Utils::read_DIMACS(const char* DIMACS_file_solution, tri_edge* edges, int num_edges, vector<tri_node>& nodes, triangle* tri, int num_triangle)
{
	if (DIMACS_file_solution == NULL ||
		edges == NULL ||
		num_edges < 3 ||
		nodes.size() < 3 ||
		tri == NULL ||
		num_triangle < 1
		)
	{
		fprintf(stderr, "read_DIMACS(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(DIMACS_file_solution, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_DIMACS(): can't open %s \n", DIMACS_file_solution);
		return -1;
	}

	char instring[256];
	char ch;
	double obj_value = 0;
	int i, tmp, end1, end2, end3, row1, col1, row2, col2, row3, col3;
	int end[3];
	double x1, y1, x2, y2, direction;
	long from, to;
	double flow = 0;
	bool flag;
	int x[3];
	int y[3];
	long* ptr_neigh = NULL;
	int num_neigh, target_edges;
	int num_nodes = nodes.size();
	/////////////////////读取注释///////////////////////////
	GET_NEXT_LINE;
	while (ch != 's' && ch)
	{
		if (ch != 'c')
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}
		GET_NEXT_LINE;
	}
	/////////////////////读取优化目标值/////////////////////
	for (i = 1; i < 81; i++)
	{
		if (isspace((int)instring[i]) > 0)
		{
			i++;
			break;
		}
	}
	if (sscanf(&(instring[i]), "%lf", &obj_value) != 1)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
		return -1;
	}
	if (obj_value < 0.0)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): this problem can't be solved(unbounded or infeasible)!\n\n");
		return -1;
	}
	////////////////////////读取MCF结果////////////////////////
	GET_NEXT_LINE;
	while (ch && ch == 'f')
	{
		if (sscanf(&(instring[2]), "%ld %ld %lf", &from, &to, &flow) != 3 ||
			flow < 0.0 || from < 0 || to < 0)
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}

		if (from > 0 &&
			from <= num_triangle &&
			to > 0 &&
			to <= num_triangle)
		{
			if (/*from > 0 &&
				from <= num_triangle &&
				to > 0 &&
				to <= num_triangle &&*/
				(tri + from - 1) != NULL &&
				(tri + to - 1) != NULL
				)
			{
				end[0] = -1;
				end[1] = -1;
				x[0] = (tri + from - 1)->p1;
				x[1] = (tri + from - 1)->p2;
				x[2] = (tri + from - 1)->p3;
				y[0] = (tri + to - 1)->p1;
				y[1] = (tri + to - 1)->p2;
				y[2] = (tri + to - 1)->p3;
				i = 0; tmp = 0; flag = false;
				while (end[0] == -1 || end[1] == -1)
				{
					if (i > 2)
					{
						if (fp) fclose(fp);
						fprintf(stderr, "read_DIMACS(): illegal Delaunay triangle!\n\n");
						return -1;
					}
					for (int j = 0; j < 3; j++)
					{
						if (x[i] == y[j])
						{
							end[tmp] = x[i];
							tmp++;
							if (i == 0) flag = true;
							break;
						}
					}
					i++;

				}
			}
			if (i == 3)
			{
				if (flag) end[2] = x[1];
				else
				{
					end[2] = x[0];
				}
			}
			else
			{
				end[2] = x[2];
			}
			if (end[0] > end[1])
			{
				end1 = end[1];
				end2 = end[0];
			}
			else
			{
				end1 = end[0];
				end2 = end[1];
			}
			end3 = end[2];




			if (end1 > 0 && end1 <= num_nodes && end2 > 0 && end2 <= num_nodes && end3 > 0 && end3 <= num_nodes)
			{
				//找到边序号target_edges
				nodes[end1 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
				for (i = 0; i < num_neigh; i++)
				{
					if ((ptr_neigh + i) != NULL && *(ptr_neigh + i) > 0 && *(ptr_neigh + i) <= num_edges)
					{
						if ((edges + *(ptr_neigh + i) - 1)->end1 == end2 || (edges + *(ptr_neigh + i) - 1)->end2 == end2)
						{
							target_edges = *(ptr_neigh + i);
						}
					}
				}

				nodes[end1 - 1].get_pos(&row1, &col1);
				nodes[end2 - 1].get_pos(&row2, &col2);
				nodes[end3 - 1].get_pos(&row3, &col3);
				x1 = double(col1 - col3);
				y1 = double(row3 - row1);
				x2 = double(col2 - col1);
				y2 = double(row1 - row2);
				direction = x1 * y2 - x2 * y1;
				if (direction > 0.0)//在目标三角形中顺残差方向
				{
					(edges + target_edges - 1)->gain = flow;
				}
				else
				{
					(edges + target_edges - 1)->gain = -flow;
				}
			}
		}
		GET_NEXT_LINE;
	}

	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int Utils::read_DIMACS(
	const char* DIMACS_file_solution,
	vector<tri_edge>& edges, 
	vector<tri_node>& nodes,
	vector<triangle>& triangle
)
{
	if (DIMACS_file_solution == NULL ||
		edges.size() < 3 ||
		nodes.size() < 3 ||
		triangle.size() < 1
		)
	{
		fprintf(stderr, "read_DIMACS(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(DIMACS_file_solution, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_DIMACS(): can't open %s \n", DIMACS_file_solution);
		return -1;
	}

	char instring[256];
	char ch;
	double obj_value = 0;
	int i, tmp, end1, end2, end3, row1, col1, row2, col2, row3, col3;
	int end[3];
	double x1, y1, x2, y2, direction;
	long from, to;
	double flow = 0;
	bool flag;
	int x[3];
	int y[3];
	long* ptr_neigh = NULL;
	int num_neigh, target_edges;
	int num_nodes = nodes.size();
	int num_triangle = triangle.size(); int num_edges = edges.size();
	/////////////////////读取注释///////////////////////////
	GET_NEXT_LINE;
	while (ch != 's' && ch)
	{
		if (ch != 'c')
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}
		GET_NEXT_LINE;
	}
	/////////////////////读取优化目标值/////////////////////
	for (i = 1; i < 81; i++)
	{
		if (isspace((int)instring[i]) > 0)
		{
			i++;
			break;
		}
	}
	if (sscanf(&(instring[i]), "%lf", &obj_value) != 1)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
		return -1;
	}
	if (obj_value < 0.0)
	{
		if (fp) fclose(fp);
		fprintf(stderr, "read_DIMACS(): this problem can't be solved(unbounded or infeasible)!\n\n");
		return -1;
	}
	////////////////////////读取MCF结果////////////////////////
	GET_NEXT_LINE;
	while (ch && ch == 'f')
	{
		if (sscanf(&(instring[2]), "%ld %ld %lf", &from, &to, &flow) != 3 ||
			flow < 0.0 || from < 0 || to < 0)
		{
			if (fp) fclose(fp);
			fprintf(stderr, "read_DIMACS(): unknown file format!\n\n");
			return -1;
		}
		//非接地边
		if (from > 0 &&
			from <= num_triangle &&
			to > 0 &&
			to <= num_triangle)
		{
			/////////寻找两个三角形的公共边//////////////
			{
				end[0] = -1;
				end[1] = -1;
				x[0] = triangle[from - 1].p1;
				x[1] = triangle[from - 1].p2;
				x[2] = triangle[from - 1].p3;
				y[0] = triangle[to - 1].p1;
				y[1] = triangle[to - 1].p2;
				y[2] = triangle[to - 1].p3;
				i = 0; tmp = 0; flag = false;
				while (end[0] == -1 || end[1] == -1)
				{
					if (i > 2)
					{
						if (fp) fclose(fp);
						fprintf(stderr, "read_DIMACS(): illegal Delaunay triangle!\n\n");
						return -1;
					}
					for (int j = 0; j < 3; j++)
					{
						if (x[i] == y[j])
						{
							end[tmp] = x[i];
							tmp++;
							if (i == 0) flag = true;
							break;
						}
					}
					i++;

				}
			}
			if (i == 3)
			{
				if (flag) end[2] = x[1];
				else
				{
					end[2] = x[0];
				}
			}
			else
			{
				end[2] = x[2];
			}
			if (end[0] > end[1])
			{
				end1 = end[1];
				end2 = end[0];
			}
			else
			{
				end1 = end[0];
				end2 = end[1];
			}
			end3 = end[2];
			/////////寻找两个三角形的公共边//////////////



			if (end1 > 0 && end1 <= num_nodes && end2 > 0 && end2 <= num_nodes && end3 > 0 && end3 <= num_nodes)
			{
				//找到边序号target_edges
				nodes[end1 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
				for (i = 0; i < num_neigh; i++)
				{
					if ((ptr_neigh + i) != NULL && *(ptr_neigh + i) > 0 && *(ptr_neigh + i) <= num_edges)
					{
						if (edges[*(ptr_neigh + i) - 1].end1 == end2 || edges[*(ptr_neigh + i) - 1].end2 == end2)
						{
							target_edges = *(ptr_neigh + i);
						}
					}
				}

				nodes[end1 - 1].get_pos(&row1, &col1);
				nodes[end2 - 1].get_pos(&row2, &col2);
				nodes[end3 - 1].get_pos(&row3, &col3);
				x1 = double(col1 - col3);
				y1 = -double(row3 - row1);
				x2 = double(col2 - col1);
				y2 = -double(row1 - row2);
				direction = x1 * y2 - x2 * y1;
				if (direction < 0.0)//在目标三角形中顺残差方向
				{
					edges[target_edges - 1].gain = -flow;
				}
				else
				{
					edges[target_edges - 1].gain = flow;
				}
			}
		}
		//接地边
		if (from == num_triangle + 1 || to == num_triangle + 1)
		{
			if (from == num_triangle + 1)
			{
				if (edges[triangle[to - 1].edge1 - 1].isBoundry)target_edges = triangle[to - 1].edge1;
				else if (edges[triangle[to - 1].edge2 - 1].isBoundry) target_edges = triangle[to - 1].edge2;
				else target_edges = triangle[to - 1].edge3;


				if (edges[target_edges - 1].end1 > edges[target_edges - 1].end2)
				{
					end1 = edges[target_edges - 1].end2;
					end2 = edges[target_edges - 1].end1;
				}
				else
				{
					end1 = edges[target_edges - 1].end1;
					end2 = edges[target_edges - 1].end2;
				}

				if (triangle[to - 1].p1 != end1 && triangle[to - 1].p1 != end2) end3 = triangle[to - 1].p1;
				else if (triangle[to - 1].p2 != end1 && triangle[to - 1].p2 != end2) end3 = triangle[to - 1].p2;
				else end3 = triangle[to - 1].p3;

				nodes[end1 - 1].get_pos(&row1, &col1);
				nodes[end2 - 1].get_pos(&row2, &col2);
				nodes[end3 - 1].get_pos(&row3, &col3);
				x1 = double(col1 - col3);
				y1 = -double(row3 - row1);
				x2 = double(col2 - col1);
				y2 = -double(row1 - row2);
				direction = x1 * y2 - x2 * y1;
				if (direction < 0.0)//在目标三角形中顺残差方向
				{
					edges[target_edges - 1].gain = flow;
				}
				else
				{
					edges[target_edges - 1].gain = -flow;
				}
			}
			else
			{
				if (edges[triangle[from - 1].edge1 - 1].isBoundry)target_edges = triangle[from - 1].edge1;
				else if (edges[triangle[from - 1].edge2 - 1].isBoundry) target_edges = triangle[from - 1].edge2;
				else target_edges = triangle[from - 1].edge3;


				if (edges[target_edges - 1].end1 > edges[target_edges - 1].end2)
				{
					end1 = edges[target_edges - 1].end2;
					end2 = edges[target_edges - 1].end1;
				}
				else
				{
					end1 = edges[target_edges - 1].end1;
					end2 = edges[target_edges - 1].end2;
				}

				if (triangle[from - 1].p1 != end1 && triangle[from - 1].p1 != end2) end3 = triangle[from - 1].p1;
				else if (triangle[from - 1].p2 != end1 && triangle[from - 1].p2 != end2) end3 = triangle[from - 1].p2;
				else end3 = triangle[from - 1].p3;

				nodes[end1 - 1].get_pos(&row1, &col1);
				nodes[end2 - 1].get_pos(&row2, &col2);
				nodes[end3 - 1].get_pos(&row3, &col3);
				x1 = double(col1 - col3);
				y1 = -double(row3 - row1);
				x2 = double(col2 - col1);
				y2 = -double(row1 - row2);
				direction = x1 * y2 - x2 * y1;
				if (direction < 0.0)//在目标三角形中顺残差方向
				{
					edges[target_edges - 1].gain = -flow;
				}
				else
				{
					edges[target_edges - 1].gain = flow;
				}
			}
		}
		GET_NEXT_LINE;
	}

	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int Utils::cvmat2bin(const char* filename, Mat& mat)
{
	int nr = mat.rows;
	int nc = mat.cols;
	if (nr < 1 || nc < 1 || mat.type() != CV_64F || mat.channels() != 1)
	{
		fprintf(stderr, "cvmat2bin(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fopen_s(&fp, filename, "wb");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open file: %s\n", filename);
		return -1;
	}
	fwrite(&nr, sizeof(int), 1, fp);
	fwrite(&nc, sizeof(int), 1, fp);
	fwrite((double*)mat.data, sizeof(double), nr * nc, fp);
	if (fp != NULL) fclose(fp);
	return 0;
}

int Utils::bin2cvmat(const char* filename, Mat& dst)
{

	FILE* fp = NULL;
	fopen_s(&fp, filename, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open file: %s\n", filename);
		return -1;
	}
	int rows, cols;
	fread(&rows, sizeof(int), 1, fp);
	fread(&cols, sizeof(int), 1, fp);
	if (rows < 1 || cols < 1)
	{
		fprintf(stderr, "unknown file format!\n");
		if (fp) fclose(fp);
		return -1;
	}
	Mat matrix(rows, cols, CV_64F, cv::Scalar::all(0));
	double* p = (double*)malloc(sizeof(double) * rows * cols);
	if (p == NULL)
	{
		fprintf(stderr, "failed to allocate memory for reading data from %s!\n", filename);
		if (fp) fclose(fp);
		return -1;
	}
	fread(p, sizeof(double), rows * cols, fp);
	std::memcpy(matrix.data, p, sizeof(double) * rows * cols);
	if (p != NULL) free(p);
	if (fp != NULL) fclose(fp);
	dst = matrix;
	return 0;
}

int Utils::multilook(ComplexMat& Master, ComplexMat& Slave, Mat& phase, int multilook_times)
{
	if (Master.GetRows() != Slave.GetRows() ||
		Master.GetCols() != Slave.GetCols() ||
		Master.type() != CV_64F ||
		Slave.type() != CV_64F ||
		Master.GetRows() < 1 ||
		Master.GetCols() < 1 ||
		multilook_times < 1 ||
		Master.GetRows() < multilook_times||
		Master.GetCols() < multilook_times)
	{
		fprintf(stderr, "multilook(): input check failed!\n\n");
		return -1;
	}
	int ret;
	if (multilook_times == 1)
	{
		ret = generate_phase(Master, Slave, phase);
		if (return_check(ret, "generate_phase(*, *, *)", error_head)) return -1;
		return 0;
	}
	ComplexMat tmp;
	ret = Master.Mul(Slave, tmp, true);
	if (return_check(ret, "Master.Mul(*, *, *)", error_head)) return -1;
	int nr = tmp.GetRows();
	int nc = tmp.GetCols();
	nr = (nr - (nr % multilook_times)) / multilook_times;
	nc = (nc - (nc % multilook_times)) / multilook_times;
	Mat real = Mat::zeros(nr, nc, CV_64F);
	Mat imag = Mat::zeros(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			real.at<double>(i, j) = cv::mean(tmp.re(Range(i * multilook_times, (i + 1) * multilook_times),
				Range(j * multilook_times, (j + 1) * multilook_times)))[0];
			imag.at<double>(i, j) = cv::mean(tmp.im(Range(i * multilook_times, (i + 1) * multilook_times),
				Range(j * multilook_times, (j + 1) * multilook_times)))[0];
		}
	}
	tmp.SetRe(real);
	tmp.SetIm(imag);
	phase = tmp.GetPhase();
	return 0;
}

int Utils::multilook(const ComplexMat& master, const ComplexMat& slave, int multilook_rg, int multilook_az, Mat& phase)
{
	if (master.GetRows() != slave.GetRows() ||
		master.GetCols() != slave.GetCols() ||
		master.type() != CV_64F ||
		slave.type() != CV_64F ||
		master.GetRows() < 1 ||
		master.GetCols() < 1 ||
		multilook_rg < 1 ||
		multilook_az < 1 ||
		master.GetRows() < multilook_az ||
		master.GetCols() < multilook_rg)
	{
		fprintf(stderr, "multilook(): input check failed!\n\n");
		return -1;
	}
	int ret;
	if (multilook_rg == 1 && multilook_az == 1)
	{
		ret = generate_phase(master, slave, phase);
		if (return_check(ret, "generate_phase(*, *, *)", error_head)) return -1;
		return 0;
	}
	ComplexMat tmp;
	ret = master.Mul(slave, tmp, true);
	if (return_check(ret, "Master.Mul(*, *, *)", error_head)) return -1;
	int nr = tmp.GetRows();
	int nc = tmp.GetCols();
	int radius_rg = multilook_rg / 2;
	int radius_az = multilook_az / 2;
	Mat real = Mat::zeros(nr, nc, CV_64F);
	Mat imag = Mat::zeros(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		int left, right, bottom, top;
		for (int j = 0; j < nc; j++)
		{
			left = j - radius_rg; left = left < 0 ? 0 : left;
			right = left + multilook_rg; right = right > nc - 1 ? nc - 1 : right;
			top = i - radius_az; top = top < 0 ? 0 : top;
			bottom = top + multilook_az; bottom = bottom > nr - 1 ? nr - 1 : bottom;
			real.at<double>(i, j) = cv::mean(tmp.re(Range(top, bottom + 1),Range(left, right + 1)))[0];
			imag.at<double>(i, j) = cv::mean(tmp.im(Range(top, bottom + 1), Range(left, right + 1)))[0];
		}
	}
	tmp.SetRe(real);
	tmp.SetIm(imag);
	tmp.GetPhase().copyTo(phase);
	return 0;
}

int Utils::phase2cos(const Mat& phase, Mat& cos, Mat& sin)
{
	if (phase.rows < 1 ||
		phase.cols < 1 ||
		phase.type() != CV_64F ||
		phase.channels() != 1)
	{
		fprintf(stderr, "phase2cos(): input check failed!\n\n");
		return -1;
	}
	int nr = phase.rows;
	int nc = phase.cols;
	Mat Cos(nr, nc, CV_64F);
	Mat Sin(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			Cos.at<double>(i, j) = std::cos(phase.at<double>(i, j));
			Sin.at<double>(i, j) = std::sin(phase.at<double>(i, j));
		}
	}
	Cos.copyTo(cos);
	Sin.copyTo(sin);
	return 0;
}

int Utils::xyz2ell(Mat xyz, Mat& llh)
{
	if (xyz.rows != 1 ||
		xyz.cols != 3 ||
		xyz.type() != CV_64F ||
		xyz.channels() != 1)
	{
		fprintf(stderr, "xyz2ell(): input check failed!\n\n");
		return -1;
	}
	double x = xyz.at<double>(0, 0);
	double y = xyz.at<double>(0, 1);
	double z = xyz.at<double>(0, 2);
	double Rad_earth_e = 6378136.49;
	double f = 1 / 298.257223563;
	double t = Rad_earth_e * (1 - f);
	double e = sqrt((Rad_earth_e * Rad_earth_e - t * t) / (Rad_earth_e * Rad_earth_e));
	double r = x * x + y * y;
	if (fabs(r) < 1e-15)
	{
		r = 1e-14;
	}
	if (fabs(x) < 1e-15)
	{
		x = 1e-14;
	}
	double lat = atan(z / r);
	double lon = atan(y / x);
	double N, height;
	double tmp = 0.0;
	double pi = 3.1415926535;
	for (int k = 0; k < 10; k++)
	{
		tmp = (sqrt(1 - e * e * sin(lat) * sin(lat)));
		if (fabs(tmp) < 1e-15)
		{
			tmp = 1e-14;
		}
		N = Rad_earth_e / tmp;
		tmp = sin(lat);
		if (fabs(tmp) < 1e-15)
		{
			tmp = 1e-14;
		}
		height = z / tmp - N * (1 - e * e);
		tmp = (sqrt(x * x + y * y) * (N * (1 - e * e) + height));
		if (fabs(tmp) < 1e-15)
		{
			tmp = 1e-14;
		}
		lat = atan(z * (N + height) / tmp);
	}
	if (x > 0 && y < 0)
	{
		lon = -lon;
	}
	else if (x < 0 && y > 0)
	{
		lon = pi + lon;
	}
	else if (x < 0 && y < 0)
	{
		lon = -pi + lon;
	}
	lat = lat * 180 / pi;
	lon = lon * 180 / pi;
	llh = Mat::zeros(1, 3, CV_64F);
	llh.at<double>(0, 0) = lat;
	llh.at<double>(0, 1) = lon;
	llh.at<double>(0, 2) = height;
	return 0;
}

int Utils::ell2xyz(Mat llh, Mat& xyz)
{
	if (llh.cols != 3 ||
		llh.rows != 1 ||
		llh.type() != CV_64F ||
		llh.channels() != 1
		)
	{
		fprintf(stderr, "ell2xyz(): input check failed!\n\n");
		return -1;
	}
	double e2 = 0.00669438003551279091;
	double lat = llh.at<double>(0, 0);
	double lon = llh.at<double>(0, 1);
	double height = llh.at<double>(0, 2);
	lat = lat / 180.0 * PI;
	lon = lon / 180.0 * PI;
	double Ea = 6378136.49;
	double N = Ea / sqrt(1 - e2 * (sin(lat) * sin(lat)));
	double Nph = N + height;
	double x = Nph * cos(lat) * cos(lon);
	double y = Nph * cos(lat) * sin(lon);
	double z = (Nph - e2 * N) * sin(lat);
	Mat tmp = Mat::zeros(1, 3, CV_64F);
	tmp.at<double>(0, 0) = x;
	tmp.at<double>(0, 1) = y;
	tmp.at<double>(0, 2) = z;
	tmp.copyTo(xyz);
	return 0;
}





int Utils::saveSLC(const char* filename, double db, ComplexMat& SLC)
{
	if (filename == NULL ||
		db < 0 ||
		SLC.GetRows() < 1 ||
		SLC.GetCols() < 1 /*||
		SLC.type() != CV_64F*/)
	{
		fprintf(stderr, "saveSLC(): input check failed!\n\n");
		return -1;
	}
	ComplexMat tmp;
	Mat mod;
	if (SLC.type() != CV_64F && SLC.type() != CV_32F)
	{
		SLC.re.convertTo(tmp.re, CV_32F);
		SLC.im.convertTo(tmp.im, CV_32F);
		mod = tmp.GetMod();
	}
	else
	{
		mod = SLC.GetMod();
	}
	int nr = mod.rows;
	int nc = mod.cols;
	double max, min;
	if (SLC.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				mod.at<double>(i, j) = 20 * log10(mod.at<double>(i, j) + 0.000001);
			}
		}

		minMaxLoc(mod, &min, &max);
		if (fabs(max - min) < 0.00000001)
		{
			fprintf(stderr, "SLC image intensity is the same for every pixel\n\n");
			return -1;
		}
		min = max - db;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if (mod.at<double>(i, j) <= min) mod.at<double>(i, j) = min;
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				mod.at<float>(i, j) = 20 * log10(mod.at<float>(i, j) + 0.000001);
			}
		}

		minMaxLoc(mod, &min, &max);
		if (fabs(max - min) < 0.00000001)
		{
			fprintf(stderr, "SLC image intensity is the same for every pixel\n\n");
			return -1;
		}
		min = max - db;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if (mod.at<float>(i, j) <= min) mod.at<float>(i, j) = min;
			}
		}
	}
	mod = (mod - min) / (max - min) * 255.0;
	mod.convertTo(mod, CV_8U);
	bool ret = cv::imwrite(filename, mod);
	if (!ret)
	{
		fprintf(stderr, "cv::imwrite(): can't write to %s\n\n", filename);
		return -1;
	}
	return 0;
}

int Utils::savephase(const char* filename, const char* colormap, Mat phase)
{
	if (filename == NULL ||
		colormap == NULL ||
		phase.rows < 1 ||
		phase.cols < 1 ||
		phase.type() != CV_64F ||
		phase.channels() != 1)
	{
		fprintf(stderr, "savephase(): input check failed!\n\n");
		return -1;
	}
	bool gray = false;
	cv::ColormapTypes type = cv::COLORMAP_PARULA;
	if (strcmp(colormap, "jet") == 0) type = cv::COLORMAP_JET;
	if (strcmp(colormap, "hsv") == 0) type = cv::COLORMAP_HSV;
	if (strcmp(colormap, "cool") == 0) type = cv::COLORMAP_COOL;
	if (strcmp(colormap, "rainbow") == 0) type = cv::COLORMAP_RAINBOW;
	if (strcmp(colormap, "spring") == 0) type = cv::COLORMAP_SPRING;
	if (strcmp(colormap, "summer") == 0) type = cv::COLORMAP_SUMMER;
	if (strcmp(colormap, "winter") == 0) type = cv::COLORMAP_WINTER;
	if (strcmp(colormap, "autumn") == 0) type = cv::COLORMAP_AUTUMN;
	if (strcmp(colormap, "gray") == 0) gray = true;

	double min, max;
	Mat tmp;
	phase.copyTo(tmp);
	cv::minMaxLoc(tmp, &min, &max);
	if (fabs(max - min) < 0.000001)
	{
		fprintf(stderr, "phase value is the same for every pixel\n\n");
		return -1;
	}
	tmp = (tmp - min) / (max - min)*255.0;
	tmp.convertTo(tmp, CV_8U);
	if (!gray)
	{
		cv::applyColorMap(tmp, tmp, type);
	}
	bool ret = cv::imwrite(filename, tmp);
	if (!ret)
	{
		fprintf(stderr, "cv::imwrite(): can't write to %s\n\n", filename);
		return -1;
	}
	return 0;
}

int Utils::resampling(const char* Src_file, const char* Dst_file, int dst_rows, int dst_cols)
{
	if (Src_file == NULL ||
		Dst_file == NULL ||
		dst_rows < 1 ||
		dst_cols < 1)
	{
		fprintf(stderr, "down_sampling(): input check failed!\n\n");
		return -1;
	}
	Mat img = cv::imread(Src_file);
	if (img.rows < 1 || img.cols < 1)
	{
		fprintf(stderr, "can't read from %s!\n\n", Dst_file);
		return -1;
	}
	cv::resize(img, img, cv::Size(dst_cols, dst_rows));
	if (cv::imwrite(Dst_file, img) == false)
	{
		fprintf(stderr, "failed to write %s\n\n", Dst_file);
		return -1;
	}
	return 0;
}

int Utils::amplitude_phase_blend(const char* amplitude_file, const char* phase_file, const char* blended_file, double SAR_ratio)
{
	if (amplitude_file == NULL ||
		phase_file == NULL ||
		blended_file == NULL
		)
	{
		fprintf(stderr, "amplitude_phase_blend(): input check failed!\n\n");
		return -1;
	}
	Mat amplitude = imread(amplitude_file);
	Mat phase = imread(phase_file);
	if (amplitude.rows < 1 ||
		amplitude.cols < 1
		)
	{
		fprintf(stderr, "amplitude_phase_blend(): can't open %s!\n\n", amplitude_file);
		return -1;
	}
	if (phase.rows < 1 ||
		phase.cols < 1
		)
	{
		fprintf(stderr, "amplitude_phase_blend(): can't open %s!\n", phase_file);
		return -1;
	}
	SAR_ratio = SAR_ratio < 0.5 ? 0.5 : SAR_ratio;
	SAR_ratio = SAR_ratio > 0.99 ? 0.99 : SAR_ratio;
	Mat blend;
	addWeighted(amplitude, SAR_ratio, phase, 1.0 - SAR_ratio, 0, blend);
	if (blend.rows < 1 || blend.cols < 1)
	{
		fprintf(stderr, "amplitude_phase_blend(): failed to blend !\n");
		return -1;
	}
	if (!imwrite(blended_file, blend))
	{
		fprintf(stderr, "amplitude_phase_blend(): failed to write to %s !\n", blended_file);
		return -1;
	}
	return 0;
}



int Utils::read_edges(const char* filename, tri_edge** edges, long* num_edges, int** neighbours, long num_nodes)
{
	if (filename == NULL ||
		num_edges == NULL ||
		num_nodes < 3||
		edges == NULL||
		neighbours == NULL)
	{
		fprintf(stderr, "read_edges(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(filename, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_edges(): can't open %s\n", filename);
		return -1;
	}
	char str[1024];
	char* ptr;
	fgets(str, 1024, fp);
	*num_edges = strtol(str, &ptr, 0);
	if (*num_edges <= 0)
	{
		fprintf(stderr, "read_edges(): %s is unknown format!\n", filename);
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	*edges = (tri_edge*)malloc(*num_edges * sizeof(tri_edge));
	if (*edges == NULL)
	{
		fprintf(stderr, "read_edges(): unreasonable number of edges, out of memory!\n");
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	memset(*edges, 0, *num_edges * sizeof(tri_edge));
	*neighbours = (int*)malloc(sizeof(int) * num_nodes);
	if (*neighbours == NULL)
	{
		fprintf(stderr, "read_edges(): unreasonable number of nodes, out of memory!\n");
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		if (*edges)
		{
			free(*edges);
			*edges = NULL;
		}
		return -1;
	}
	memset(*neighbours, 0, sizeof(int) * num_nodes);
	long end1, end2, edges_number, boundry_marker;
	for (int i = 0; i < *num_edges; i++)
	{
		fgets(str, 1024, fp);
		edges_number = strtol(str, &ptr, 0);
		end1 = strtol(ptr, &ptr, 0);
		end2 = strtol(ptr, &ptr, 0);
		boundry_marker = strtol(ptr, &ptr, 0);
		(*edges + i)->end1 = end1;
		(*edges + i)->end2 = end2;
		(*edges + i)->num = i + 1;
		(*edges + i)->gain = 0;
		(*edges + i)->isResidueEdge = false;
		(*edges + i)->isBoundry = (boundry_marker == 1);
		if (end1 < 1 ||
			end1 > num_nodes ||
			end2 < 1 ||
			end2 > num_nodes)
		{
			fprintf(stderr, "read_edges(): endpoints exceed 1~num_nodes!\n");
			if (fp)
			{
				fclose(fp);
				fp = NULL;
			}
			if (*edges)
			{
				free(*edges);
				*edges = NULL;
			}
			if (*neighbours)
			{
				free(*neighbours);
				*neighbours = NULL;
			}
			return -1;
		}
		*(*neighbours + end1 - 1) = *(*neighbours + end1 - 1) + 1;//统计每个节点有多少邻接边
		*(*neighbours + end2 - 1) = *(*neighbours + end2 - 1) + 1;
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int Utils::read_edges(const char* edge_file, vector<tri_edge>& edges, std::vector<int>& node_neighbours, long num_nodes)
{
	if (edge_file == NULL ||
		num_nodes < 3)
	{
		fprintf(stderr, "read_edges(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(edge_file, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_edges(): can't open %s\n", edge_file);
		return -1;
	}
	char str[1024];
	char* ptr;
	fgets(str, 1024, fp);
	long long num_edges = 0;
	num_edges = strtol(str, &ptr, 0);
	if (num_edges <= 0)
	{
		fprintf(stderr, "read_edges(): %s is unknown format!\n", edge_file);
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	edges.clear(); node_neighbours.clear();
	edges.resize(num_edges);
	node_neighbours.resize(num_nodes);
	long end1, end2, edges_number, boundry_marker;
	for (int i = 0; i < num_edges; i++)
	{
		fgets(str, 1024, fp);
		edges_number = strtol(str, &ptr, 0);
		end1 = strtol(ptr, &ptr, 0);
		end2 = strtol(ptr, &ptr, 0);
		boundry_marker = strtol(ptr, &ptr, 0);
		edges[i].end1 = end1;
		edges[i].end2 = end2;
		edges[i].num = i + 1;
		edges[i].gain = 0;
		edges[i].isResidueEdge = false;
		edges[i].phase_diff = 0.0;
		edges[i].isBoundry = (boundry_marker == 1);
		if (end1 < 1 ||
			end1 > num_nodes ||
			end2 < 1 ||
			end2 > num_nodes)
		{
			fprintf(stderr, "read_edges(): endpoints exceed 1~num_nodes!\n");
			if (fp)
			{
				fclose(fp);
				fp = NULL;
			}
			return -1;
		}
		node_neighbours[end1 - 1] += 1;//统计每个节点有多少邻接边
		node_neighbours[end2 - 1] += 1;//统计每个节点有多少邻接边
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int Utils::init_tri_node(vector<tri_node>& node_array, Mat& phase, Mat& mask, tri_edge* edges, long num_edges, int* num_neighbour, int num_nodes)
{
	if (phase.rows < 2 ||
		phase.cols < 2 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F ||
		mask.rows != phase.rows ||
		mask.cols != phase.cols ||
		mask.channels() != 1 ||
		mask.type() != CV_32S ||
		edges == NULL ||
		num_edges < 3 ||
		num_neighbour == NULL ||
		num_nodes < 3)
	{
		fprintf(stderr, "init_tri_node(): input check failed!\n\n");
		return -1;
	}
	int sum = cv::countNonZero(mask);;
	if (sum != num_nodes)
	{
		fprintf(stderr, "init_tri_node(): mask and num_nodes mismatch!\n\n");
		return -1;
	}
	int rows = phase.rows;
	int cols = phase.cols;
	int count = 0;
	tri_node* ptr = NULL;
	double Phase;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				Phase = phase.at<double>(i, j);
				ptr = new tri_node(i, j, *(num_neighbour + count), Phase);
				if (mask.at<int>(i, j) > 1)//邻接已解缠节点标记
				{
					ptr->set_status(true);
				}
				node_array.push_back(*ptr);
				delete ptr;
				ptr = NULL;
				count++;
				
			}
		}
	}

	long* neighbour_ptr = NULL;
	int dummy, ret;
	tri_edge tmp;
	for (int i = 0; i < num_edges; i++)
	{
		tmp = *(edges + i);
		if (tmp.end1 < 0 ||
			tmp.end2 < 0 ||
			tmp.end1 > node_array.size() ||
			tmp.end2 > node_array.size())
		{
			fprintf(stderr, "init_tri_node(): edges' endpoint exceed legal value!\n\n");
			return -1;
		}
		ret = node_array[tmp.end1 - 1].get_neigh_ptr(&neighbour_ptr, &dummy);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;

		ret = node_array[tmp.end2 - 1].get_neigh_ptr(&neighbour_ptr, &dummy);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;
	}
	return 0;
}

int Utils::init_tri_node(
	vector<tri_node>& node_array,
	const Mat& phase,
	const Mat& mask,
	const vector<tri_edge>& edges,
	const vector<int>& node_neighbours,
	int num_nodes
)
{
	if (phase.rows < 2 ||
		phase.cols < 2 ||
		phase.channels() != 1 ||
		phase.type() != CV_64F ||
		mask.rows != phase.rows ||
		mask.cols != phase.cols ||
		mask.channels() != 1 ||
		mask.type() != CV_32S ||
		edges.size() < 3 ||
		(node_neighbours.size() - num_nodes) != 0 ||
		num_nodes < 3)
	{
		fprintf(stderr, "init_tri_node(): input check failed!\n\n");
		return -1;
	}
	int sum = cv::countNonZero(mask);;
	if (sum != num_nodes)
	{
		fprintf(stderr, "init_tri_node(): mask and num_nodes mismatch!\n\n");
		return -1;
	}
	long long num_edges = edges.size();
	node_array.clear();
	node_array.resize(num_nodes);
	int rows = phase.rows;
	int cols = phase.cols;
	int count = 0;
	tri_node* ptr = NULL;
	double Phase;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				Phase = phase.at<double>(i, j);
				ptr = new tri_node(i, j, node_neighbours[count], Phase);
				if (mask.at<int>(i, j) > 1)//邻接已解缠节点标记
				{
					ptr->set_status(true);
				}
				node_array[count] = *ptr;
				delete ptr;
				ptr = NULL;
				count++;

			}
		}
	}

	long* neighbour_ptr = NULL;
	int dummy, ret;
	tri_edge tmp;
	for (int i = 0; i < num_edges; i++)
	{
		tmp = edges[i];
		if (tmp.end1 < 0 ||
			tmp.end2 < 0 ||
			tmp.end1 > node_array.size() ||
			tmp.end2 > node_array.size())
		{
			fprintf(stderr, "init_tri_node(): edges' endpoint exceed legal value!\n\n");
			return -1;
		}
		ret = node_array[tmp.end1 - 1].get_neigh_ptr(&neighbour_ptr, &dummy);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;

		ret = node_array[tmp.end2 - 1].get_neigh_ptr(&neighbour_ptr, &dummy);
		if (return_check(ret, "tri_node::get_neigh_ptr(*, *)", error_head)) return -1;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;
	}
	return 0;
}

int Utils::init_edge_phase_diff(vector<tri_edge>& edges, const vector<tri_node>& node_array)
{
	if (edges.size() < 3 || node_array.size() < 3)
	{
		fprintf(stderr, "init_edge_phase_diff(): input check failed!\n");
		return -1;
	}
	size_t num_nodes = node_array.size();
	size_t num_edges = edges.size();
	size_t end1, end2;
	double phi1, phi2;
	for (size_t i = 0; i < num_edges; i++)
	{
		end1 = edges[i].end1 < edges[i].end2 ? edges[i].end1 : edges[i].end2;
		end2 = edges[i].end1 < edges[i].end2 ? edges[i].end2 : edges[i].end1;
		node_array[end1 - 1].get_phase(&phi1);
		node_array[end2 - 1].get_phase(&phi2);
		phi1 = phi2 - phi1;
		phi1 = atan2(sin(phi1), cos(phi1));
		edges[i].phase_diff = phi1;
	}
	return 0;
}

int Utils::init_edges_quality(Mat& quality, tri_edge* edges, int num_edges, vector<tri_node>& nodes)
{
	if (quality.rows < 2 ||
		quality.cols < 2 ||
		quality.type() != CV_64F ||
		quality.channels() != 1 ||
		edges == NULL ||
		num_edges < 3||
		nodes.size() < 3
		)
	{
		fprintf(stderr, "init_edges_quality(): input check failed!\n\n");
		return -1;
	}
	int rows, cols;
	double qual = 0.0;
	for (int i = 0; i < num_edges; i++)
	{
		qual = 0.0;
		nodes[(edges + i)->end1 - 1].get_pos(&rows, &cols);
		qual += quality.at<double>(rows, cols);
		nodes[(edges + i)->end2 - 1].get_pos(&rows, &cols);
		qual += quality.at<double>(rows, cols);
		(edges + i)->quality = qual / 2.0;
	}
	return 0;
}

int Utils::init_edges_quality(const Mat& quality_map, vector<tri_edge>& edges, const vector<tri_node>& nodes)
{
	if (quality_map.rows < 2 ||
		quality_map.cols < 2 ||
		quality_map.type() != CV_64F ||
		quality_map.channels() != 1 ||
		edges.size() < 3 ||
		nodes.size() < 3
		)
	{
		fprintf(stderr, "init_edges_quality(): input check failed!\n\n");
		return -1;
	}
	int rows, cols;
	double qual = 0.0;
	size_t num_edges = edges.size();
	for (int i = 0; i < num_edges; i++)
	{
		qual = 0.0;
		nodes[edges[i].end1 - 1].get_pos(&rows, &cols);
		qual += quality_map.at<double>(rows, cols);
		nodes[edges[i].end2 - 1].get_pos(&rows, &cols);
		qual += quality_map.at<double>(rows, cols);
		edges[i].quality = qual / 2.0;
	}
	return 0;
}

int Utils::read_triangle(
	const char* ele_file,
	const char* neigh_file,
	triangle** tri,
	int* num_triangle,
	vector<tri_node>& nodes,
	tri_edge* edges,
	int num_edges
)
{
	if (ele_file == NULL ||
		neigh_file == NULL ||
		tri == NULL ||
		num_triangle == NULL||
		nodes.size() < 3||
		edges == NULL||
		num_edges < 3
		)
	{
		fprintf(stderr, "read_triangle(): input check failed!\n\n");
		return -1;
	}
	FILE* fp_ele, *fp_neigh;
	fp_ele = NULL;
	fp_neigh = NULL;
	fp_ele = fopen(ele_file, "rt");
	if (fp_ele == NULL)
	{
		fprintf(stderr, "read_triangle(): can't open %s\n", ele_file);
		return -1;
	}
	fp_neigh = fopen(neigh_file, "rt");
	if (fp_neigh == NULL)
	{
		fprintf(stderr, "read_triangle(): can't open %s\n", neigh_file);
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		return -1;
	}

	char str[INPUTMAXSIZE];
	char* ptr;
	fgets(str, INPUTMAXSIZE, fp_ele);
	*num_triangle = strtol(str, &ptr, 0);
	if (*num_triangle < 1)
	{
		fprintf(stderr, "read_triangle(): number of triangles exceed legal range!\n");
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		if (fp_neigh)
		{
			fclose(fp_neigh);
			fp_neigh = NULL;
		}
		return -1;
	}
	fgets(str, INPUTMAXSIZE, fp_neigh);

	*tri = (triangle*)malloc(*num_triangle * sizeof(triangle));
	if (*tri == NULL)
	{
		fprintf(stderr, "read_triangle(): out of memory!\n");
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		if (fp_neigh)
		{
			fclose(fp_neigh);
			fp_neigh = NULL;
		}
		return -1;
	}
	memset(*tri, 0, sizeof(triangle) * (*num_triangle));

	int p1, p2, p3, neigh1, neigh2, neigh3, num1, num2;
	for (int i = 0; i < *num_triangle; i++)
	{
		fgets(str, INPUTMAXSIZE, fp_ele);
		num1 = strtol(str, &ptr, 0);
		p1 = strtol(ptr, &ptr, 0);
		p2 = strtol(ptr, &ptr, 0);
		p3 = strtol(ptr, &ptr, 0);

		fgets(str, INPUTMAXSIZE, fp_neigh);
		num2 = strtol(str, &ptr, 0);
		neigh1 = strtol(ptr, &ptr, 0);
		neigh2 = strtol(ptr, &ptr, 0);
		neigh3 = strtol(ptr, &ptr, 0);
		(*tri + i)->p1 = p1;
		(*tri + i)->p2 = p2;
		(*tri + i)->p3 = p3;
		(*tri + i)->neigh1 = neigh1;
		(*tri + i)->neigh2 = neigh2;
		(*tri + i)->neigh3 = neigh3;
		(*tri + i)->num = num1;
	}
	if (fp_ele)
	{
		fclose(fp_ele);
		fp_ele = NULL;
	}
	if (fp_neigh)
	{
		fclose(fp_neigh);
		fp_neigh = NULL;
	}
	//获取三角形的边序号
	long* ptr_neigh = NULL;
	int num_neigh, count;
	int edge[3];
	memset(edge, 0, sizeof(int) * 3);
	for (int j = 0; j < *num_triangle; j++)
	{
		count = 0;
		nodes[(*tri + j)->p1 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges + *(ptr_neigh + i) - 1)->end1 == (*tri + j)->p2 ||
				(edges + *(ptr_neigh + i) - 1)->end1 == (*tri + j)->p3 ||
				(edges + *(ptr_neigh + i) - 1)->end2 == (*tri + j)->p2 ||
				(edges + *(ptr_neigh + i) - 1)->end2 == (*tri + j)->p3
				)
			{
				edge[count] = *(ptr_neigh + i);
				count++;
			}
		}
		nodes[(*tri + j)->p2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges + *(ptr_neigh + i) - 1)->end1 == (*tri + j)->p3 ||
				(edges + *(ptr_neigh + i) - 1)->end2 == (*tri + j)->p3)
			{
				edge[count] = *(ptr_neigh + i);
				//count++;
			}
		}
		(*tri + j)->edge1 = edge[0];
		(*tri + j)->edge2 = edge[1];
		(*tri + j)->edge3 = edge[2];
	}
	
	
	
	return 0;
}

int Utils::read_triangle(
	const char* ele_file,
	const char* neigh_file, 
	vector<triangle>& triangle,
	vector<tri_node>& nodes,
	vector<tri_edge>& edges
)
{
	if (ele_file == NULL ||
		neigh_file == NULL ||
		nodes.size() < 3 ||
		edges.size() < 3
		)
	{
		fprintf(stderr, "read_triangle(): input check failed!\n\n");
		return -1;
	}
	FILE* fp_ele, * fp_neigh;
	fp_ele = NULL;
	fp_neigh = NULL;
	fp_ele = fopen(ele_file, "rt");
	if (fp_ele == NULL)
	{
		fprintf(stderr, "read_triangle(): can't open %s\n", ele_file);
		return -1;
	}
	fp_neigh = fopen(neigh_file, "rt");
	if (fp_neigh == NULL)
	{
		fprintf(stderr, "read_triangle(): can't open %s\n", neigh_file);
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		return -1;
	}

	char str[INPUTMAXSIZE];
	char* ptr;
	fgets(str, INPUTMAXSIZE, fp_ele);
	long num_triangle = strtol(str, &ptr, 0);
	if (num_triangle < 1)
	{
		fprintf(stderr, "read_triangle(): number of triangles exceed legal range!\n");
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		if (fp_neigh)
		{
			fclose(fp_neigh);
			fp_neigh = NULL;
		}
		return -1;
	}
	fgets(str, INPUTMAXSIZE, fp_neigh);
	triangle.clear();
	triangle.resize(num_triangle);

	int p1, p2, p3, neigh1, neigh2, neigh3, num1, num2;
	for (int i = 0; i < num_triangle; i++)
	{
		fgets(str, INPUTMAXSIZE, fp_ele);
		num1 = strtol(str, &ptr, 0);
		p1 = strtol(ptr, &ptr, 0);
		p2 = strtol(ptr, &ptr, 0);
		p3 = strtol(ptr, &ptr, 0);

		fgets(str, INPUTMAXSIZE, fp_neigh);
		num2 = strtol(str, &ptr, 0);
		neigh1 = strtol(ptr, &ptr, 0);
		neigh2 = strtol(ptr, &ptr, 0);
		neigh3 = strtol(ptr, &ptr, 0);
		triangle[i].p1 = p1;
		triangle[i].p2 = p2;
		triangle[i].p3 = p3;
		triangle[i].neigh1 = neigh1;
		triangle[i].neigh2 = neigh2;
		triangle[i].neigh3 = neigh3;
		triangle[i].num = num1;
	}
	if (fp_ele)
	{
		fclose(fp_ele);
		fp_ele = NULL;
	}
	if (fp_neigh)
	{
		fclose(fp_neigh);
		fp_neigh = NULL;
	}
	//获取三角形的边序号
	long* ptr_neigh = NULL;
	int num_neigh, count;
	int edge[3];
	memset(edge, 0, sizeof(int) * 3);
	for (int j = 0; j < num_triangle; j++)
	{
		count = 0;
		nodes[triangle[j].p1 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p2) ||
				(edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p3) ||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p2)||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p3)
				)
			{
				edge[count] = *(ptr_neigh + i);
				count++;
			}
		}
		nodes[triangle[j].p2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p3) ||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p3))
			{
				edge[count] = *(ptr_neigh + i);
				//count++;
			}
		}
		triangle[j].edge1 = edge[0];
		triangle[j].edge2 = edge[1];
		triangle[j].edge3 = edge[2];
	}

	return 0;
}

int Utils::gen_delaunay(const char* filename, const char* exe_path)
{
	if (filename == NULL ||
		exe_path == NULL
		)
	{
		fprintf(stderr, "gen_delaunay(): input check failed!\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(filename, "rt");
	if (!fp)
	{
		fprintf(stderr, "gen_delaunay(): can't open %s!\n", filename);
		return -1;
	}
	else
	{
		fclose(fp);
		fp = NULL;
	}
	USES_CONVERSION;
	LPWSTR szCommandLine = new TCHAR[256];
	wcscpy(szCommandLine, A2W(exe_path));
	wcscat(szCommandLine, L"\\delaunay.exe -en ");
	wcscat(szCommandLine, A2W(filename));

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
		HANDLE hd = CreateJobObjectA(NULL, "delaunay");
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
		fprintf(stderr, "gen_triangle(): create triangle.exe process failed!\n\n");
		if (szCommandLine != NULL) delete[] szCommandLine;
		return -1;
	}


	return 0;
}


int Utils::write_node_file(const char* filename, const Mat& mask)
{
	if (filename == NULL ||
		mask.rows < 2 ||
		mask.cols < 2 ||
		mask.channels() != 1 ||
		mask.type() != CV_32S
		)
	{
		fprintf(stderr, "write_node_file(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(filename, "wt");
	if (fp == NULL)
	{
		fprintf(stderr, "write_node_file(): can't open %s\n", filename);
		return -1;
	}
	int nonzero = cv::countNonZero(mask);
	if (nonzero <= 0)
	{
		fprintf(stderr, "write_node_file(): no node exist!\n");
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	int dim, attr1, attr2, count;
	dim = 2;
	attr1 = 0;
	attr2 = 0;
	fprintf(fp, "%d %d %d %d\n", nonzero, dim, attr1, attr2);
	int rows = mask.rows;
	int cols = mask.cols;
	count = 1;
	double x, y;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				x = double(i + 1);
				y = double(j + 1);
				fprintf(fp, "%d %lf %lf\n", count, x, y);
				count++;
			}
		}
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int Utils::PS_amp_dispersion(const vector<Mat>& amplitude, double thresh, Mat& mask)
{
	if (
		amplitude.size() < 3 ||
		thresh < 0.0
		)
	{
		fprintf(stderr, "PS_amp_dispersion(): input check failed!\n\n");
		return -1;
	}
	if (amplitude[0].rows < 1 || amplitude[0].cols < 1 || amplitude[0].channels() != 1 || amplitude[0].type() != CV_64F)
	{
		fprintf(stderr, "PS_amp_dispersion(): input check failed!\n\n");
		return -1;
	}
	int nr = amplitude[0].rows;
	int nc = amplitude[0].cols;
	int num_images = amplitude.size();
	Mat tmp1 = Mat::zeros(nr, nc, CV_64F);
	tmp1.copyTo(mask);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double mean, tmp;
		for (int j = 0; j < nc; j++)
		{
			mean = 0;
			tmp = 0;
			for (int k = 0; k < num_images; k++)
			{
				mean += amplitude[k].at<double>(i, j);
				//mean += (amplitude + k)->at<double>(i, j);
			}
			mean = mean / num_images;
			for (int m = 0; m < num_images; m++)
			{
				tmp += (amplitude[m].at<double>(i, j) - mean) * (amplitude[m].at<double>(i, j) - mean);
				//tmp += ((amplitude + m)->at<double>(i, j) - mean) * ((amplitude + m)->at<double>(i, j) - mean);
			}
			tmp = sqrt(tmp / num_images);
			if (mean < 1e-8) mean = 1e-8;
			tmp = tmp / mean;
			if (tmp <= thresh)
			{
				mask.at<double>(i, j) = 1.0;
			}

		}
	}
	return 0;
}

int Utils::butter_lowpass(int grid_size, int n_win, double low_pass_wavelength, Mat& lowpass)
{
	if (grid_size < 3 ||
		n_win < 3||
		low_pass_wavelength < 1
		)
	{
		fprintf(stderr, "butter_lowpass(): input check failed!\n\n");
		return -1;
	}
	double freq0 = 1 / low_pass_wavelength;
	double start = -(n_win) / grid_size / n_win / 2;
	double interval = 1 / grid_size / n_win;
	double end = (n_win - 2) / grid_size / n_win / 2;
	Mat freq_i = Mat::zeros(1, n_win, CV_64F);
	for (int i = 0; i < n_win; i++)
	{
		freq_i.at<double>(0, i) = start + i * interval;
	}
	freq_i = freq_i / freq0;
	cv::pow(freq_i, 10, freq_i);
	freq_i = freq_i + 1.0;
	freq_i = 1 / freq_i;
	Mat tmp;
	cv::transpose(freq_i, tmp);
	tmp = tmp * freq_i;
	int ret = fftshift2(tmp);
	if (return_check(ret, "fftshift2(*)", error_head)) return -1;
	tmp.copyTo(lowpass);
	//freq_i = -(n_win) / grid_size / n_win / 2:1 / grid_size / n_win : (n_win - 2) / grid_size / n_win / 2;
	//butter_i = 1. / (1 + (freq_i / freq0). ^ (2 * 5));
	//low_pass = butter_i'*butter_i;
	//	low_pass = fftshift(low_pass);
	return 0;
}

int Utils::circshift(Mat& out, const cv::Point& delta)
{
	Size sz = out.size();
	if (sz.height <= 0 ||
		sz.width <= 0
		)
	{
		fprintf(stderr, "circshift(): input check failed!\n\n");
		return -1;
	}
	if ((sz.height == 1 && sz.width == 1)||
		(delta.x == 0 && delta.y == 0)
		)
	{
		return 0;
	}
	int x = delta.x;
	int y = delta.y;
	if (x > 0) x = x % sz.width;
	if (y > 0) y = y % sz.height;
	if (x < 0) x = x % sz.width + sz.width;
	if (y < 0) y = y % sz.height + sz.height;
	vector<Mat> planes;
	split(out, planes);
	for (int i = 0; i < planes.size(); i++)
	{
		Mat tmp0, tmp1, tmp2, tmp3;
		Mat q0(planes[i], Rect(0, 0, sz.width, sz.height - y));
		Mat q1(planes[i], Rect(0, sz.height - y, sz.width, y));
		q0.copyTo(tmp0);
		q1.copyTo(tmp1);
		tmp0.copyTo(planes[i](Rect(0, y, sz.width, sz.height - y)));
		tmp1.copyTo(planes[i](Rect(0, 0, sz.width, y)));


		Mat q2(planes[i], Rect(0, 0, sz.width - x, sz.height));
		Mat q3(planes[i], Rect(sz.width - x, 0, x, sz.height));
		q2.copyTo(tmp2);
		q3.copyTo(tmp3);
		tmp2.copyTo(planes[i](Rect(x, 0, sz.width - x, sz.height)));
		tmp3.copyTo(planes[i](Rect(0, 0, x, sz.height)));
	}
	merge(planes, out);
	return 0;
}

int Utils::fftshift2(Mat& out)
{
	Size sz = out.size();
	Point pt(0, 0);
	pt.x = (int)floor(sz.width / 2.0);
	pt.y = (int)floor(sz.height / 2.0);
	int ret = circshift(out, pt);
	if (return_check(ret, "circshift(*)", error_head)) return -1;
	return 0;
}

int Utils::ifftshift(Mat& out)
{
	Size sz = out.size();
	Point pt(0, 0);
	pt.x = (int)ceil(sz.width / 2.0);
	pt.y = (int)ceil(sz.height / 2.0);
	int ret = circshift(out, pt);
	if (return_check(ret, "circshift(*)", error_head)) return -1;
	return 0;
}

int Utils::fft2(Mat& Src, Mat& Dst)
{
	if (Src.rows < 1 ||
		Src.cols < 1 ||
		Src.channels() != 1 ||
		Src.type() != CV_64F)
	{
		fprintf(stderr, "fft2(): input check failed!\n\n");
		return -1;
	}
	Mat planes[] = { Mat_<double>(Src), Mat::zeros(Src.size(), CV_64F) };
	Mat complexImg;
	merge(planes, 2, complexImg);
	dft(complexImg, Dst, DFT_COMPLEX_OUTPUT);
	return 0;
}

int Utils::fft2(ComplexMat& src, ComplexMat& Dst)
{
	if (src.GetCols() < 1 ||
		src.GetRows() < 1 ||
		src.type() != CV_64F
		)
	{
		fprintf(stderr, "fft2(): input check failed!\n\n");
		return -1;
	}
	Mat re, im;
	src.GetIm().copyTo(im);
	src.GetRe().copyTo(re);
	Mat planes[] = { re, im };
	Mat complexImg;
	merge(planes, 2, complexImg);
	dft(complexImg, complexImg, DFT_COMPLEX_INPUT);
	split(complexImg, planes);
	Dst.SetRe(planes[0]);
	Dst.SetIm(planes[1]);
	return 0;
}

int Utils::ifft2(ComplexMat& src, ComplexMat& dst)
{
	if (src.GetCols() < 1 ||
		src.GetRows() < 1 ||
		src.type() != CV_64F
		)
	{
		fprintf(stderr, "ifft2(): input check failed!\n\n");
		return -1;
	}
	Mat re, im;
	im = src.GetIm();
	re = src.GetRe();
	Mat planes[] = { re, im };
	Mat complexImg;
	merge(planes, 2, complexImg);
	idft(complexImg, complexImg, DFT_COMPLEX_OUTPUT);
	split(complexImg, planes);
	dst.SetRe(planes[0]);
	dst.SetIm(planes[1]);
	return 0;
}

int Utils::std(const Mat& input, double* std)
{
	if (input.rows < 1 ||
		input.cols < 1 ||
		input.channels() != 1 ||
		input.type() != CV_64F ||
		std == NULL
		)
	{
		fprintf(stderr, "std(): input check failed!\n\n");
		return -1;
	}
	Mat data;
	double mean_v = cv::mean(input)[0];
	data = input - mean_v;
	data = data.mul(data);
	double sum_v = cv::sum(data)[0];
	int x = data.rows * data.cols - 1;
	if (x > 0)
	{
		*std = sqrt(sum_v / double(x));
	}
	else
	{
		*std = 0.0;
	}
	
	return 0;
}


int Utils::PS_deflat(
	vector<Mat>& interf_phase,
	Mat& interf_combination,
	vector<Mat>& pos, 
	vector<Mat>& gcps,
	Mat& start_row,
	Mat& start_col, 
	int mode,
	double lambda
)
{
	if (interf_phase.size() < 5 ||
		gcps.size() != pos.size() ||
		gcps.size() < 5 ||
		interf_combination.channels() != 1 ||
		interf_combination.type() != CV_32S ||
		interf_combination.cols != 2 ||
		start_row.rows != start_col.rows ||
		start_row.cols != start_col.cols ||
		//start_row.rows != gcps.size() ||
		start_col.cols != 1 ||
		start_col.channels() != 1 ||
		start_row.channels() != 1 ||
		start_col.type() != CV_32S ||
		start_row.type() != CV_32S ||
		mode < 1 ||
		mode > 2 ||
		lambda <= 0.0
		)
	{
		fprintf(stderr, "PS_deflat(): input check failed!\n\n");
		return -1;
	}
	int n_images = pos.size();
	int n_interf = interf_phase.size();
	if (n_interf != interf_combination.rows ||
		start_row.rows != n_images
		)
	{
		fprintf(stderr, "PS_deflat(): parameter check failed!\n\n");
		return -1;
	}
	volatile bool parallel_flag = true;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n_interf; i++)
	{
		if (!parallel_flag) continue;
		Mat GCPS, POS1, POS2, phase;
		int row_start, col_start, master_ix, slave_ix, ret;
		master_ix = interf_combination.at<int>(i, 0);
		slave_ix = interf_combination.at<int>(i, 1);
		if (master_ix < 1 ||
			master_ix > n_images ||
			slave_ix < 1 ||
			slave_ix > n_images
			)
		{
			fprintf(stderr, "PS_deflat(): image index out of range!\n");
			parallel_flag = false;
			continue;
		}
		GCPS = gcps[master_ix - 1];
		POS1 = pos[master_ix - 1];
		POS2 = pos[slave_ix - 1];
		phase = interf_phase[i];
		row_start = start_row.at<int>(master_ix - 1, 0);
		col_start = start_col.at<int>(slave_ix - 1, 0);
		ret = _PS_deflat(phase, POS1, POS2, GCPS, 1, 1, mode, lambda);
		if (ret < 0)
		{
			parallel_flag = false;
			continue;
		}
		interf_phase[i] = phase;
	}
	if (parallel_check(parallel_flag, "PS_deflat()", parallel_error_head)) return -1;
	return 0;
}

int Utils::_PS_deflat(
	Mat& phase,
	Mat& pos1,
	Mat& pos2,
	Mat& gcps,
	int start_row,
	int start_col,
	int mode,
	double lambda
)
{
	if (
		phase.rows < 1||
		phase.cols < 1||
		phase.channels() != 1||
		phase.type() != CV_64F||
		pos1.rows < phase.rows||
		pos1.cols != 6||
		pos1.channels() != 1||
		pos1.type() != CV_64F||
		pos2.rows < phase.rows ||
		pos2.cols != 6 ||
		pos2.channels() != 1 ||
		pos2.type() != CV_64F ||
		gcps.channels() != 1||
		gcps.type() != CV_64F||
		gcps.cols != 5||
		gcps.rows < 4||
		start_col < 1||
		start_row < 1||
		start_row > pos1.rows||
		mode > 2||
		mode < 1||
		lambda <= 0.0
		)
	{
		fprintf(stderr, "_PS_deflat(): input check failed!\n\n");
		return -1;
	}
	double C = 4.0 * PI;
	if (mode == 1)
	{
		C = 4.0 * PI;
	}
	else
	{
		C = 2.0 * PI;
	}

	int ret;
	int num_gcps = gcps.rows;
	Mat R_M = Mat::zeros(num_gcps, 1, CV_64F);
	Mat R_S = Mat::zeros(num_gcps, 1, CV_64F);
	
	/*计算每个地面控制点的斜距*/
#pragma omp parallel for schedule(guided)
	for (int j = 0; j < num_gcps; j++)
	{
		Mat llh = Mat::zeros(1, 3, CV_64F);
		Mat xyz0, range_sat_tar1, range_sat_tar2;
		Mat fdcs1 = Mat::zeros(1, pos1.rows, CV_64F);
		Mat fdcs2 = Mat::zeros(1, pos2.rows, CV_64F);
		double range_sat_tar_norm1, range_sat_tar_norm2;

		llh.at<double>(0, 0) = gcps.at<double>(j, 0);
		llh.at<double>(0, 1) = gcps.at<double>(j, 1);
		llh.at<double>(0, 2) = gcps.at<double>(j, 2);
		ell2xyz(llh, xyz0);
		Mat xyz1 = Mat::zeros(pos1.rows, 3, CV_64F);
		Mat xyz2 = Mat::zeros(pos2.rows, 3, CV_64F);
		for (int k = 0; k < xyz1.rows; k++)
		{
			xyz0.copyTo(xyz1(Range(k, k + 1), Range(0, 3)));
		}
		for (int k = 0; k < xyz2.rows; k++)
		{
			xyz0.copyTo(xyz2(Range(k, k + 1), Range(0, 3)));
		}
		range_sat_tar1 = xyz1 - pos1(Range(0, pos1.rows), Range(0, 3));
		range_sat_tar2 = xyz2 - pos2(Range(0, pos2.rows), Range(0, 3));
		for (int k = 0; k < range_sat_tar1.rows; k++)
		{
			range_sat_tar1(Range(k, k + 1), Range(0, 3)).copyTo(xyz1);
			range_sat_tar_norm1 = sqrt(sum(xyz1.mul(xyz1))[0]);
			transpose(xyz1, xyz1);
			xyz1 = 2.0 * pos1(Range(k, k + 1), Range(3, 6)) * xyz1 / lambda / (range_sat_tar_norm1 + 1e-12);
			fdcs1.at<double>(0, k) = xyz1.at<double>(0, 0);
		}
		for (int k = 0; k < range_sat_tar2.rows; k++)
		{
			range_sat_tar2(Range(k, k + 1), Range(0, 3)).copyTo(xyz1);
			range_sat_tar_norm2 = sqrt(sum(xyz1.mul(xyz1))[0]);
			transpose(xyz1, xyz1);
			xyz1 = 2.0 * pos2(Range(k, k + 1), Range(3, 6)) * xyz1 / lambda / (range_sat_tar_norm2 + 1e-12);
			fdcs2.at<double>(0, k) = xyz1.at<double>(0, 0);
		}
		Point p1, p2;
		fdcs1 = abs(fdcs1);
		fdcs2 = abs(fdcs2);
		minMaxLoc(fdcs1, NULL, NULL, &p1, NULL);
		minMaxLoc(fdcs2, NULL, NULL, &p2, NULL);
		Mat Sat_pos_M, Sat_pos_S;
		pos1(Range(p1.x, p1.x + 1), Range(0, 3)).copyTo(Sat_pos_M);
		pos2(Range(p2.x, p2.x + 1), Range(0, 3)).copyTo(Sat_pos_S);
		Mat Slantrange_M = xyz0 - Sat_pos_M;
		Mat Slantrange_S = xyz0 - Sat_pos_S;
		R_M.at<double>(j, 0) = sqrt(sum(Slantrange_M.mul(Slantrange_M))[0]);
		R_S.at<double>(j, 0) = sqrt(sum(Slantrange_S.mul(Slantrange_S))[0]);
	}
	Mat A = Mat::zeros(num_gcps, 3, CV_64F);
	for (int j = 0; j < num_gcps; j++)
	{
		A.at<double>(j, 0) = 1.0;
		A.at<double>(j, 1) = gcps.at<double>(j, 3);
		A.at<double>(j, 2) = gcps.at<double>(j, 4);
	}
	Mat A_t;
	transpose(A, A_t);
	A = A_t * A;
	Mat b = A_t * (R_M - R_S);
	Mat x;
	if (!solve(A, b, x, DECOMP_LU))
	{
		fprintf(stderr, "_PS_deflat(): can't solve least square problem!\n");
		return -1;
	}

	Mat flat_phase = Mat::zeros(phase.rows, phase.cols, CV_64F);
	int nr = phase.rows;
	int nc = phase.cols;
#pragma omp parallel for schedule(guided)
	for (int j = 0; j < nr; j++)
	{
		for (int k = 0; k < nc; k++)
		{
			double xx = C * (x.at<double>(0, 0) + (double(j) + double(start_row)) * x.at<double>(1, 0) +
				(double(k) + double(start_col)) * x.at<double>(2, 0)) / lambda;
			phase.at<double>(j, k) = phase.at<double>(j, k) + xx;
		}
	}
	wrap(phase, phase);
	return 0;
}

int Utils::stack_coregistration(
	vector<ComplexMat>& SAR_images,
	Mat& offset,
	int Master_index,
	int coh_method,
	int interp_times,
	int blocksize
)
{
	if (SAR_images.size() < 2||
		Master_index < 1 ||
		Master_index > SAR_images.size() ||
		coh_method < 0 ||
		coh_method > 1 ||
		interp_times < 1||
		blocksize < 4
		)
	{
		fprintf(stderr, "stack_coregistration(): input check failed!\n\n");
		return -1;
	}
	vector<ComplexMat> SAR_images_out;
	SAR_images_out.resize(SAR_images.size());
	Registration coregis;
	int nr = SAR_images[Master_index - 1].GetRows();
	int nc = SAR_images[Master_index - 1].GetCols();
	int n_images = SAR_images.size();
	Mat move_row = Mat::zeros(1, n_images, CV_32S);
	Mat move_col = Mat::zeros(1, n_images, CV_32S);
	Mat Slave_indx = Mat::zeros(1, n_images - 1, CV_32S);
	Mat offset_topleft = Mat::zeros(n_images, 2, CV_32S);
	int count = 0;
	for (int i = 0; i < n_images; i++)
	{
		if (i != Master_index - 1)
		{
			Slave_indx.at<int>(0, count) = i;
			count++;
		}
	}
	////////////粗配准///////////////
	ComplexMat master;
	master = SAR_images[Master_index - 1];
	volatile bool parallel_flag = true;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n_images - 1; i++)
	{
		if (!parallel_flag) continue;
		ComplexMat slave;
		int ret, move_r, move_c;
		int slave_img = Slave_indx.at<int>(0, i);
		slave = SAR_images[slave_img];
		ret = coregis.real_coherent(master, slave, &move_r, &move_c);
		if (ret < 0)
		{
			parallel_flag = false;
			continue;
		}
		move_row.at<int>(0, slave_img) = move_r;
		move_col.at<int>(0, slave_img) = move_c;
	}
	if (parallel_check(parallel_flag, "stack_coregistration()", parallel_error_head)) return -1;
	////////////粗配准公共部分裁剪///////////////
	int move_r_min, move_r_max, move_c_min, move_c_max;
	move_r_min = 100000000;
	move_r_max = -100000000;
	move_c_min = 100000000;
	move_c_max = -100000000;

	int cut_rows;
	for (int i = 0; i < n_images; i++)
	{
		move_r_min = move_r_min < move_row.at<int>(0, i) ? move_r_min : move_row.at<int>(0, i);
		move_c_min = move_c_min < move_col.at<int>(0, i) ? move_c_min : move_col.at<int>(0, i);
		move_r_max = move_r_max > move_row.at<int>(0, i) ? move_r_max : move_row.at<int>(0, i);
		move_c_max = move_c_max > move_col.at<int>(0, i) ? move_c_max : move_col.at<int>(0, i);
	}
	if (abs(move_r_max) >= nr ||
		abs(move_r_min) >= nr ||
		abs(move_c_min) >= nc ||
		abs(move_c_max) >= nc ||
		abs(move_r_max - move_r_min) >= nr||
		abs(move_c_max - move_c_min) >= nc
		)
	{
		fprintf(stderr, "stack_coregistration(): SAR images have no common area!\n\n");
		return -1;
	}
	if (move_r_min >= 0)
	{
		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(0, nr - move_r_max), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = 0;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = move_row.at<int>(0, slave_ix);
			int row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
			int col_start = 0;
			int col_end = SAR_images[slave_ix].GetCols();
			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix);
		}
	}
	if (move_r_max <= 0 && move_r_min < 0)
	{
		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(-move_r_min, nr), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = move_row.at<int>(0, slave_ix) - move_r_min;
			int row_end = nr + move_row.at<int>(0, slave_ix);
			int col_start = 0;
			int col_end = SAR_images[slave_ix].GetCols();
			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
		}
	}
	if (move_r_max > 0 && move_r_min < 0)
	{
		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(-move_r_min, nr - move_r_max), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = move_row.at<int>(0, slave_ix) - move_r_min;
			int row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
			int col_start = 0;
			int col_end = SAR_images[slave_ix].GetCols();
			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
		}
	}


	if (move_c_min >= 0)
	{
		cut_rows = SAR_images_out[Master_index - 1].GetRows();
		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(0, nc - move_c_max));
		offset_topleft.at<int>(Master_index - 1, 1) = 0;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = 0;
			int row_end = SAR_images_out[slave_ix].GetRows();
			int col_start = move_col.at<int>(0, slave_ix);
			int col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix);
		}
	}
	if (move_c_max <= 0 && move_c_min < 0)
	{
		cut_rows = SAR_images_out[Master_index - 1].GetRows();
		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(-move_c_min, nc));
		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = 0;
			int row_end = SAR_images_out[slave_ix].GetRows();
			int col_start = move_col.at<int>(0, slave_ix) - move_c_min;
			int col_end = nc + move_col.at<int>(0, slave_ix);
			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
		}
	}
	if (move_c_max > 0 && move_c_min < 0)
	{
		cut_rows = SAR_images_out[Master_index - 1].GetRows();
		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(-move_c_min, nc - move_c_max));
		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_images - 1; i++)
		{
			int slave_ix = Slave_indx.at<int>(0, i);
			int row_start = 0;
			int row_end = SAR_images_out[slave_ix].GetRows();
			int col_start = move_col.at<int>(0, slave_ix) - move_c_min;
			int col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
		}
	}
	offset_topleft.copyTo(offset);

	////////////////////////精配准////////////////////////////
	SAR_images[Master_index - 1] = SAR_images_out[Master_index - 1];
	master = SAR_images_out[Master_index - 1];
	//这里并行加速，因为子函数已经进行了加速
//#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n_images - 1; i++)
	{
		ComplexMat slave;
		int slave_ix = Slave_indx.at<int>(0, i);
		slave = SAR_images_out[slave_ix];
		int ret;
		ret = coregis.coregistration_subpixel(master, slave, blocksize, interp_times);
		if (return_check(ret, "registration_subpixel(*, *, *, *)", error_head)) return -1;
		SAR_images[slave_ix] = slave;
		//cout << i + 1 << "/" << n_images - 1 << "\n";
	}
	return 0;
}

int Utils::stack_coregistration(vector<string>& SAR_images, vector<string>& SAR_images_out, Mat& offset, int Master_index, int interp_times, int blocksize)
{
	if (SAR_images.size() < 2 ||
		Master_index < 1 ||
		Master_index > SAR_images.size() ||
		SAR_images_out.size() != SAR_images.size() ||
		interp_times < 1 ||
		blocksize < 16
		)
	{
		fprintf(stderr, "stack_coregistration(): input check failed!\n\n");
		return -1;
	}

	Registration coregis; FormatConversion conversion;
	int n_images = SAR_images.size(), ret;
	Mat move_row = Mat::zeros(1, n_images, CV_32S);
	Mat move_col = Mat::zeros(1, n_images, CV_32S);
	Mat Slave_indx = Mat::zeros(1, n_images - 1, CV_32S);
	Mat offset_topleft = Mat::zeros(n_images, 2, CV_32S);
	int count = 0;
	//创建配准后输出的h5文件
	for (int i = 0; i < n_images; i++)
	{
		ret = conversion.creat_new_h5(SAR_images_out[i].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	}
	for (int i = 0; i < n_images; i++)
	{
		if (i != Master_index - 1)
		{
			Slave_indx.at<int>(0, count) = i;
			count++;
		}
	}
	//检查并确保SAR图像尺寸大小一致
	int min_row = 10000000, min_col = 10000000;
	Mat azimuth_len, range_len;
	for (int i = 0; i < n_images; i++)
	{
		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "azimuth_len", azimuth_len);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "range_len", range_len);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		min_row = azimuth_len.at<int>(0, 0) > min_row ? min_row : azimuth_len.at<int>(0, 0);
		min_col = range_len.at<int>(0, 0) > min_col ? min_col : range_len.at<int>(0, 0);
	}
	if (min_col < 1 || min_row < 1)
	{
		fprintf(stderr, "invalide SAR images size\n");
		return -1;
	}
	
	////////////粗配准///////////////
	ComplexMat master, slave, master_out, slave_out;
	ret = conversion.read_slc_from_h5(SAR_images[Master_index - 1].c_str(), master);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	master = master(cv::Range(0, min_row), cv::Range(0, min_col));
	if (master.type() != CV_64F) master.convertTo(master, CV_64F);
	int nr = master.GetRows();
	int nc = master.GetCols();
	for (int i = 0; i < n_images - 1; i++)
	{
		
		int move_r, move_c;
		int slave_img = Slave_indx.at<int>(0, i);
		ret = conversion.read_slc_from_h5(SAR_images[slave_img].c_str(), slave);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		slave = slave(cv::Range(0, min_row), cv::Range(0, min_col));
		if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
		ret = coregis.real_coherent(master, slave, &move_r, &move_c);
		if (return_check(ret, "real_coherent()", error_head)) return -1;
		move_row.at<int>(0, slave_img) = move_r;
		move_col.at<int>(0, slave_img) = move_c;
	}
	////////////粗配准公共部分裁剪///////////////
	int move_r_min, move_r_max, move_c_min, move_c_max;
	move_r_min = 100000000;
	move_r_max = -100000000;
	move_c_min = 100000000;
	move_c_max = -100000000;

	int cut_rows;
	for (int i = 0; i < n_images; i++)
	{
		move_r_min = move_r_min < move_row.at<int>(0, i) ? move_r_min : move_row.at<int>(0, i);
		move_c_min = move_c_min < move_col.at<int>(0, i) ? move_c_min : move_col.at<int>(0, i);
		move_r_max = move_r_max > move_row.at<int>(0, i) ? move_r_max : move_row.at<int>(0, i);
		move_c_max = move_c_max > move_col.at<int>(0, i) ? move_c_max : move_col.at<int>(0, i);
	}
	if (abs(move_r_max) >= nr ||
		abs(move_r_min) >= nr ||
		abs(move_c_min) >= nc ||
		abs(move_c_max) >= nc ||
		abs(move_r_max - move_r_min) >= nr ||
		abs(move_c_max - move_c_min) >= nc
		)
	{
		fprintf(stderr, "stack_coregistration(): SAR images have no common area!\n\n");
		return -1;
	}

	//主图像裁剪
	if (move_r_min >= 0)
	{
		master = master(Range(0, nr - move_r_max), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = 0;
	}
	if (move_r_max <= 0 && move_r_min < 0)
	{
		master = master(Range(-move_r_min, nr), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
	}
	if (move_r_max > 0 && move_r_min < 0)
	{
		master = master(Range(-move_r_min, nr - move_r_max), Range(0, nc));
		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
	}
	if (move_c_min >= 0)
	{
		cut_rows = master.GetRows();
		master = master(Range(0, cut_rows), Range(0, nc - move_c_max));
		offset_topleft.at<int>(Master_index - 1, 1) = 0;
	}
	if (move_c_max <= 0 && move_c_min < 0)
	{
		cut_rows = master.GetRows();
		master = master(Range(0, cut_rows), Range(-move_c_min, nc));
		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
	}
	if (move_c_max > 0 && move_c_min < 0)
	{
		cut_rows = master.GetRows();
		master = master(Range(0, cut_rows), Range(-move_c_min, nc - move_c_max));
		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
	}
	ret = conversion.write_slc_to_h5(SAR_images_out[Master_index - 1].c_str(), master);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	azimuth_len.at<int>(0, 0) = master.GetRows();
	range_len.at<int>(0, 0) = master.GetCols();
	ret = conversion.write_array_to_h5(SAR_images_out[Master_index - 1].c_str(), "azimuth_len", azimuth_len);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = conversion.write_array_to_h5(SAR_images_out[Master_index - 1].c_str(), "range_len", range_len);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//辅图像裁剪
	for (int i = 0; i < n_images - 1; i++)
	{
		int slave_ix, row_start, row_end, col_start, col_end;
		slave_ix = Slave_indx.at<int>(0, i);
		ret = conversion.read_slc_from_h5(SAR_images[slave_ix].c_str(), slave);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
		if (move_r_min >= 0)
		{
			
			row_start = move_row.at<int>(0, slave_ix);
			row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
			col_start = 0;
			col_end = slave.GetCols();
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix);
		}
		if (move_r_max <= 0 && move_r_min < 0)
		{
			row_start = move_row.at<int>(0, slave_ix) - move_r_min;
			row_end = nr + move_row.at<int>(0, slave_ix);
			col_start = 0;
			col_end = slave.GetCols();
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
		}
		if (move_r_max > 0 && move_r_min < 0)
		{
			row_start = move_row.at<int>(0, slave_ix) - move_r_min;
			row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
			col_start = 0;
			col_end = slave.GetCols();
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
		}


		if (move_c_min >= 0)
		{
			row_start = 0;
			row_end = slave.GetRows();
			col_start = move_col.at<int>(0, slave_ix);
			col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix);
		}
		if (move_c_max <= 0 && move_c_min < 0)
		{
			row_start = 0;
			row_end = slave.GetRows();
			col_start = move_col.at<int>(0, slave_ix) - move_c_min;
			col_end = nc + move_col.at<int>(0, slave_ix);
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
		}
		if (move_c_max > 0 && move_c_min < 0)
		{
			row_start = 0;
			row_end = slave.GetRows();
			col_start = move_col.at<int>(0, slave_ix) - move_c_min;
			col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
		}

		//精配准
		ret = coregis.coregistration_subpixel(master, slave, blocksize, interp_times);
		if (return_check(ret, "coregistration_subpixel()", error_head)) return -1;
		//写出
		ret = conversion.write_slc_to_h5(SAR_images_out[slave_ix].c_str(), slave);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		azimuth_len.at<int>(0, 0) = slave.GetRows();
		range_len.at<int>(0, 0) = slave.GetCols();
		ret = conversion.write_array_to_h5(SAR_images_out[slave_ix].c_str(), "azimuth_len", azimuth_len);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(SAR_images_out[slave_ix].c_str(), "range_len", range_len);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}
	offset_topleft.copyTo(offset);


	return 0;
}

int Utils::PS_gen_interferogram(
	vector<ComplexMat>& SAR_images,
	vector<Mat>& interf_phase,
	Mat& time_baseline, 
	Mat& spatial_baseline,
	double time_thresh,
	double spatial_thresh,
	int Mast_index,
	Mat& interf_combination
)
{
	int n_images = SAR_images.size();
	if (SAR_images.size() < 5 ||
		time_thresh < 0.0 ||
		spatial_thresh < 0.0||
		Mast_index > n_images||
		time_baseline.rows != time_baseline.cols||
		time_baseline.type() != CV_64F||
		time_baseline.channels() != 1||
		spatial_baseline.rows != spatial_baseline.cols ||
		spatial_baseline.type() != CV_64F ||
		spatial_baseline.channels() != 1 ||
		spatial_baseline.rows != time_baseline.rows||
		n_images != spatial_baseline.rows
		)
	{
		fprintf(stderr, "PS_gen_interferogram(): input check failed!\n\n");
		return -1;
	}
	vector<int> master;
	vector<int> slave;
	/*自由组合干涉*/
	if (Mast_index <= 0)
	{
		vector<double> temporal_Baseline;
		vector<double> spatial_Baseline;
		Mat phase;
		for (int i = 0; i < n_images; i++)
		{
			for (int j = i + 1; j < n_images; j++)
			{
				if (fabs(time_baseline.at<double>(i, j)) <= time_thresh &&
					fabs(spatial_baseline.at<double>(i, j)) < spatial_thresh
					)
				{
					temporal_Baseline.push_back(time_baseline.at<double>(i, j));
					spatial_Baseline.push_back(spatial_baseline.at<double>(i, j));
					phase = (SAR_images[i] * (SAR_images[j].conj())).GetPhase();
					interf_phase.push_back(phase);
					master.push_back(i + 1);
					slave.push_back(j + 1);
				}
			}
		}
		Mat t, s;
		t = Mat::zeros(1, temporal_Baseline.size(), CV_64F);
		s = Mat::zeros(1, temporal_Baseline.size(), CV_64F);
		for (int i = 0; i < temporal_Baseline.size(); i++)
		{
			t.at<double>(0, i) = temporal_Baseline[i];
			s.at<double>(0, i) = spatial_Baseline[i];
		}
		t.copyTo(time_baseline);
		s.copyTo(spatial_baseline);
	}
	/*公共主图像干涉*/
	else
	{
		vector<double> temporal_Baseline;
		vector<double> spatial_Baseline;
		Mat phase;
		for (int i = 0; i < n_images; i++)
		{
			if (i != Mast_index - 1)
			{
				temporal_Baseline.push_back(time_baseline.at<double>(Mast_index - 1, i));
				spatial_Baseline.push_back(spatial_baseline.at<double>(Mast_index - 1, i));
				phase = (SAR_images[Mast_index - 1] * (SAR_images[i].conj())).GetPhase();
				interf_phase.push_back(phase);
				master.push_back(Mast_index);
				slave.push_back(i + 1);
			}
		}
		Mat t, s;
		t = Mat::zeros(1, temporal_Baseline.size(), CV_64F);
		s = Mat::zeros(1, temporal_Baseline.size(), CV_64F);
		for (int i = 0; i < temporal_Baseline.size(); i++)
		{
			t.at<double>(0, i) = temporal_Baseline[i];
			s.at<double>(0, i) = spatial_Baseline[i];
		}
		t.copyTo(time_baseline);
		s.copyTo(spatial_baseline);
	}
	Mat combin = Mat::zeros(master.size(), 2, CV_32S);
	for (int i = 0; i < combin.rows; i++)
	{
		combin.at<int>(i, 0) = master[i];
		combin.at<int>(i, 1) = slave[i];
	}
	combin.copyTo(interf_combination);
	return 0;
}

int Utils::PS_spatialbaseline_est(vector<Mat>& pos, Mat gcps, double lambda, Mat& MB_effect)
{
	if (pos.size() < 3||
		gcps.rows < 1 ||
		gcps.cols != 3 ||
		gcps.channels() != 1 ||
		gcps.type() != CV_64F||
		lambda <= 0.0
		)
	{
		fprintf(stderr, "PS_spatialbaseline_est(): input check failed!\n\n");
		return -1;
	}
	int n_images = pos.size();
	int Num_gcp = gcps.rows;
	Mat MB_effect_tmp = Mat::zeros(n_images, n_images, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n_images; i++)
	{
		int Mast_ind = i;
		int count = 0;
		Mat Slave_ind = Mat::zeros(1, n_images - 1, CV_32S);
		for (int j = 0; j < n_images; j++)
		{
			if (j != i)
			{
				Slave_ind.at<int>(0, count) = j;
				count++;
			}
		}
		Mat VecR_S_T;
		Mat B_effect = Mat::zeros(Num_gcp, n_images, CV_64F);
		Mat ms_theta = Mat::zeros(Num_gcp, n_images, CV_64F);
		Mat Baselinespan = Mat::zeros(Num_gcp, n_images - 1, CV_64F);
		int rows = pos[Mast_ind].rows;
		for (int j = 0; j < Num_gcp; j++)
		{
			///////////////////计算主星成像位置///////////////////////
			Mat xyz, xyz0;
			Mat llh = Mat::zeros(1, 3, CV_64F);
			llh.at<double>(0, 0) = gcps.at<double>(j, 0);
			llh.at<double>(0, 1) = gcps.at<double>(j, 1);
			llh.at<double>(0, 2) = gcps.at<double>(j, 2);
			ell2xyz(llh, xyz);
			xyz.copyTo(xyz0);
			//cout << xyz0 << "\n";
			Mat xyz_rep = Mat::zeros(rows, 3, CV_64F);
			for (int k = 0; k < rows; k++)
			{
				xyz.copyTo(xyz_rep(Range(k, k + 1), Range(0, 3)));
			}
			Mat fdcs = Mat::zeros(1, rows, CV_64F);
			Mat range_sat_tar;
			range_sat_tar = xyz_rep - pos[Mast_ind](Range(0, rows), Range(0, 3));
			double norm_range_sat_tar;
			for (int k = 0; k < rows; k++)
			{
				range_sat_tar(Range(k, k + 1), Range(0, 3)).copyTo(xyz);
				norm_range_sat_tar = sqrt(sum(xyz.mul(xyz))[0]);
				transpose(xyz, xyz);
				xyz = 2 * pos[Mast_ind](Range(k, k + 1), Range(3, 6)) * xyz / lambda / (norm_range_sat_tar + 1e-9);
				fdcs.at<double>(0, k) = xyz.at<double>(0, 0);
			}
			Point p;
			fdcs = abs(fdcs);
			minMaxLoc(fdcs, NULL, NULL, &p, NULL);
			Mat M_Sat_pos;
			pos[Mast_ind](Range(p.x, p.x + 1), Range(0, 3)).copyTo(M_Sat_pos);
			//cout << M_Sat_pos << "\n";
			VecR_S_T = xyz0 - M_Sat_pos;
			Mat Sat_velocity;
			pos[Mast_ind](Range(p.x, p.x + 1), Range(3, 6)).copyTo(Sat_velocity);
			//cout << Sat_velocity << "\n";
			Mat effect_dir;
			cross(VecR_S_T, Sat_velocity, effect_dir);
			//cout << effect_dir << "\n";
			double norm_effect_dir = sqrt(sum(effect_dir.mul(effect_dir))[0]);
			if (sum(effect_dir.mul(xyz0))[0] < 0.0)
			{
				effect_dir = -effect_dir / (norm_effect_dir + 1e-10);
			}
			else
			{
				effect_dir = effect_dir / (norm_effect_dir + 1e-10);
			}

			///////////////////计算辅星成像位置///////////////////////
			int rows_slave;
			int slave_images_idx;
			for (int k = 0; k < n_images - 1; k++)
			{
				slave_images_idx = Slave_ind.at<int>(0, k);
				rows_slave = pos[slave_images_idx].rows;
				fdcs = Mat::zeros(1, rows_slave, CV_64F);
				xyz_rep = Mat::zeros(rows_slave, 3, CV_64F);
				for (int kk = 0; kk < rows_slave; kk++)
				{
					xyz0.copyTo(xyz_rep(Range(kk, kk + 1), Range(0, 3)));
				}
				range_sat_tar = xyz_rep - pos[slave_images_idx](Range(0, rows_slave), Range(0, 3));
				for (int kk = 0; kk < rows_slave; kk++)
				{
					range_sat_tar(Range(kk, kk + 1), Range(0, 3)).copyTo(xyz);
					norm_range_sat_tar = sqrt(sum(xyz.mul(xyz))[0]);
					transpose(xyz, xyz);
					xyz = 2 * pos[slave_images_idx](Range(kk, kk + 1), Range(3, 6)) * xyz / lambda / (norm_range_sat_tar + 1e-9);
					fdcs.at<double>(0, kk) = xyz.at<double>(0, 0);
				}
				fdcs = abs(fdcs);
				minMaxLoc(fdcs, NULL, NULL, &p, NULL);
				Mat S_Sat_pos;
				pos[slave_images_idx](Range(p.x, p.x + 1), Range(0, 3)).copyTo(S_Sat_pos);
				Mat Rang_S_T = xyz0 - S_Sat_pos;
				Mat Baselines = S_Sat_pos - M_Sat_pos;
				double effect_symb;
				if (sum(Baselines.mul(effect_dir))[0] < 0.0) effect_symb = -1.0;
				else
				{
					effect_symb = 1.0;
				}
				Baselinespan.at<double>(j, k) = sqrt(sum(Baselines.mul(Baselines))[0]);

				double norm_VecR_S_T = sqrt(sum(VecR_S_T.mul(VecR_S_T))[0]);
				double norm_Rang_S_T = sqrt(sum(Rang_S_T.mul(Rang_S_T))[0]);
				ms_theta.at<double>(j, slave_images_idx) = acos((sum(Baselines.mul(Baselines))[0] +
					norm_VecR_S_T * norm_VecR_S_T - norm_Rang_S_T * norm_Rang_S_T) /
					(2 * Baselinespan.at<double>(j, k) * norm_VecR_S_T + 1e-10));
				B_effect.at<double>(j, slave_images_idx) = effect_symb * Baselinespan.at<double>(j, k) * 
					sin(ms_theta.at<double>(j, slave_images_idx));

			}

		}
		for (int j = 0; j < n_images; j++)
		{
			MB_effect_tmp.at<double>(i, j) = sum(B_effect(Range(0, Num_gcp), Range(j, j + 1)))[0] / double(Num_gcp);
		}
	}
	MB_effect_tmp.copyTo(MB_effect);
	return 0;
}

int Utils::Clap(ComplexMat& phase, ComplexMat& phase_filter, double alpha, double beta, int n_win, int n_pad, Mat& lowpass)
{
	if (phase.GetCols() < 3 ||
		phase.GetRows() < 3 ||
		phase.type() != CV_64F ||
		alpha <= 0 ||
		beta < 0 ||
		n_win < 5 ||
		n_pad < 0 ||
		lowpass.cols < 3 ||
		lowpass.rows < 3 ||
		lowpass.channels() != 1 ||
		lowpass.type() != CV_64F)
	{
		fprintf(stderr, "Clap(): input check failed!\n\n");
		return -1;
	}
	if (lowpass.cols != (n_win + n_pad) ||
		lowpass.rows != (n_win + n_pad))
	{
		fprintf(stderr, "Clap(): n_win + n_pad doesn't match lowpass.size, please check!\n\n");
		return -1;
	}
	if (n_win % 2 == 1)
	{
		n_win -= 1;
		n_pad = lowpass.rows - n_win;
	}
	int n_i = phase.GetRows();
	int n_j = phase.GetCols();
	ComplexMat ph;
	//Mat cos, sin;
	int ret;
	//ret = phase2cos(phase, cos, sin);
	//if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	//ph.SetIm(sin);
	//ph.SetRe(cos);
	ph = phase;
	ComplexMat ph_out(n_i, n_j);
	int n_inc = floor(n_win / 4);
	int n_win_i = ceil(n_i / n_inc) - 3;
	int n_win_j = ceil(n_j / n_inc) - 3;
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
				fprintf(stderr, "Clap(): wf2.size and n_win mismatch, please check to make sure n_win is even!\n\n");
				return -1;
			}
			temp = ph(cv::Range(i1 - 1, i2), cv::Range(j1 - 1, j2));
			ret = ph_bit.SetValue(cv::Range(0, n_win), cv::Range(0, n_win), temp);
			if (return_check(ret, "ComplexMat::SetValue(*, *, *)", error_head)) return -1;
			ret = fft2(ph_bit, fft_out);
			if (return_check(ret, "fft2(*, *)", error_head)) return -1;
			H = fft_out.GetMod();
			ret = fftshift2(H);
			
			if (return_check(ret, "fftshift2(*, *)", error_head)) return -1;
			filter2D(H, H, -1, guasswin, Point(-1, -1), 0.0, BORDER_CONSTANT);
			
			ret = ifftshift(H);
			
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
			H = H - 1;
			for (int i = 0; i < H.rows; i++)
			{
				for (int j = 0; j < H.cols; j++)
				{
					if (H.at<double>(i, j) < 0.0)
					{
						H.at<double>(i, j) = 0.0;
					}
				}
			}
			H = H * beta + lowpass;
			
			fft_out = fft_out * H;
			ret = ifft2(fft_out, temp);
			if (return_check(ret, "ifft2(*, *)", error_head)) return -1;
			Mat re, im;
			re = temp.GetRe();
			im = temp.GetIm();
			temp1 = temp(Range(0, n_win), Range(0, n_win));
			ph_filt =  temp1 * wf2;
			temp = ph_out(Range(i1 - 1, i2), Range(j1 - 1, j2)) + ph_filt;
			ret = ph_out.SetValue(Range(i1 - 1, i2), Range(j1 - 1, j2), temp);
			if (return_check(ret, "ph_out.SetValue(*, *, *)", error_head)) return -1;
		}

	}
	phase_filter = ph_out;
	return 0;
}

int Utils::PS_topofit(Mat& dif_phase, Mat& bperp0, double num_trial_wraps, double* K0, double* C0, double* coh0, Mat& phase_residue)
{
	if (dif_phase.cols != 1 ||
		dif_phase.rows != bperp0.rows ||
		bperp0.cols != 1 ||
		dif_phase.rows < 3 ||
		dif_phase.type() != CV_64F||
		dif_phase.channels() != 1||
		bperp0.channels() != 1||
		bperp0.type() != CV_64F||
		num_trial_wraps < 0.0||
		K0 == NULL||
		C0 == NULL||
		coh0 == NULL
		)
	{
		fprintf(stderr, "PS_topofit(): input check failed!\n\n");
		return -1;
	}
	ComplexMat cpxphase, trial_phase_mat, phaser;
	Mat re, im, trial_phase;
	int ret, nonzero, count, n_trials, lower;
	double min, max, bperp_range;
	nonzero = 0;
	for (int i = 0; i < dif_phase.rows; i++)
	{
		if (fabs(dif_phase.at<double>(i, 0)) > DBL_EPSILON) nonzero++;
	}
	if (nonzero == 0)
	{
		fprintf(stderr, "PS_topofit(): dif_phase is a zero matrix!\n");
		return -1;
	}
	Mat phase = Mat::zeros(nonzero, 1, CV_64F);
	Mat bperp = Mat::zeros(nonzero, 1, CV_64F);
	count = 0;
	for (int i = 0; i < dif_phase.rows; i++)
	{
		if (fabs(dif_phase.at<double>(i, 0)) > DBL_EPSILON)
		{
			phase.at<double>(count, 0) = dif_phase.at<double>(i, 0);
			bperp.at<double>(count, 0) = bperp0.at<double>(i, 0);
			count++;
		}
		
	}
	ret = phase2cos(phase, re, im);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	cpxphase.SetRe(re);
	cpxphase.SetIm(im);
	minMaxLoc(bperp, &min, &max);
	bperp_range = max - min;
	if (fabs(bperp_range) <= 1e-10)
	{
		fprintf(stderr, "PS_topofit(): bperp_range too small!\n");
		return -1;
	}
	lower = int(ceil(num_trial_wraps * 8));
	n_trials = lower * 2 + 1;
	Mat trial_mult = Mat::zeros(1, n_trials, CV_64F);
	for (int i = 0; i < n_trials; i++)
	{
		trial_mult.at<double>(0, i) = double(i - lower);
	}
	trial_phase = bperp / bperp_range * PI / 4;
	trial_phase = 0 - trial_phase * trial_mult;
	ret = phase2cos(trial_phase, re, im);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	trial_phase_mat.SetRe(re);
	trial_phase_mat.SetIm(im);
	ComplexMat cpxphase_mat(cpxphase.GetRows(), n_trials);
	for (int i = 0; i < n_trials; i++)
	{
		ret = cpxphase_mat.SetValue(Range(0, cpxphase_mat.GetRows()), Range(i, i + 1), cpxphase);
		if (return_check(ret, "cpxphase_mat.SetValue(*, *, *)", error_head)) return -1;
	}
	ret = trial_phase_mat.Mul(cpxphase_mat, phaser, false);
	if (return_check(ret, "ComplexMat::Mul(*, *, *)", error_head)) return -1;
	ComplexMat phaser_sum = phaser.sum();
	Mat C_trial = phaser_sum.GetPhase();
	Mat coh_trial = phaser_sum.GetMod() / nonzero;
	Point coh_high_max_ix;
	minMaxLoc(coh_trial, NULL, NULL, NULL, &coh_high_max_ix);
	*K0 = PI / 4 / bperp_range * trial_mult.at<double>(0, coh_high_max_ix.x);
	*C0 = C_trial.at<double>(0, coh_high_max_ix.x);
	*coh0 = coh_trial.at<double>(0, coh_high_max_ix.x);

	ComplexMat resphase, offset_phase;
	trial_phase = -(*K0) * bperp;
	ret = phase2cos(trial_phase, re, im);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	resphase.SetRe(re);
	resphase.SetIm(im);
	ret = cpxphase.Mul(resphase, resphase, false);
	if (return_check(ret, "ComplexMat::Mul(*, *, *)", error_head)) return -1;
	offset_phase = resphase.sum();
	offset_phase = offset_phase.conj();
	resphase = resphase * offset_phase;
	Mat resphase_angle;
	resphase_angle = resphase.GetPhase();
	Mat A, b;
	bperp.copyTo(A);
	resphase_angle.copyTo(b);
	double a = (A.dot(A));
	if (a < 1e-9) a = 1e-9;
	//transpose(A, A);
	a = 1 / a * (A.dot(b));
	*K0 += a;
	trial_phase = -(*K0) * bperp;
	ret = phase2cos(trial_phase, re, im);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	resphase.SetRe(re);
	resphase.SetIm(im);
	ret = cpxphase.Mul(resphase, resphase, false);
	if (return_check(ret, "ComplexMat::Mul(*, *, *)", error_head)) return -1;
	resphase.GetPhase().copyTo(phase_residue);
	resphase = resphase.sum();
	*C0 = resphase.GetPhase().at<double>(0, 0);
	*coh0 = resphase.GetMod().at<double>(0, 0) / nonzero;
	return 0;
}

int Utils::hist(Mat& input, double lowercenter, double uppercenter, double interval, Mat& out)
{
	if (input.rows != 1 ||
		input.cols < 3 ||
		input.channels() != 1 ||
		input.type() != CV_64F ||
		interval >(uppercenter - lowercenter) ||
		lowercenter > uppercenter ||
		interval <= 0.0
		)
	{
		fprintf(stderr, "hist(): input check failed!\n\n");
		return -1;
	}
	Mat tmp;
	input.copyTo(tmp);
	cv::sort(tmp, tmp, SORT_ASCENDING + SORT_EVERY_ROW);
	int len = int(floor((uppercenter - lowercenter) / interval)) + 1;
	int nc = input.cols;
	Mat H = Mat::zeros(1, len, CV_64F);
	double cmp = lowercenter + interval*0.5;
	int count = 0;
	int idx = 0;
	int k;
	for (k = 0; k < nc; k++)
	{
		if (tmp.at<double>(0, k) < cmp) count++;
		else
		{
			H.at<double>(0, idx) = double(count);
			idx++;
			count = 1;
			cmp += interval;
			if (idx >= len - 1) break;
		}
	}
	H.at<double>(0, len - 1) = H.at<double>(0, len - 1) + double(nc - k);
	H.copyTo(out);
	return 0;
}

int Utils::llh2local(const Mat& llh0, double lon0, double lat0, Mat& xy)
{
	if (llh0.cols < 1 ||
		llh0.rows != 2 ||
		llh0.type() != CV_64F ||
		llh0.channels() != 1 ||
		lon0 > 180.0 ||
		lon0 < 0.0 ||
		lat0 > 90.0 ||
		lat0 < 0.0
		)
	{
		fprintf(stderr, "llh2local(): input check failed!\n\n");
		return -1;
	}
	int ret;
	double a = 6378137.0;
	double e = 0.08209443794970;
	Mat llh = llh0 / 180.0 * PI;
	lon0 = lon0 / 180 * PI;
	lat0 = lat0 / 180 * PI;
	Mat dlambda = llh(Range(0, 1), Range(0, llh.cols)) - lon0;
	double const1, const2, const3, const4, M0;
	Mat M, N, E, tmp1, tmp2, tmp3;
	M = Mat::zeros(1, llh.cols, CV_64F);
	const1 = 1 - e * e / 4 - 3 * e * e * e * e / 64 - 5 * e * e * e * e * e * e / 256;
	const2 = 3 * e * e / 8 + 3 * e * e * e * e / 32 + 45 * e * e * e * e * e * e / 1024;
	const3 = 15 * e * e * e * e / 256 + 45 * e * e * e * e * e * e / 1024;
	const4 = 35 * e * e * e * e * e * e / 3072;
	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	tmp1 = tmp1 * const1;
	tmp1.copyTo(M);

	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	tmp1 = tmp1 * 2;
	ret = phase2cos(tmp1, tmp2, tmp3);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	M = M - const2 * tmp3;

	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	tmp1 = tmp1 * 4;
	ret = phase2cos(tmp1, tmp2, tmp3);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	M = M + const3 * tmp3;

	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	tmp1 = tmp1 * 6;
	ret = phase2cos(tmp1, tmp2, tmp3);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	M = M - const4 * tmp3;

	M = M * a;

	M0 = a * (const1 * lat0 - const2 * sin(2 * lat0) + const3 * sin(4 * lat0) - const4 * sin(6 * lat0));

	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	ret = phase2cos(tmp1, tmp2, tmp1);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	tmp1 = tmp1.mul(tmp1);
	tmp1 = tmp1 * e * e;
	tmp1 = 1 - tmp1;
	cv::pow(tmp1, 1.0 / 2.0, tmp1);
	N = a / tmp1;

	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	ret = phase2cos(tmp1, tmp2, tmp1);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	E = dlambda.mul(tmp1);

	Mat xy1 = Mat::zeros(llh.rows, llh.cols, CV_64F);
	llh(Range(1, 2), Range(0, llh.cols)).copyTo(tmp1);
	for (int i = 0; i < llh.cols; i++)
	{
		if (fabs(llh.at<double>(1, i)) > 1e-10)
		{
			xy1.at<double>(0, i) = N.at<double>(0, i) * cos(tmp1.at<double>(0, i)) 
				/ sin(tmp1.at<double>(0, i)) * sin(E.at<double>(0, i));
			xy1.at<double>(1, i) = M.at<double>(0, i) - M0 + N.at<double>(0, i) * cos(tmp1.at<double>(0, i))
				/ sin(tmp1.at<double>(0, i)) * (1.0 - cos(E.at<double>(0, i)));
		}
		else
		{
			xy1.at<double>(0, i) = a * dlambda.at<double>(0, i);
			xy1.at<double>(1, i) = -M0;
		}
	}
	xy1 = xy1 / 1000;
	xy1.copyTo(xy);
	return 0;
}

int Utils::gen_grid(const Mat& xy, double grid_size, Mat& grid_ij, int* rows, int* cols)
{
	if (xy.rows < 3 ||
		xy.cols != 3 ||
		xy.type() != CV_64F ||
		xy.channels() != 1 ||
		grid_size < 1 ||
		rows == NULL ||
		cols == NULL
		)
	{
		fprintf(stderr, "gen_grid(): input check failed!\n\n");
		return -1;
	}
	int nr = xy.rows;
	double min2, min3, max2, max3;
	Mat grid = Mat::zeros(xy.rows, 2, CV_32S);
	Mat tmp;
	xy(Range(0, xy.rows), Range(1, 2)).copyTo(tmp);
	minMaxLoc(tmp, &min2, &max2);
	xy(Range(0, xy.rows), Range(2, 3)).copyTo(tmp);
	minMaxLoc(tmp, &min3, &max3);
	for (int i = 0; i < nr; i++)
	{
		grid.at<int>(i, 0) = int(ceil((xy.at<double>(i, 2) - min3 ) / grid_size + 1e-6));
		grid.at<int>(i, 1) = int(ceil((xy.at<double>(i, 1) - min2) / grid_size + 1e-6));
	}
	int Max1, Max2;
	Max1 = -1; Max2 = -1;
	for (int i = 0; i < nr; i++)
	{
		Max1 = Max1 >= grid.at<int>(i, 0) ? Max1 : grid.at<int>(i, 0);
		Max2 = Max2 >= grid.at<int>(i, 1) ? Max2 : grid.at<int>(i, 1);
	}
	for (int i = 0; i < nr; i++)
	{
		if (grid.at<int>(i, 0) == Max1) grid.at<int>(i, 0) = Max1 - 1;
		if (grid.at<int>(i, 1) == Max2) grid.at<int>(i, 1) = Max2 - 1;
	}
	grid.copyTo(grid_ij);
	*rows = Max1 - 1;
	*cols = Max2 - 1;
	return 0;
}

int Utils::gen_CohDistributionOfRandomPhase(int nrand, int n_ifg, const Mat& bperp, double num_trial_wraps, Mat& Nr, int* Nr_max_nz_ix)
{
	if (nrand < 0 ||
		n_ifg < 3 ||
		bperp.rows != n_ifg ||
		bperp.cols != 1 ||
		bperp.channels() != 1 ||
		bperp.type() != CV_64F ||
		num_trial_wraps < 0.0 ||
		Nr_max_nz_ix == NULL
		)
	{
		fprintf(stderr, "gen_CohDistributionOfRandomPhase(): input check failed!\n\n");
		return -1;
	}
	nrand = nrand < 10000 ? 10000 : nrand;
	//srand(201);
	Mat rand_ifg = Mat::zeros(nrand, n_ifg, CV_64F);
//#pragma omp parallel for schedule(guided)
	//for (int i = 0; i < nrand; i++)
	//{
	//	for (int j = 0; j < n_ifg; j++)
	//	{
	//		rand_ifg.at<double>(i, j) = 2 * PI * (double(rand()) / double(RAND_MAX) - 0.5);
	//	}
	//}
	cv::randu(rand_ifg, cv::Scalar(0.0), cv::Scalar(1.0));
	rand_ifg = (rand_ifg - 0.5) * 2.0 * PI;
	Mat coh_rand = Mat::zeros(1, nrand, CV_64F);
	Mat tmp, phase_residue;
	double K0, C0, coh0, ret;
	Mat bperp0;
	bperp.copyTo(bperp0);
	for (int i = 0; i < nrand; i++)
	{
		rand_ifg(Range(i, i + 1), Range(0, n_ifg)).copyTo(tmp);
		transpose(tmp, tmp);
		ret = PS_topofit(tmp, bperp0, num_trial_wraps, &K0, &C0, &coh0, phase_residue);
		if (return_check(ret, "PS_topofit(*)", error_head)) return -1;
		coh_rand.at<double>(0, i) = coh0;
	}
	Mat pdf;
	ret = hist(coh_rand, 0.005, 0.995, 0.01, pdf);
	int x = pdf.cols - 1;
	while (fabs(pdf.at<double>(0, x)) < 1e-10)
	{
		x--;
	}
	*Nr_max_nz_ix = x;
	pdf.copyTo(Nr);
	return 0;
}

int Utils::PS_est_gamma_quick(
	const Mat& dif_phase,
	const Mat& D_A,
	const Mat& grid_ij,
	const Mat& bperp_mat,
	Mat& lowpass,
	Mat& Nr,
	double alpha,
	double beta,
	double n_trial_wraps,
	int gamma_max_iterations,
	int Nr_max_nz_ix,
	double gamma_change_convergence,
	Mat& K_ps,
	Mat& C_ps,
	Mat& coh_ps,
	int n_win
)
{
	if (dif_phase.rows != D_A.rows ||
		dif_phase.rows != bperp_mat.rows ||
		dif_phase.cols != bperp_mat.cols ||
		dif_phase.rows < 10 ||
		dif_phase.cols < 3 ||
		dif_phase.type() != CV_64F ||
		dif_phase.type() != D_A.type() ||
		dif_phase.type() != bperp_mat.type() ||
		dif_phase.channels() != 1 ||
		D_A.channels() != 1 ||
		bperp_mat.channels() != 1 ||
		grid_ij.rows != dif_phase.rows ||
		grid_ij.cols != 2 ||
		grid_ij.channels() != 1 ||
		grid_ij.type() != CV_32S ||
		Nr.rows != 1 ||
		Nr.cols != 100 ||
		Nr.channels() != 1 ||
		Nr.type() != CV_64F ||
		lowpass.rows < 7 ||
		lowpass.cols < 7 ||
		lowpass.channels() != 1 ||
		lowpass.type() != CV_64F ||
		alpha < 0.0 ||
		beta < 0.0 ||
		gamma_change_convergence < 0.00001 ||
		gamma_max_iterations < 1 ||
		n_trial_wraps < 0.0 ||
		n_win % 2 == 1 ||
		Nr_max_nz_ix < 20
		)
	{
		fprintf(stderr, "PS_est_gamma_quick(): input check failed!\n\n");
		return -1;
	}
	int n_ps = grid_ij.rows;
	int n_ifg = dif_phase.cols;
	int n_i, n_j;
	int loop_end_sw = 0;
	n_i = -1; n_j = -1;
	for (int i = 0; i < n_ps; i++)
	{
		n_i = n_i > grid_ij.at<int>(i, 0) ? n_i : grid_ij.at<int>(i, 0);
		n_j = n_j > grid_ij.at<int>(i, 1) ? n_j : grid_ij.at<int>(i, 1);
	}
	if (n_i < 1 || n_j < 1)
	{
		fprintf(stderr, "PS_est_gamma_quick(): n_i or n_j < 1!\n\n");
		return -1;
	}
	//vector<Mat> ph_grid, ph_filt;
	vector<ComplexMat> ph_grid, ph_filt;
	Mat tmp1, ph_weight_angle, ph_grid_angle;
	ComplexMat zeromat(n_i, n_j);
	ComplexMat ph_patch(n_ps, n_ifg);
	ComplexMat ph, ph_filt_angle;
	ComplexMat ph_weight, tmp2;
	Mat re, im;
	dif_phase.copyTo(tmp1);
	int ret = phase2cos(tmp1, re, im);
	if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
	ph.SetRe(re);
	ph.SetIm(im);

	Mat weighting;
	D_A.copyTo(tmp1);
	tmp1 = tmp1 + 1e-6;
	tmp1 = 1 / tmp1;
	tmp1.copyTo(weighting);
	for (int i = 0; i < n_ifg - 1; i++)
	{
		hconcat(weighting, tmp1, weighting);
	}
	
	Mat K_ps_mat = Mat::zeros(n_ps, 1, CV_64F);
	Mat C_ps_mat = Mat::zeros(n_ps, 1, CV_64F);
	Mat Coh_ps_mat = Mat::zeros(n_ps, 1, CV_64F);
	Mat ph_res = Mat::zeros(n_ps, n_ifg, CV_64F);
	for (int i = 0; i < n_ifg; i++)
	{
		ph_filt.push_back(zeromat);
	}
	int ii, jj;

	ComplexMat temp0, temp1, temp2;
	Mat bperp_tmp, phase_temp, ph_residual, Na, Prand, ph_angle;
	Mat coh_ps_save = Mat::zeros(n_ps, 1, CV_64F);
	double Kopt, Copt, cohopt, gamma_change_rms, gamma_change_change, gamma_change_save;
	double Na_sum, Nr_sum;
	gamma_change_save = 0.0;
	int low_coh_thresh = 31;
	int ret_p;
	int i_loop = 1;
	int idx;

	while (loop_end_sw == 0)
	{
		ph_grid.clear();
		for (int i = 0; i < n_ifg; i++)
		{
			ph_grid.push_back(zeromat);
		}
		K_ps_mat.copyTo(tmp1);
		for (int i = 0; i < n_ifg - 1; i++)
		{
			hconcat(tmp1, K_ps_mat, tmp1);
		}
		tmp1 = tmp1.mul(bperp_mat);
		tmp1 = -tmp1;
		ret = phase2cos(tmp1, re, im);
		if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
		tmp2.SetRe(re);
		tmp2.SetIm(im);
		ph_weight = ph * tmp2;
		ph_weight = ph_weight * weighting;
//#pragma omp parallel for schedule(guided)
		for (int j = 0; j < n_ifg; j++)
		{
			int ii, jj;
			ComplexMat tmp2;
			for (int i = 0; i < n_ps; i++)
			{
				ii = grid_ij.at<int>(i, 0) - 1;
				jj = grid_ij.at<int>(i, 1) - 1;
				if (ii >= 0 && jj >= 0)
				{
					tmp2 = ph_grid[j](Range(ii, ii + 1), Range(jj, jj + 1)) + ph_weight(Range(i, i + 1), Range(j, j + 1));
					ph_grid[j].SetValue(Range(ii, ii + 1), Range(jj, jj + 1), tmp2);
				}
			}
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_ifg; i++)
		{
			ComplexMat ph_filt_angle;
			int ret;
			ret = Clap(ph_grid[i], ph_filt_angle, alpha, beta, int(0.75 * n_win), (n_win - int(0.75 * n_win)), lowpass);
			ph_filt[i].SetValue(Range(0, n_i), Range(0, n_j), ph_filt_angle);
		}

#pragma omp parallel for schedule(guided)
		for (int j = 0; j < n_ifg; j++)
		{
			int ii, jj;
			ComplexMat tmp2;
			for (int i = 0; i < n_ps; i++)
			{
				ii = grid_ij.at<int>(i, 0) - 1;
				jj = grid_ij.at<int>(i, 1) - 1;
				if (ii > -1 && jj > -1)
				{
					tmp2 = ph_filt[j](Range(ii, ii + 1), Range(jj, jj + 1));
					ph_patch.SetValue(Range(i, i + 1), Range(j, j + 1), tmp2);
				}

			}
		}

		ph_angle = ph_patch.GetPhase();
		ret = phase2cos(ph_angle, re, im);
		if (return_check(ret, "phase2cos(*, *, *)", error_head)) return -1;
		ph_patch.SetRe(re);
		ph_patch.SetIm(im);


		////////////////////////Now estimate topo error//////////////////////////
		ph_patch = ph_patch.conj();

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < n_ps; i++)
		{
			ComplexMat temp0, temp1;
			Mat bperp_tmp, phase_temp, ph_residual;
			double Kopt, Copt, cohopt;
			Kopt = 0.0; Copt = 0.0; cohopt = 0.0;
			int ret_p;
			temp0 = ph(Range(i, i + 1), Range(0, n_ifg));
			temp1 = ph_patch(Range(i, i + 1), Range(0, n_ifg));
			temp0 = temp0 * temp1;
			phase_temp = temp0.GetPhase();
			transpose(phase_temp, phase_temp);
			bperp_mat(Range(i, i + 1), Range(0, n_ifg)).copyTo(bperp_tmp);
			transpose(bperp_tmp, bperp_tmp);
			if (temp0.countNonzero() > 0)
			{
				ret_p = PS_topofit(phase_temp, bperp_tmp, n_trial_wraps, &Kopt, &Copt, &cohopt, ph_residual);
			}
			K_ps_mat.at<double>(i, 0) = Kopt;
			C_ps_mat.at<double>(i, 0) = Copt;
			Coh_ps_mat.at<double>(i, 0) = cohopt;
		}
		coh_ps_save = Coh_ps_mat - coh_ps_save;
		coh_ps_save = coh_ps_save.mul(coh_ps_save);
		gamma_change_rms = sum(coh_ps_save)[0];
		gamma_change_rms = gamma_change_rms / double(n_ps);
		gamma_change_rms = sqrt(gamma_change_rms);
		gamma_change_change = gamma_change_rms - gamma_change_save;
		gamma_change_save = gamma_change_rms;
		Coh_ps_mat.copyTo(coh_ps_save);

		if (abs(gamma_change_change) < gamma_change_convergence || i_loop >= gamma_max_iterations)
		{
			loop_end_sw = 1;
		}
		else
		{
			i_loop++;
			transpose(Coh_ps_mat, Coh_ps_mat);
			ret = hist(Coh_ps_mat, 0.005, 0.995, 0.01, Na);
			if (return_check(ret, "hist()", error_head)) return -1;
			Na_sum = sum(Na(Range(0, 1), Range(0, low_coh_thresh - 1)))[0];
			Nr_sum = sum(Nr(Range(0, 1), Range(0, low_coh_thresh - 1)))[0];
			Nr = Nr * (Na_sum / (Nr_sum + 1e-6));
			for (int i = 0; i < Na.cols; i++)
			{
				if (Na.at<double>(0, i) < 0.5) Na.at<double>(0, i) = 1.0;
			}
			Prand = Nr / Na;
			for (int i = 0; i < low_coh_thresh; i++)
			{
				Prand.at<double>(0, i) = 1.0;
			}
			for (int i = Nr_max_nz_ix; i < Prand.cols; i++)
			{
				Prand.at<double>(0, i) = 0.0;
			}
			for (int i = 0; i < Prand.cols; i++)
			{
				Prand.at<double>(0, i) = Prand.at<double>(0, i) > 1.0 ? 1.0 : Prand.at<double>(0, i);
			}
			GaussianBlur(Prand, Prand, cv::Size(7, 7), 1);
			resize(Prand, Prand, Size(0, 0), 10.0, 1.0, cv::INTER_LINEAR);
#pragma omp parallel for schedule(guided)
			for (int i = 0; i < n_ps; i++)
			{
				int idx;
				for (int j = 0; j < n_ifg; j++)
				{
					idx = (int)(floor(Coh_ps_mat.at<double>(0, i) * 1000) + 1);
					idx = idx < 0 ? 0 : idx;
					idx = idx > Prand.cols - 1 ? Prand.cols - 1 : idx;
					weighting.at<double>(i, j) = (1 - Prand.at<double>(0, idx)) * (1 - Prand.at<double>(0, idx));
				}
			}
			transpose(Coh_ps_mat, Coh_ps_mat);
		}
	}
	K_ps_mat.copyTo(K_ps);
	C_ps_mat.copyTo(C_ps);
	Coh_ps_mat.copyTo(coh_ps);
	return 0;
}

int Utils::coh_thresh_fit(const Mat& D_A, Mat& coh_ps, Mat& Nr, Mat& coh_thresh_D_A, int low_coh_thresh, double false_alarm)
{
	if (D_A.rows < 3 ||
		D_A.cols != 1 ||
		D_A.type() != CV_64F ||
		D_A.channels() != 1 ||
		coh_ps.rows != D_A.rows ||
		coh_ps.cols != D_A.cols ||
		coh_ps.channels() != 1 ||
		coh_ps.type() != CV_64F ||
		Nr.rows != 1 ||
		Nr.cols != 100 ||
		Nr.channels() != 1 ||
		Nr.type() != CV_64F ||
		low_coh_thresh < 1 ||
		low_coh_thresh > 100 ||
		false_alarm < 0.001 ||
		false_alarm > 0.2
		)
	{
		fprintf(stderr, "coh_thresh_fit(): input check failed!\n\n");
		return -1;
	}
	Mat D_A_sort, D_A_max, min_coh, D_A_mean, Nr_dist;
	int bin_size;
	int n_ps = D_A.rows;
	if (n_ps < 10000)
	{
		fprintf(stderr, "coh_thresh_fit(): not enough(10000) PS candidates!\n\n");
		return -1;
	}
	cv::sort(D_A, D_A_sort, SORT_ASCENDING + SORT_EVERY_COLUMN);
	if (n_ps < 50000)
	{
		bin_size = 2000;
	}
	else
	{
		bin_size = 10000;
	}
	int D_A_max_size = int(ceil(n_ps / bin_size)) + 1;
	D_A_max = Mat::zeros(D_A_max_size, 1, CV_64F);
	Mat D_A_max_idx = Mat::zeros(D_A_max_size, 1, CV_32S);
	for (int i = 1; i < D_A_max_size - 1; i++)
	{
		D_A_max.at<double>(i, 0) = D_A_sort.at<double>(i * bin_size, 0);
		D_A_max_idx.at<int>(i, 0) = i * bin_size;
	}
	D_A_max.at<double>(0, 0) = 0.0;
	D_A_max_idx.at<int>(0, 0) = 0;
	D_A_max.at<double>(D_A_max_size - 1, 0) = D_A_sort.at<double>(D_A_sort.rows - 1, 0);
	D_A_max_idx.at<int>(D_A_max_size - 1, 0) = D_A_sort.rows;

	min_coh = Mat::zeros(D_A_max_size - 1, 1, CV_64F);
	D_A_mean = Mat::zeros(D_A_max_size - 1, 1, CV_64F);
	Nr.copyTo(Nr_dist);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < D_A_max_size - 1; i++)
	{
		Mat coh_chunk = Mat::zeros(1, bin_size, CV_64F);
		Mat D_A_tmp, Na, Nr, percent_rand, a, A, A_t, x, b;
		double low, up;
		int count = 0;
		int min_fit_ix, max_fit_ix;
		bool flag = false;
		D_A_sort(Range(D_A_max_idx.at<int>(i, 0), D_A_max_idx.at<int>(i + 1, 0)), Range(0, 1)).copyTo(D_A_tmp);
		D_A_mean.at<double>(i, 0) = mean(D_A_tmp)[0];
		low = D_A_max.at<double>(i, 0);
		up = D_A_max.at<double>(i + 1, 0);
		for (int j = 0; j < n_ps; j++)
		{
			if (D_A.at<double>(j, 0) > low && D_A.at<double>(j, 0) <= up && count < bin_size)
			{
				coh_chunk.at<double>(0, count) = coh_ps.at<double>(j, 0);
				count++;
			}
		}
		hist(coh_chunk, 0.005, 0.995, 0.01, Na);
		//cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\Na.bin", Na);
		Nr = Nr_dist * (sum(Na(Range(0, 1), Range(0, low_coh_thresh)))[0] / (sum(Nr_dist(Range(0, 1), Range(0, low_coh_thresh)))[0] + 1e-6));
		for (int j = 0; j < Na.cols; j++)//avoid divide by zero
		{
			if (fabs(Na.at<double>(0, j)) < 1e-9) Na.at<double>(0, j) = 1.0;
		}
		flip(Nr, Nr, 1);
		cumsum(Nr, 2);
		flip(Na, Na, 1);
		cumsum(Na, 2);
		percent_rand = Nr / Na * 100.0;
		flip(percent_rand, percent_rand, 1);
		//cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\percent_rand.bin", percent_rand);
		for (int k = 0; k < percent_rand.cols; k++)
		{
			if (percent_rand.at<double>(0, k) <= false_alarm && !flag)
			{
				min_fit_ix = k - 2;
				flag = true;
			}
		}
		if (min_fit_ix > 1 && min_fit_ix <= 95)
		{
			max_fit_ix = min_fit_ix + 5;
			percent_rand(Range(0, 1), Range(min_fit_ix, max_fit_ix)).copyTo(a);
			double std, mean;
			mean = cv::mean(a)[0];
			this->std(a, &std);
			a = (a - mean) / (std + 1e-12);
			b = Mat::zeros(5, 1, CV_64F);
			for (int k = min_fit_ix; k < max_fit_ix; k++)
			{
				b.at<double>(k - min_fit_ix, 0) = double(0.01 * k);
			}
			A = Mat::ones(5, 4, CV_64F);
			for (int k = 0; k < 5; k++)
			{
				for (int kk = 1; kk < 4; kk++)
				{
					A.at<double>(k, kk) = std::pow(a.at<double>(0, k), kk);
				}
			}
			transpose(A, A_t);
			b = A_t * b;
			A = A_t * A;
			double false_alarm_tmp;
			if (solve(A, b, x, cv::DECOMP_LU))
			{

				false_alarm_tmp = (false_alarm - mean) / (std + 1e-12);
				min_coh.at<double>(i, 0) = x.at<double>(0, 0) + x.at<double>(1, 0) * false_alarm_tmp +
					x.at<double>(2, 0) * pow(false_alarm_tmp, 2) + x.at<double>(3, 0) * pow(false_alarm_tmp, 3) + mean;
			}
			else
			{
				min_coh.at<double>(i, 0) = 0.0;
			}
		}
		else
		{
			min_coh.at<double>(i, 0) = 0.0;
		}
		

	}
	int count = 0;
	for (int i = 0; i < min_coh.rows; i++)
	{
		if (min_coh.at<double>(i, 0) > 1e-9) count++;
	}
	//cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\min_coh.bin", min_coh);
	Mat temp_ones = Mat::ones(n_ps, 1, CV_64F);
	if (count >= 2)
	{
		Mat D_A_mean_new = Mat::zeros(count, 1, CV_64F);
		Mat min_coh_new = Mat::zeros(count, 1, CV_64F);
		int count1 = 0;
		for (int i = 0; i < min_coh.rows; i++)
		{
			if (min_coh.at<double>(i, 0) > 1e-9 && count1 < count)
			{
				D_A_mean_new.at<double>(count1, 0) = D_A_mean.at<double>(i, 0);
				min_coh_new.at<double>(count1, 0) = min_coh.at<double>(i, 0);
				count1++;
			}
		}
		Mat A = Mat::ones(count, 2, CV_64F);
		Mat b = Mat::zeros(count, 1, CV_64F);
		Mat A_t;
		for (int i = 0; i < count; i++)
		{
			A.at<double>(i, 1) = D_A_mean_new.at<double>(i, 0);
			b.at<double>(i, 0) = min_coh_new.at<double>(i, 0);
		}
		transpose(A, A_t);
		b = A_t * b;
		A = A_t * A;
		Mat x;
		
		if (solve(A, b, x, DECOMP_LU))
		{
			double a0 = x.at<double>(0, 0);
			double a1 = x.at<double>(1, 0);
			if (a1 > 0.0)//拟合斜率为正
			{
				hconcat(temp_ones, D_A, temp_ones);
				coh_thresh_D_A = temp_ones * x;
			}
			else
			{
				coh_thresh_D_A = temp_ones - 0.65;
			}
		}
		else
		{
			coh_thresh_D_A = temp_ones - 0.65;
		}
		
	}
	else
	{
		coh_thresh_D_A = temp_ones - 0.65;
	}
	return 0;
}

int Utils::PS_topovel_fit_search(ComplexMat& ph, Mat& coef_delta_vel, Mat& coef_delta_height, double radius_delta_vel, double radius_delta_height, double* MC, double* delta_vel, double* delta_height)
{
	if (ph.GetCols() != 1 ||
		ph.GetRows() < 3 ||
		ph.type() != CV_64F ||
		coef_delta_vel.rows != ph.GetRows() ||
		coef_delta_height.rows != ph.GetRows() ||
		coef_delta_vel.cols != 1 ||
		coef_delta_height.cols != 1 ||
		coef_delta_height.channels() != 1 ||
		coef_delta_vel.channels() != 1 ||
		coef_delta_height.type() != CV_64F ||
		coef_delta_vel.type() != CV_64F ||
		radius_delta_vel < 0.0 ||
		radius_delta_height < 0.0 ||
		MC == NULL ||
		delta_height == NULL ||
		delta_vel == NULL
		)
	{
		fprintf(stderr, "PS_topovel_fit_search(): input check failed!\n\n");
		return -1;
	}
	int M = ph.GetRows();
	//搜索区间间隔(1mm/y, 2m),可以适当调整
	double interval_v = 0.001;
	double interval_height = 2.0;
	int nr = int(floor(2 * radius_delta_height / interval_height)) + 1;
	int nc = int(floor(2 * radius_delta_vel / interval_v)) + 1;
	if (nr < 2 || nc < 2)
	{
		fprintf(stderr, "PS_topovel_fit_search(): search radius too small!\n\n");
		return -1;
	}

	Mat MC_model(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double delta_vel_tmp, delta_height_tmp;
		ComplexMat temp1;
		Mat tmp1, tmp2;
		for (int j = 0; j < nc; j++)
		{
			delta_vel_tmp = double(j) * interval_v - radius_delta_vel;
			delta_height_tmp = double(i) * interval_height - radius_delta_height;
			tmp1 = coef_delta_vel * delta_vel_tmp;
			tmp1 = tmp1 + coef_delta_height * delta_height_tmp;
			phase2cos(tmp1, tmp1, tmp2);
			temp1.SetRe(tmp1);
			temp1.SetIm(tmp2);
			temp1 = ph * (temp1.conj());
			temp1 = temp1.sum();
			MC_model.at<double>(i, j) = temp1.GetMod().at<double>(0, 0) / double(M);
		}
	}
	cv::Point p;
	double max_val;
	//cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\MC_model1.bin", MC_model);
	minMaxLoc(MC_model, NULL, &max_val, NULL, &p);
	if (max_val <= 0.0)
	{
		fprintf(stderr, "PS_topovel_fit_search(): failed to locate the maximum value!\n\n");
		return -1;
	}
	*delta_vel = double(p.x) * interval_v - radius_delta_vel;
	*delta_height = double(p.y) * interval_height - radius_delta_height;
	*MC = max_val;
	return 0;
}

int Utils::PS_topovel_fit_search_2(ComplexMat& ph, Mat& coef_delta_vel, Mat& coef_delta_height, double rough_delta_vel_center, double rough_delta_height_center, double radius_delta_vel, double radius_delta_height, double* MC, double* delta_vel, double* delta_height)
{
	if (ph.GetCols() != 1 ||
		ph.GetRows() < 3 ||
		ph.type() != CV_64F ||
		coef_delta_vel.rows != ph.GetRows() ||
		coef_delta_height.rows != ph.GetRows() ||
		coef_delta_vel.cols != 1 ||
		coef_delta_height.cols != 1 ||
		coef_delta_height.channels() != 1 ||
		coef_delta_vel.channels() != 1 ||
		coef_delta_height.type() != CV_64F ||
		coef_delta_vel.type() != CV_64F ||
		radius_delta_vel < 0.0 ||
		radius_delta_height < 0.0 ||
		MC == NULL ||
		delta_height == NULL ||
		delta_vel == NULL
		)
	{
		fprintf(stderr, "PS_topovel_fit_search(): input check failed!\n\n");
		return -1;
	}
	int M = ph.GetRows();
	//搜索区间间隔(0.1mm/y, 0.1m)
	double interval_v = 0.0001;
	double interval_height = 0.2;
	int nr = int(floor(2 * radius_delta_height / interval_height)) + 1;
	int nc = int(floor(2 * radius_delta_vel / interval_v)) + 1;
	if (nr < 2 || nc < 2)
	{
		fprintf(stderr, "PS_topovel_fit_search(): search radius too small!\n\n");
		return -1;
	}

	Mat MC_model(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		double delta_vel_tmp, delta_height_tmp;
		ComplexMat temp1;
		Mat tmp1, tmp2;
		for (int j = 0; j < nc; j++)
		{
			delta_vel_tmp = double(j) * interval_v - radius_delta_vel + rough_delta_vel_center;
			delta_height_tmp = double(i) * interval_height - radius_delta_height + rough_delta_height_center;
			tmp1 = coef_delta_vel * delta_vel_tmp;
			tmp1 = tmp1 + coef_delta_height * delta_height_tmp;
			phase2cos(tmp1, tmp1, tmp2);
			temp1.SetRe(tmp1);
			temp1.SetIm(tmp2);
			temp1 = ph * (temp1.conj());
			temp1 = temp1.sum();
			MC_model.at<double>(i, j) = temp1.GetMod().at<double>(0, 0) / double(M);
		}
	}
	cv::Point p;
	double max_val;
	//cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\MC_model.bin", MC_model);
	minMaxLoc(MC_model, NULL, &max_val, NULL, &p);
	if (max_val <= 0.0)
	{
		fprintf(stderr, "PS_topovel_fit_search(): failed to locate the maximum value!\n\n");
		return -1;
	}
	*delta_vel = double(p.x) * interval_v - radius_delta_vel + rough_delta_vel_center;
	*delta_height = double(p.y) * interval_height - radius_delta_height + rough_delta_height_center;
	*MC = max_val;
	return 0;
}

int Utils::Dump_unqualified_pairnodes(vector<tri_node>& nodes, tri_edge* edges, int num_edges, double MC_thresh)
{
	if (nodes.size() < 3 ||
		edges == NULL ||
		num_edges < 3
		)
	{
		fprintf(stderr, "Dump_unqualified_pairnodes(): input check failed!\n\n");
		return -1;
	}
	if (MC_thresh > 1.0)
	{
		MC_thresh = 0.65;
	}
	int num_nodes = nodes.size();
	long* neigh_ptr = NULL;
	int num_neigh;
	for (int i = 0; i < num_nodes; i++)
	{
		int count = 0;
		nodes[i].get_neigh_ptr(&neigh_ptr, &num_neigh);
		for (int j = 0; j < num_neigh; j++)
		{
			if ((edges + j)->MC < MC_thresh) count++;
		}
		if ((double(count) / (double(num_neigh) + 1e-10)) > 0.5)
		{
			nodes[i].set_balance(false);
		}
	}
	//for (int i = 0; i < num_edges; i++)
	//{
	//	if ((edges + i)->MC < MC_thresh)
	//	{
	//		nodes[(edges + i)->end1 - 1].set_balance(false);
	//		nodes[(edges + i)->end2 - 1].set_balance(false);
	//	}
	//}
	return 0;
}

int Utils::PS_topovel_incre(vector<tri_node>& nodes, tri_edge* edges, int num_edges, double distance_thresh, double MC_thresh)
{
	if (nodes.size() < 3 ||
		edges == NULL ||
		num_edges < 3
		)
	{
		fprintf(stderr, "PS_topovel_incre(): input check failed!\n\n");
		return -1;
	}
	if (distance_thresh < 1.0) distance_thresh = 1.0;
	if (MC_thresh > 0.9) MC_thresh = 0.9;
	//找到增量积分起始点
	int ix = 0;
	double MC = -1.0;
	for (int i = 0; i < num_edges; i++)
	{
		if ((edges + i)->MC > MC)
		{
			ix = i;
			MC = (edges + i)->MC;
		}
	}
	int start = (edges + ix)->end1;

	//采用类似质量图法解缠的算法进行增量积分集成
	edge_index tmp;
	priority_queue<edge_index> que;
	nodes[start - 1].set_vel(0.0);//起始点形变速率和高程误差设置为0，后续可根据参考点进行校正
	nodes[start - 1].set_height(0.0);
	nodes[start - 1].set_status(true);
	long* ptr_neigh = NULL;
	int num_neigh, end2, number, row1, col1, row2, col2;
	double distance, vel, height, delta_v, delta_h, MC_total, vel_total, height_total;
	
	nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	for (int i = 0; i < num_neigh; i++)
	{
		end2 = (edges + *(ptr_neigh + i) - 1)->end1 == start ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= distance_thresh &&
			(edges + *(ptr_neigh + i) - 1)->MC > MC_thresh
			)
		{
			tmp.num = *(ptr_neigh + i);
			tmp.quality = 1 - (edges + *(ptr_neigh + i) - 1)->MC;
			que.push(tmp);
		}
	}

	while (que.size() != 0)
	{
		tmp = que.top();
		que.pop();
		if (nodes[(edges + tmp.num - 1)->end1 - 1].get_status())
		{
			number = (edges + tmp.num - 1)->end1;
			end2 = (edges + tmp.num - 1)->end2;
		}
		else
		{
			number = (edges + tmp.num - 1)->end2;
			end2 = (edges + tmp.num - 1)->end1;
		}
		MC_total = 1e-10;
		vel_total = 0.0;
		height_total = 0.0;
		if (!nodes[end2 - 1].get_status())
		{
			nodes[end2 - 1].get_pos(&row2, &col2);
			nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			for (int i = 0; i < num_neigh; i++)
			{
				number = (edges + *(ptr_neigh + i) - 1)->end1 == end2 ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
				if (nodes[number - 1].get_status())
				{
					vel = nodes[number - 1].get_vel();
					height = nodes[number - 1].get_height();
					nodes[number - 1].get_pos(&row1, &col1);
					if (row1 > row2)
					{
						delta_h = -(edges + *(ptr_neigh + i) - 1)->delta_height;
						delta_v = -(edges + *(ptr_neigh + i) - 1)->delta_vel;
					}
					if (row1 < row2)
					{
						delta_h = (edges + *(ptr_neigh + i) - 1)->delta_height;
						delta_v = (edges + *(ptr_neigh + i) - 1)->delta_vel;
					}
					if (row1 == row2)
					{
						if (col1 > col2)
						{
							delta_h = -(edges + *(ptr_neigh + i) - 1)->delta_height;
							delta_v = -(edges + *(ptr_neigh + i) - 1)->delta_vel;
						}
						else
						{
							delta_h = (edges + *(ptr_neigh + i) - 1)->delta_height;
							delta_v = (edges + *(ptr_neigh + i) - 1)->delta_vel;
						}
					}
					vel_total += (delta_v + vel) * (edges + *(ptr_neigh + i) - 1)->MC;
					height_total += (delta_h + height) * (edges + *(ptr_neigh + i) - 1)->MC;
					MC_total += (edges + *(ptr_neigh + i) - 1)->MC;
				}
			}
			nodes[end2 - 1].set_height(height_total / MC_total);
			nodes[end2 - 1].set_vel(vel_total / MC_total);
			nodes[end2 - 1].set_status(true);
			nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			number = end2;
			for (int i = 0; i < num_neigh; i++)
			{

				end2 = (edges + *(ptr_neigh + i) - 1)->end1 == number ? (edges + *(ptr_neigh + i) - 1)->end2 : (edges + *(ptr_neigh + i) - 1)->end1;
				nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
				if (!nodes[end2 - 1].get_status() &&
					distance <= distance_thresh &&
					(edges + *(ptr_neigh + i) - 1)->MC > MC_thresh
					)
				{
					tmp.num = *(ptr_neigh + i);
					tmp.quality = 1 - (edges + *(ptr_neigh + i) - 1)->MC;
					que.push(tmp);
				}
			}
		}
	}
	return 0;
}

int Utils::stateVec_interp(Mat& stateVec, double time_interval, Mat& stateVec_interp)
{
	if (stateVec.empty() ||
		stateVec.cols != 7||
		stateVec.rows < 7||
		time_interval < 0.0)
	{
		fprintf(stderr, "stateVec_interp(): input check failed!\n");
		return -1;
	}
	Mat statevec;
	if (stateVec.type() != CV_64F) stateVec.convertTo(statevec, CV_64F);
	else stateVec.copyTo(statevec);

	int rows = statevec.rows; int cols = statevec.cols;
	Mat time; statevec(cv::Range(0, rows), cv::Range(0, 1)).copyTo(time);
	time = time - time.at<double>(0, 0);
	Mat A = Mat::ones(rows, 6, CV_64F);
	Mat temp, b;
	//拟合x
	Mat x; statevec(cv::Range(0, rows), cv::Range(1, 2)).copyTo(x);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * x;
	A = temp * A;
	if (!cv::solve(A, b, x, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//拟合y
	A = Mat::ones(rows, 6, CV_64F);
	Mat y; statevec(cv::Range(0, rows), cv::Range(2, 3)).copyTo(y);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * y;
	A = temp * A;
	if (!cv::solve(A, b, y, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//拟合z

	A = Mat::ones(rows, 6, CV_64F);
	Mat z; statevec(cv::Range(0, rows), cv::Range(3, 4)).copyTo(z);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * z;
	A = temp * A;
	if (!cv::solve(A, b, z, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//拟合vx

	A = Mat::ones(rows, 6, CV_64F);
	Mat vx; statevec(cv::Range(0, rows), cv::Range(4, 5)).copyTo(vx);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * vx;
	A = temp * A;
	if (!cv::solve(A, b, vx, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//拟合vy

	A = Mat::ones(rows, 6, CV_64F);
	Mat vy; statevec(cv::Range(0, rows), cv::Range(5, 6)).copyTo(vy);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * vy;
	A = temp * A;
	if (!cv::solve(A, b, vy, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//拟合vz

	A = Mat::ones(rows, 6, CV_64F);
	Mat vz; statevec(cv::Range(0, rows), cv::Range(6, 7)).copyTo(vz);
	time.copyTo(A(cv::Range(0, rows), cv::Range(1, 2)));
	temp = time.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(2, 3)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(3, 4)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(4, 5)));
	temp = temp.mul(time);
	temp.copyTo(A(cv::Range(0, rows), cv::Range(5, 6)));
	transpose(A, temp);
	b = temp * vz;
	A = temp * A;
	if (!cv::solve(A, b, vz, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "stateVec_interp(): matrix defficiency!\n");
		return -1;
	}

	//插值

	int count = 1;
	double t = 0;
	while (t <= time.at<double>(rows - 1, 0))
	{
		count++;
		t += time_interval;
	}
	stateVec_interp.create(count, 7, CV_64F);
	Mat tt = Mat::ones(count, 6, CV_64F);
	tt(cv::Range(0, count), cv::Range(0, 1)).copyTo(stateVec_interp(cv::Range(0, count), cv::Range(0, 1)));
	t = 0.0;
	for (int i = 0; i < count; i++)
	{
		tt.at<double>(i, 1) = t;
		tt.at<double>(i, 2) = t * t;
		tt.at<double>(i, 3) = t * t * t;
		tt.at<double>(i, 4) = t * t * t * t;
		tt.at<double>(i, 5) = t * t * t * t * t;
		t += time_interval;
	}
	x = tt * x;
	y = tt * y;
	z = tt * z;
	vx = tt * vx;
	vy = tt * vy;
	vz = tt * vz;
	x.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(1, 2)));
	y.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(2, 3)));
	z.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(3, 4)));
	vx.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(4, 5)));
	vy.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(5, 6)));
	vz.copyTo(stateVec_interp(cv::Range(0, count), cv::Range(6, 7)));
	return 0;
}

int Utils::get_AOI_from_h5SLC(const char* h5_file, double lon_topleft, double lat_topleft, double lon_bottomright, double lat_bottomright, ComplexMat& slc, int* offset_row, int* offset_col)
{
	if (h5_file == NULL ||
		fabs(lon_topleft) > 180.0 ||
		fabs(lon_bottomright) > 180.0 ||
		fabs(lat_topleft) > 90.0 ||
		fabs(lat_bottomright) > 90.0)
	{
		fprintf(stderr, "get_AOI_from_h5SLC(): input check failed!\n");
		return -1;
	}
	int ret, row_start, row_end, col_start, col_end;double lon_start, lat_start, lon_end, lat_end;
	FormatConversion conversion;
	Mat row_coef, col_coef, lon_coef, lat_coef;

	//获取经纬度范围
	int rows, cols;
	Mat tmp;
	ret = conversion.read_array_from_h5(h5_file, "azimuth_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	rows = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "range_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	cols = tmp.at<int>(0, 0);

	ret = conversion.read_array_from_h5(h5_file, "lon_coefficient", lon_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5_file, "lat_coefficient", lat_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	Mat tmp_row = Mat::zeros(1, 1, CV_64F); tmp_row.at<double>(0, 0) = 1;
	Mat tmp_col = Mat::zeros(1, 1, CV_64F); tmp_col.at<double>(0, 0) = 1;
	Mat tmp_out;
	ret = coord_conversion(lon_coef, tmp_row, tmp_col, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	lon_start = (tmp_out.at<double>(0, 0));
	ret = coord_conversion(lat_coef, tmp_row, tmp_col, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	lat_start = (tmp_out.at<double>(0, 0));
	tmp_row.at<double>(0, 0) = rows; tmp_col.at<double>(0, 0) = cols;
	ret = coord_conversion(lon_coef, tmp_row, tmp_col, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	lon_end = (tmp_out.at<double>(0, 0));
	ret = coord_conversion(lat_coef, tmp_row, tmp_col, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	lat_end = (tmp_out.at<double>(0, 0));

	double x;
	if (lon_start > lon_end)
	{
		x = lon_start;
		lon_start = lon_end;
		lon_end = x;
	}
	if (lat_start > lat_end)
	{
		x = lat_start;
		lat_start = lat_end;
		lat_end = x;
	}
	if (lon_topleft > lon_bottomright)
	{
		x = lon_topleft;
		lon_topleft = lon_bottomright;
		lon_bottomright = x;
	}
	if (lat_topleft < lat_bottomright)
	{
		x = lat_topleft;
		lat_topleft = lat_bottomright;
		lat_bottomright = x;
	}

	lon_topleft = lon_topleft < lon_start ? lon_start : lon_topleft;
	lon_topleft = lon_topleft < lon_end ? lon_topleft : lon_start;

	lat_topleft = lat_topleft < lat_end ? lat_topleft : lat_end;
	lat_topleft = lat_topleft < lat_start ? lat_end : lat_topleft;

	lon_bottomright = lon_bottomright < lon_end ? lon_bottomright : lon_end;
	lon_bottomright = lon_bottomright < lon_start ? lon_end : lon_bottomright;

	lat_bottomright = lat_bottomright < lat_start ? lat_start : lat_bottomright;
	lat_bottomright = lat_bottomright < lat_end ? lat_bottomright : lat_start;

	ret = conversion.read_array_from_h5(h5_file, "row_coefficient", row_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5_file, "col_coefficient", col_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	Mat tmp_lon = Mat::zeros(1, 1, CV_64F); tmp_lon.at<double>(0, 0) = lon_topleft;
	Mat tmp_lat = Mat::zeros(1, 1, CV_64F); tmp_lat.at<double>(0, 0) = lat_topleft;
	
	ret = coord_conversion(row_coef, tmp_lon, tmp_lat, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	row_start = (int)floor(tmp_out.at<double>(0, 0));
	ret = coord_conversion(col_coef, tmp_lon, tmp_lat, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	col_start = (int)floor(tmp_out.at<double>(0, 0));

	tmp_lon.at<double>(0, 0) = lon_bottomright; tmp_lat.at<double>(0, 0) = lat_bottomright;
	ret = coord_conversion(row_coef, tmp_lon, tmp_lat, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	row_end = (int)floor(tmp_out.at<double>(0, 0));
	ret = coord_conversion(col_coef, tmp_lon, tmp_lat, tmp_out);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	col_end = (int)floor(tmp_out.at<double>(0, 0));

	//读取slc数据
	
	row_start = row_start < 0 ? 0 : row_start;
	row_start = row_start > (rows - 1) ? (rows - 1) : row_start;
	row_end = row_end < 0 ? 0 : row_end;
	row_end = row_end > (rows - 1) ? (rows - 1) : row_end;
	col_start = col_start < 0 ? 0 : col_start;
	col_start = col_start > (cols - 1) ? (cols - 1) : col_start;
	col_end = col_end < 0 ? 0 : col_end;
	col_end = col_end > (cols - 1) ? (cols - 1) : col_end;

	int start_r, end_r, start_c, end_c;
	start_r = row_start < row_end ? row_start : row_end;
	end_r = row_start < row_end ? row_end : row_start;
	start_c = col_start < col_end ? col_start : col_end;
	end_c = col_start < col_end ? col_end : col_start;
	if(offset_row) *offset_row = start_r;
	if (offset_col)*offset_col = start_c;
	ret = conversion.read_subarray_from_h5(h5_file, "s_re", start_r, start_c, (end_r - start_r), (end_c - start_c), slc.re);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_subarray_from_h5(h5_file, "s_im", start_r, start_c, (end_r - start_r), (end_c - start_c), slc.im);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	return 0;
}

int Utils::get_AOI_from_h5slc(const char* h5_file, double lon_center, double lat_center, double width, double height, ComplexMat& slc, int* offset_row, int* offset_col)
{
	if (h5_file == NULL ||
		fabs(lon_center) > 180.0 ||
		fabs(lat_center) > 90.0 ||
		width < 0.0 ||
		height < 0.0)
	{
		fprintf(stderr, "get_AOI_from_h5slc(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion;
	int ret, row_center, col_center, total_rows, total_cols, AOI_rows, AOI_cols, offset_r, offset_c;
	double inc_center, range_spacing, azimuth_spacing;
	Mat row_coef, col_coef, lon, lat, tmp;
	lon = Mat::zeros(1, 1, CV_64F); lon.at<double>(0, 0) = lon_center;
	lat = Mat::zeros(1, 1, CV_64F); lat.at<double>(0, 0) = lat_center;
	//读取图像行列总数
	ret = conversion.read_array_from_h5(h5_file, "azimuth_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	total_rows = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "range_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	total_cols = tmp.at<int>(0, 0);
	//读取采样间隔和下视角
	ret = conversion.read_array_from_h5(h5_file, "azimuth_spacing", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	azimuth_spacing = tmp.at<double>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "range_spacing", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	range_spacing = tmp.at<double>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "incidence_center", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	inc_center = tmp.at<double>(0, 0);
	//确定AOI中心图像坐标
	ret = conversion.read_array_from_h5(h5_file, "row_coefficient", row_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5_file, "col_coefficient", col_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = coord_conversion(row_coef, lon, lat, tmp);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	row_center = (int)floor(tmp.at<double>(0, 0));
	ret = coord_conversion(col_coef, lon, lat, tmp);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	col_center = (int)floor(tmp.at<double>(0, 0));
	if (row_center < 1 || row_center > total_rows - 1 || col_center < 1 || col_center > total_cols - 1)
	{
		fprintf(stderr, "get_AOI_from_h5slc(): AOI not found in %s!\n", h5_file);
		return -1;
	}

	AOI_rows = (int)floor(height / azimuth_spacing);
	AOI_cols = (int)floor(width / range_spacing * sin(inc_center / 180.0 * PI));

	AOI_rows = AOI_rows > total_rows ? total_rows : AOI_rows;
	AOI_cols = AOI_cols > total_cols ? total_cols : AOI_cols;

	offset_r = row_center - (int)AOI_rows / 2;
	offset_r = offset_r < 0 ? 0 : offset_r;
	offset_c = col_center - (int)AOI_cols / 2;
	offset_c = offset_c < 0 ? 0 : offset_c;
	AOI_rows = (total_rows - offset_r) > AOI_rows ? AOI_rows : (total_rows - offset_r);
	AOI_cols = (total_cols - offset_c) > AOI_cols ? AOI_cols : (total_cols - offset_c);
	if (offset_row) *offset_row = offset_r;
	if (offset_col)*offset_col = offset_c;
	ret = conversion.read_subarray_from_h5(h5_file, "s_re", offset_r, offset_c, AOI_rows, AOI_cols, slc.re);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_subarray_from_h5(h5_file, "s_im", offset_r, offset_c, AOI_rows, AOI_cols, slc.im);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	return 0;
}

int Utils::get_AOI_size(const char* h5_file, double lon_center, double lat_center, double width, double height, int* rows, int* cols, int* offset_row, int* offset_col)
{
	if (h5_file == NULL || rows == NULL || cols == NULL||
		fabs(lon_center) > 180.0 ||
		fabs(lat_center) > 180.0 ||
		width < 0.0 ||
		height < 0.0)
	{
		fprintf(stderr, "get_AOI_from_h5slc(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion;
	int ret, row_center, col_center, total_rows, total_cols, AOI_rows, AOI_cols, offset_r, offset_c;
	double inc_center, range_spacing, azimuth_spacing;
	Mat row_coef, col_coef, lon, lat, tmp;
	lon = Mat::zeros(1, 1, CV_64F); lon.at<double>(0, 0) = lon_center;
	lat = Mat::zeros(1, 1, CV_64F); lat.at<double>(0, 0) = lat_center;
	//读取图像行列总数
	ret = conversion.read_array_from_h5(h5_file, "azimuth_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	total_rows = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "range_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	total_cols = tmp.at<int>(0, 0);
	//读取采样间隔和下视角
	ret = conversion.read_array_from_h5(h5_file, "azimuth_spacing", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	azimuth_spacing = tmp.at<double>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "range_spacing", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	range_spacing = tmp.at<double>(0, 0);
	ret = conversion.read_array_from_h5(h5_file, "incidence_center", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	inc_center = tmp.at<double>(0, 0);
	//确定AOI中心图像坐标
	ret = conversion.read_array_from_h5(h5_file, "row_coefficient", row_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5_file, "col_coefficient", col_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = coord_conversion(row_coef, lon, lat, tmp);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	row_center = (int)floor(tmp.at<double>(0, 0));
	ret = coord_conversion(col_coef, lon, lat, tmp);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	col_center = (int)floor(tmp.at<double>(0, 0));
	if (row_center < 1 || row_center > total_rows - 1 || col_center < 1 || col_center > total_cols - 1)
	{
		fprintf(stderr, "get_AOI_from_h5slc(): AOI not found in %s!\n", h5_file);
		return -1;
	}

	AOI_rows = (int)floor(height / azimuth_spacing);
	AOI_cols = (int)floor(width / range_spacing * sin(inc_center / 180.0 * PI));

	AOI_rows = AOI_rows > total_rows ? total_rows : AOI_rows;
	AOI_cols = AOI_cols > total_cols ? total_cols : AOI_cols;

	offset_r = row_center - (int)AOI_rows / 2;
	offset_r = offset_r < 0 ? 0 : offset_r;
	offset_c = col_center - (int)AOI_cols / 2;
	offset_c = offset_c < 0 ? 0 : offset_c;
	AOI_rows = (total_rows - offset_r) > AOI_rows ? AOI_rows : (total_rows - offset_r);
	AOI_cols = (total_cols - offset_c) > AOI_cols ? AOI_cols : (total_cols - offset_c);
	*rows = AOI_rows;
	*cols = AOI_cols;
	if (offset_row) *offset_row = offset_r;
	if (offset_col)*offset_col = offset_c;
	return 0;
}

int Utils::coord_conversion(Mat& coefficient, Mat& coord_in_1, Mat& coord_in_2, Mat& coord_out)
{
	if (coefficient.cols != 32 ||
		coefficient.rows != 1 ||
		coefficient.type() != CV_64F ||
		coord_in_1.empty() ||
		coord_in_2.empty() ||
		coord_in_1.type() != CV_64F ||
		coord_in_2.type() != CV_64F||
		coord_in_1.rows != coord_in_2.rows||
		coord_in_1.cols != coord_in_2.cols)
	{
		fprintf(stderr, "coord_conversion(): input check failed!\n");
		return -1;
	}
	double offset_out, scale_out, offset_in_1, offset_in_2, scale_in_1, scale_in_2;
	double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11,
		a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24;
	offset_out = coefficient.at<double>(0, 0);
	scale_out = coefficient.at<double>(0, 1);
	offset_in_1 = coefficient.at<double>(0, 2);
	scale_in_1 = coefficient.at<double>(0, 3);
	offset_in_2 = coefficient.at<double>(0, 4);
	scale_in_2 = coefficient.at<double>(0, 5);

	a0 = coefficient.at<double>(0, 6);
	a1 = coefficient.at<double>(0, 7);
	a2 = coefficient.at<double>(0, 8);
	a3 = coefficient.at<double>(0, 9);
	a4 = coefficient.at<double>(0, 10);
	a5 = coefficient.at<double>(0, 11);
	a6 = coefficient.at<double>(0, 12);
	a7 = coefficient.at<double>(0, 13);
	a8 = coefficient.at<double>(0, 14);
	a9 = coefficient.at<double>(0, 15);
	a10 = coefficient.at<double>(0, 16);
	a11 = coefficient.at<double>(0, 17);
	a12 = coefficient.at<double>(0, 18);
	a13 = coefficient.at<double>(0, 19);
	a14 = coefficient.at<double>(0, 20);
	a15 = coefficient.at<double>(0, 21);
	a16 = coefficient.at<double>(0, 22);
	a17 = coefficient.at<double>(0, 23);
	a18 = coefficient.at<double>(0, 24);
	a19 = coefficient.at<double>(0, 25);
	a20 = coefficient.at<double>(0, 26);
	a21 = coefficient.at<double>(0, 27);
	a22 = coefficient.at<double>(0, 28);
	a23 = coefficient.at<double>(0, 29);
	a24 = coefficient.at<double>(0, 30);
	int rows = coord_in_1.rows; int cols = coord_in_1.cols;
	Mat coord_out1(rows, cols, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		double in1, in2, out;
		for (int j = 0; j < cols; j++)
		{
			in1 = (coord_in_1.at<double>(i, j) - offset_in_1) / scale_in_1;
			in2 = (coord_in_2.at<double>(i, j) - offset_in_2) / scale_in_2;
			out = a0 + a1 * in1 + a2 * in1 * in1 + a3 * in1 * in1 * in1 + a4 * in1 * in1 * in1 * in1 +
				a5 * in2 + a6 * in2 * in1 + a7 * in2 * in1 * in1 + a8 * in2 * in1 * in1 * in1 + a9 * in2 * in1 * in1 * in1 * in1 +
				a10 * in2 * in2 + a11 * in2 * in2 * in1 + a12 * in2 * in2 * in1 * in1 + a13 * in2 * in2 * in1 * in1 * in1 + a14 * in2 * in2 * in1 * in1 * in1 * in1 +
				a15 * in2 * in2 * in2 + a16 * in2 * in2 * in2 * in1 + a17 * in2 * in2 * in2 * in1 * in1 + a18 * in2 * in2 * in2 * in1 * in1 * in1 + a19 * in2 * in2 * in2 * in1 * in1 * in1 * in1 +
				a20 * in2 * in2 * in2 * in2 + a21 * in2 * in2 * in2 * in2 * in1 + a22 * in2 * in2 * in2 * in2 * in1 * in1 + a23 * in2 * in2 * in2 * in2 * in1 * in1 * in1 + a24 * in2 * in2 * in2 * in2 * in1 * in1 * in1 * in1;
			out = out * scale_out + offset_out;
			coord_out1.at<double>(i, j) = out;
		}
	}
	coord_out1.copyTo(coord_out);
	return 0;
}

int Utils::baseline_estimation(
	const Mat& stateVec1,
	const Mat& stateVec2, 
	const Mat& lon_coef,
	const Mat& lat_coef, 
	int offset_row,
	int offset_col,
	int scene_height, 
	int scene_width,
	double time_interval,
	double time_interval2,
	double* B_effect,
	double* B_parallel,
	double* sigma_B_effect, 
	double* sigma_B_parallel
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.type() != CV_64F ||
		lon_coef.cols != 32 ||
		lon_coef.rows != 1 ||
		lon_coef.type() != CV_64F ||
		lat_coef.cols != 32 ||
		lat_coef.rows != 1 ||
		lat_coef.type() != CV_64F ||
		offset_row < 0 ||
		offset_col < 0 ||
		B_effect == NULL ||
		B_parallel == NULL ||
		time_interval < 0.0 ||
		time_interval2 < 0.0 ||
		scene_height < 7 ||
		scene_width < 0
		)
	{
		fprintf(stderr, "baseline_estimation(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion;
	int ret;
	int rows = scene_height; int cols = scene_width;
	/*
	* 轨道插值
	*/
	Mat state_vec1, state_vec2;
	stateVec1.copyTo(state_vec1);
	stateVec2.copyTo(state_vec2);
	ret = stateVec_interp(state_vec1, time_interval, state_vec1);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	ret = stateVec_interp(state_vec2, time_interval2, state_vec2);
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
	ret = coord_conversion(lon_coefficient, row, col, lon);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;
	ret = coord_conversion(lat_coefficient, row, col, lat);
	if (return_check(ret, "coord_conversion()", error_head)) return -1;


	/*
	* 图像1成像点位置计算
	*/

	Mat sate1 = Mat::zeros(rows, 3, CV_64F);
	Mat sate2 = Mat::zeros(rows, 3, CV_64F);
	Mat satev1 = Mat::zeros(rows, 3, CV_64F);
	Mat satev2 = Mat::zeros(rows, 3, CV_64F);
	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = 0;
		ell2xyz(tmp, xyz);

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
			sate1_v(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(satev1(Range(j, j + 1), Range(0, 3)));
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
		ell2xyz(tmp, xyz);

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
			sate2_v(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(satev2(Range(j, j + 1), Range(0, 3)));
		}
	}

	/*
	*估计基线 
	*/
	Mat B_effe(rows, 1, CV_64F); Mat B_para(rows, 1, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		Mat R, B, tmp, xyz, effect_dir; double r;
		tmp = Mat::zeros(1, 3, CV_64F);
		tmp.at<double>(0, 0) = lat.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 1) = lon.at<double>(i, (int)cols / 2);
		tmp.at<double>(0, 2) = 0;
		ell2xyz(tmp, xyz);

		R = xyz - sate1(Range(i, i + 1), Range(0, 3));
		r = sqrt(sum(R.mul(R))[0]);
		R = R / r;
		B = sate2(Range(i, i + 1), Range(0, 3)) - sate1(Range(i, i + 1), Range(0, 3));
		B_para.at<double>(i, 0) = sum(R.mul(B))[0];//平行基线

		tmp = satev1(Range(i, i + 1), Range(0, 3));
		cross(tmp, R, effect_dir);
		r = sqrt(sum(effect_dir.mul(effect_dir))[0]);
		effect_dir = effect_dir / r;
		r = sqrt(sum(xyz.mul(xyz))[0]);
		xyz = xyz / r;
		r = sum(xyz.mul(effect_dir))[0];
		effect_dir = r < 0.0 ? -effect_dir : effect_dir;
		B_effe.at<double>(i, 0) = sum(effect_dir.mul(B))[0];//平行基线
	}
	*B_effect = sum(B_effe)[0] / (double)rows;
	*B_parallel = sum(B_para)[0] / (double)rows;
	if (sigma_B_effect)
	{
		this->std(B_effe, sigma_B_effect);
	}
	if (sigma_B_parallel)
	{
		this->std(B_para, sigma_B_parallel);
	}
	return 0;
}

int Utils::homogeneous_test(const Mat& pixel1, const Mat& pixel2, int* homo_flag, double alpha, const char* method)
{
	if (pixel1.cols != 1 ||
		pixel1.rows < 5 ||
		pixel1.type() != CV_64F ||
		pixel1.rows != pixel2.rows ||
		pixel1.cols != pixel2.cols ||
		pixel2.type() != CV_64F ||
		homo_flag == NULL ||
		method == NULL
		)
	{
		fprintf(stderr, "homogeneous_test(): input check failed!\n");
		return -1;
	}

	//Kolmogorov-Smirnov检验
	/*
	 alpha      0.20	0.15	0.10	0.05	0.025	0.01	0.005	0.001
     c(alpha)   1.073	1.138	1.224	1.358	1.48	1.628	1.731	1.949
	*/
	if (strcmp(method, "KS") == 0)
	{
		Mat p1, p2, cdf1, cdf2;
		double thresh;
		pixel1.copyTo(p1); pixel2.copyTo(p2);
		int N = p1.rows;
		cv::sort(p1, p1, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
		cv::sort(p2, p2, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
		if (p1.at<double>(N - 2, 0) <= p2.at<double>(0, 0) || p2.at<double>(N - 2, 0) <= p1.at<double>(0, 0))
		{
			*homo_flag = -1;
			return 0;
		}
		//确定threshold
		if (fabs(alpha - 0.2) < 0.01)
		{
			thresh = sqrt(2 / (double)N) * 1.073;
		}
		else if (fabs(alpha - 0.15) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.138;
		}
		else if (fabs(alpha - 0.1) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.224;
		}
		else if (fabs(alpha - 0.05) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.358;
		}
		else if (fabs(alpha - 0.025) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.48;
		}
		else if (fabs(alpha - 0.01) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.628;
		}
		else if (fabs(alpha - 0.005) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.731;
		}
		else if (fabs(alpha - 0.001) < 0.0001)
		{
			thresh = sqrt(2 / (double)N) * 1.949;
		}
		else
		{
			thresh = sqrt(2 / (double)N) * 1.358;
		}
		//计算C.D.F最大间距
		double Dmax = 0.0, tmp;
		int front_1 = 0, front_2 = 0;
		if (p1.at<double>(0, 0) > p2.at<double>(0, 0))
		{
			for (int i = 0; i < N - 1;i++)
			{
				if (p2.at<double>(i, 0) <= p1.at<double>(0, 0) && p2.at<double>(i + 1, 0) >= p1.at<double>(0, 0))
				{
					front_2 = i;
					break;
				}
			}
		}
		else
		{
			for (int i = 0; i < N - 1; i++)
			{
				if (p1.at<double>(i, 0) <= p2.at<double>(0, 0) && p1.at<double>(i + 1, 0) >= p2.at<double>(0, 0))
				{
					front_1 = i;
					break;
				}
			}
		}
		while (front_1 < N && front_2 < N)
		{
			tmp = fabs((double)front_1 / (double)N - (double)front_2 / (double)N);
			Dmax = Dmax > tmp ? Dmax : tmp;
			if (front_1 >= N - 1 || front_2 >= N - 1)
			{
				break;
			}
			if (p1.at<double>(front_1 + 1, 0) < p2.at<double>(front_2 + 1, 0)) front_1++;
			else if (p1.at<double>(front_1 + 1, 0) > p2.at<double>(front_2 + 1, 0)) front_2++;
			else
			{
				front_1++; front_2++;
			}
		}
		if (Dmax > thresh) *homo_flag = -1;
		else *homo_flag = 0;
	}

	//Anderson-Darling 检验
	else if (strcmp(method, "AD")== 0)
	{
		/*
		alpha = 0.01, AD_inf = 3.857;
		alpha = 0.05, AD_inf = 2.492;
		alpha = 0.1,  AD_inf = 1.933;
		*/
		Mat p1, p2, p;
		double thresh;
		pixel1.copyTo(p1); pixel2.copyTo(p2);
		cv::vconcat(p1, p2, p);
		int n = p1.rows; int N = 2 * n;
		cv::sort(p1, p1, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
		cv::sort(p2, p2, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
		cv::sort(p, p, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
		if (fabs(alpha - 0.05) < 0.0001)
		{
			thresh = (2.492 - 1) * (1 - 1.55 / (double)N) + 1;
		}
		else if (fabs(alpha - 0.01) < 0.0001)
		{
			thresh = (3.857 - 1) * (1 - 1.55 / (double)N) + 1;
		}
		else
		{
			thresh = (1.933 - 1) * (1 - 1.55 / (double)N) + 1;
		}
		int c = 0; double sum = 0.0, sentinel = -1.0;
		for (int i = 1; i < N; i++)
		{
			while (c <= n - 1)
			{
				if (p.at<double>(i - 1, 0) <= p1.at<double>(c, 0)) break;
				c++;
			}
			sum += double((N * c - n * i) * (N * c - n * i)) / double(i * (N - i));
		}
		sum /= (double)(n * n);
		if (sum > thresh) *homo_flag = -1;
		else *homo_flag = 0;

	}

	return 0;
}

int Utils::HermitianEVD(const ComplexMat& input, Mat& eigenvalue, ComplexMat& eigenvector)
{
	if (input.GetCols() != input.GetRows() ||
		input.GetCols() <= 1 || 
		input.type() != CV_64F
		)
	{
		fprintf(stderr, "HermitianEVD(): input check failed!\n");
		return -1;
	}
	int rows = input.GetRows(); int cols = input.GetCols();
	eigenvalue.create(cols, 1, CV_64F);
	eigenvector.re.create(rows, cols, CV_64F);
	eigenvector.im.create(rows, cols, CV_64F);
	Eigen::MatrixXcd x(rows, cols);
	complex<double> d;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			d.real(input.re.at<double>(i, j));
			d.imag(input.im.at<double>(i, j));
			x(i, j) = d;
		}
	}
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
	solver.compute(x, true);
	if (solver.info() != Eigen::Success)
	{
		fprintf(stderr, "ComplexEVD(): EVD failed!\n");
		return -1;
	}
	Eigen::MatrixXcd eigenvectors = solver.eigenvectors();
	Eigen::MatrixXcd eigenvalues = solver.eigenvalues();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			eigenvector.re.at<double>(i, j) = eigenvectors(i, j).real();
			eigenvector.im.at<double>(i, j) = eigenvectors(i, j).imag();
			eigenvalue.at<double>(i, 0) = eigenvalues(i, 0).real();
		}
	}
	Mat idx, value; ComplexMat vec; int xx, xxx;
	eigenvalue.copyTo(value); vec = eigenvector;
	cv::sortIdx(eigenvalue, idx, cv::SORT_EVERY_COLUMN + cv::SORT_DESCENDING);
	for (int i = 0; i < rows; i++)
	{
		xx = idx.at<int>(i, 0);
		for (int j = 0; j < cols; j++)
		{
			xxx = idx.at<int>(j, 0);
			eigenvalue.at<double>(i, 0) = value.at<double>(xx, 0);
			eigenvector.re.at<double>(i, j) = vec.re.at<double>(i, xxx);
			eigenvector.im.at<double>(i, j) = vec.im.at<double>(i, xxx);
		}
	}
	return 0;
}

int Utils::coherence_matrix_estimation(const vector<ComplexMat>& slc_series, ComplexMat& coherence_matrix, int est_window_width, int est_window_height,  int ref_row, int ref_col, bool b_homogeneous_test, bool b_normalize)
{
	if (slc_series.size() < 5||
		est_window_width < 3||
		est_window_height < 3||
		est_window_height % 2 != 1||
		est_window_width % 2 != 1||
		ref_col < 0||
		ref_row < 0
		)
	{
		fprintf(stderr, "coherence_estimation(): input check failed!\n");
		return -1;
	}
	if (slc_series[0].type() != CV_64F || slc_series[0].isempty() || ref_row > slc_series[0].GetRows() - 1 || ref_col > slc_series[0].GetCols() - 1
		)
	{
		fprintf(stderr, "coherence_estimation(): input check failed!\n");
		return -1;
	}
	int n_images = slc_series.size(), ret;
	int rows = slc_series[0].GetRows(); int cols = slc_series[0].GetCols();
	int radius_width = (est_window_width - 1) / 2;
	int radius_height = (est_window_height - 1) / 2;
	int left, right, bottom, top;
	left = (ref_col - radius_width) < 0 ? 0 : (ref_col - radius_width);
	right = (ref_col + radius_width) > cols - 1 ? cols - 1 : (ref_col + radius_width);
	bottom = (ref_row + radius_height) > rows - 1 ? rows - 1 : (ref_row + radius_height);
	top = (ref_row - radius_height) < 0 ? 0 : (ref_row - radius_height);
	rows = bottom - top + 1;
	cols = right - left + 1;
	if (b_homogeneous_test)
	{
		//统计同质检验
		ComplexMat pix1(n_images, 1); ComplexMat pix2(n_images, 1);
		Mat pix1_amp, pix2_amp;
		Mat mask = Mat::zeros(rows, cols, CV_32S);
		mask.at<int>(ref_row - top, ref_col - left) = 1;
		int b_homo, count = 1;
		for (int i = 0; i < n_images; i++)
		{
			pix1.re.at<double>(i, 0) = slc_series[i].re.at<double>(ref_row, ref_col);
			pix1.im.at<double>(i, 0) = slc_series[i].im.at<double>(ref_row, ref_col);
		}
		pix1_amp = pix1.GetMod();
		for (int i = 0; i < rows; i++)
		{
			
			for (int j = 0; j < cols; j++)
			{
				if (i == (ref_row - top) && j == (ref_col - left)) continue;
				for (int k = 0; k < n_images; k++)
				{
					pix2.re.at<double>(k, 0) = slc_series[k].re.at<double>(i + top, j + left);
					pix2.im.at<double>(k, 0) = slc_series[k].im.at<double>(i + top, j + left);
				}
				pix2_amp = pix2.GetMod();
				ret = homogeneous_test(pix1_amp, pix2_amp, &b_homo, 0.1);
				//ret = homogeneous_test(pix1_amp, pix2_amp, &b_homo, 0.1, "AD");
				if (return_check(ret, "homogeneous_test()", error_head)) return -1;
				if (b_homo == 0) { mask.at<int>(i, j) = 1; count++; }
			}
		}
		//Mat c; mask.convertTo(c, CV_64F);
		//cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\mask.bin", c);
		if (count < 2) 
		{
			//fprintf(stderr, "coherence_matrix_estimation(): no homogenous pixels inside estimation window!\n");
			return -1;
		}
		//估计相关矩阵
		int count2 = count; count = 0;
		ComplexMat Covariance;
		Mat sum(n_images, 1, CV_64F), A(n_images, count2, CV_64F), B(n_images, count2, CV_64F), C, A_t, B_t;
		double s;
		for (int k = 0; k < n_images; k++)
		{
			s = 0.0;
			count = 0;
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					if (mask.at<int>(i, j) > 0)
					{
						A.at<double>(k, count) = slc_series[k].re.at<double>(i + top, j + left);
						B.at<double>(k, count) = slc_series[k].im.at<double>(i + top, j + left);
						count++;
						s += A.at<double>(k, count - 1) * A.at<double>(k, count - 1)
							+ B.at<double>(k, count - 1) * B.at<double>(k, count - 1);
					}
				}
			}
			sum.at<double>(k, 0) = s;
		}
		cv::transpose(A, A_t); cv::transpose(B, B_t);
		C = A * A_t + B * B_t;
		C.copyTo(Covariance.re);
		C = B * A_t - A * B_t;
		C.copyTo(Covariance.im);
		double denum;
		if (b_normalize)
		{

			for (int i = 0; i < n_images; i++)
			{
				for (int j = 0; j < n_images; j++)
				{
					denum = sqrt(sum.at<double>(i, 0) * sum.at<double>(j, 0));
					Covariance.re.at<double>(i, j) = Covariance.re.at<double>(i, j) / (denum + 1e-10);
					Covariance.im.at<double>(i, j) = Covariance.im.at<double>(i, j) / (denum + 1e-10);
				}
			}


		}
		else
		{
			Covariance = Covariance * (1 / (double)count2);
		}
		coherence_matrix = Covariance;
	}
	else
	{

	}
	return 0;
}

int Utils::MB_phase_estimation(
	vector<string> coregis_slc_files,
	vector<string> phase_files, 
	vector<string> coherence_files,
	int master_indx, 
	int blocksize_row, 
	int blocksize_col, 
	Mat& out_mask,
	bool b_coh_est,
	int homogeneous_test_wnd,
	double thresh_c1_to_c2,
	bool b_flat,
	bool b_normalize
)
{
	if (coregis_slc_files.size() < 2 ||
		phase_files.size() != coregis_slc_files.size() ||
		coherence_files.size() != phase_files.size() ||
		master_indx < 1 ||
		master_indx > phase_files.size() ||
		blocksize_row < 100 ||
		blocksize_col < 100 ||
		thresh_c1_to_c2 < 0.0 ||
		thresh_c1_to_c2 > 1.0 ||
		homogeneous_test_wnd % 2 == 0
		)
	{
		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
		return -1;
	}
	int homotest_radius = (homogeneous_test_wnd - 1) / 2;
	if (blocksize_row <= homotest_radius || blocksize_col <= homotest_radius)
	{
		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
		return -1;
	}
	int nr, nc, ret, n_images; Mat tmp;
	n_images = coregis_slc_files.size();
	FormatConversion conversion; Deflat flat;
	ret = conversion.read_array_from_h5(coregis_slc_files[0].c_str(), "azimuth_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nr = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(coregis_slc_files[0].c_str(), "range_len", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nc = tmp.at<int>(0, 0);
	Mat mask = Mat::zeros(nr, nc, CV_32S); mask.copyTo(out_mask); mask.release();
	//检查输入SAR图像尺寸是否相同
	for (int i = 1; i < n_images; i++)
	{
		ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "azimuth_len", tmp);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		if (tmp.at<int>(0, 0) != nr)
		{
			fprintf(stderr, "%s images size mismatch!\n", coregis_slc_files[i].c_str());
			return -1;
		}
		ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "range_len", tmp);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		if (tmp.at<int>(0, 0) != nc)
		{
			fprintf(stderr, "%s images size mismatch!\n", coregis_slc_files[i].c_str());
			return -1;
		}
	}

	//去平地相位
	Mat stateVec1, prf, lon_coef, lat_coef, 
		carrier_frequency, stateVec2, prf2, phase,
		phase_deflat, flat_phase_coef, azimuth_len, range_len;
	phase = Mat::zeros(nr, nc, CV_64F);
	azimuth_len = Mat::zeros(1, 1, CV_32S); range_len = Mat::zeros(1, 1, CV_32S);
	int offset_row, offset_col;
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "offset_row", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_row = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "offset_col", tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_col = tmp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "state_vec", stateVec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "prf", prf);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "lon_coefficient", lon_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "lat_coefficient", lat_coef);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "carrier_frequency", carrier_frequency);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	for (int i = 0; i < n_images; i++)
	{
		//预先填充相关系数
		if (b_coh_est)
		{
			ret = conversion.creat_new_h5(coherence_files[i].c_str());
			if (return_check(ret, "creat_new_h5()", error_head)) return -1;
			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "coherence", phase);
			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
			azimuth_len.at<int>(0, 0) = phase.rows; range_len.at<int>(0, 0) = phase.cols;
			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "azimuth_len", azimuth_len);
			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "range_len", range_len);
			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		}
		//预先填充干涉相位
		ret = conversion.creat_new_h5(phase_files[i].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "phase", phase);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		azimuth_len.at<int>(0, 0) = phase.rows; range_len.at<int>(0, 0) = phase.cols;
		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "azimuth_len", azimuth_len);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "range_len", range_len);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		if (b_flat)
		{
			ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "state_vec", stateVec2);
			if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
			ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "prf", prf2);
			if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
			ret = flat.deflat(stateVec1, stateVec2, lon_coef, lat_coef, phase, offset_row, offset_col, 0, 1 / prf.at<double>(0, 0),
				1 / prf2.at<double>(0, 0), 1, 3e8 / carrier_frequency.at<double>(0, 0), phase, flat_phase_coef);
			ret = conversion.write_array_to_h5(phase_files[i].c_str(), "flat_phase_coefficient", flat_phase_coef);
			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		}
		fprintf(stdout, "去平地进度：%d/%d\n", i, n_images - 1);
	}

	//分块读取、计算和储存

	int left, right, top, bottom, block_num_row, block_num_col, left_pad, right_pad, top_pad, bottom_pad;
	vector<ComplexMat> slc_series, slc_series_filter;
	vector<Mat> coherence_series; coherence_series.resize(n_images);
	ComplexMat slc, slc2, temp;
	Mat flat_phase, ph, zeromat;
	if (nr % blocksize_row == 0) block_num_row = nr / blocksize_row;
	else block_num_row = int(floor((double)nr / (double)blocksize_row)) + 1;
	if (nc % blocksize_col == 0) block_num_col = nc / blocksize_col;
	else block_num_col = int(floor((double)nc / (double)blocksize_col)) + 1;
	for (int i = 0; i < block_num_row; i++)
	{
		for (int j = 0; j < block_num_col; j++)
		{
			top = i * blocksize_row;
			top_pad = top - homotest_radius; top_pad = top_pad < 0 ? 0 : top_pad;
			bottom = top + blocksize_row; bottom = bottom > nr ? nr : bottom;
			bottom_pad = bottom + homotest_radius; bottom_pad = bottom_pad > nr ? nr : bottom_pad;
			left = j * blocksize_col;
			left_pad = left - homotest_radius; left_pad = left_pad < 0 ? 0 : left_pad;
			right = left + blocksize_col; right = right > nc ? nc : right;
			right_pad = right + homotest_radius; right_pad = right_pad > nc ? nc : right_pad;

			//读取数据
			for (int k = 0; k < n_images; k++)
			{
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_re",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.re);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_im",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.im);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				if (b_flat)
				{
					if (k != master_indx - 1)
					{
						flat_phase.create(bottom_pad - top_pad, right_pad - left_pad, CV_64F);
						ret = conversion.read_array_from_h5(phase_files[k].c_str(), "flat_phase_coefficient", flat_phase_coef);
						if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
#pragma omp parallel for schedule(guided)
						for (int ii = top_pad; ii < bottom_pad; ii++)
						{
							Mat tempp(1, 6, CV_64F);
							for (int jj = left_pad; jj < right_pad; jj++)
							{
								tempp.at<double>(0, 0) = 1.0;
								tempp.at<double>(0, 1) = ii;
								tempp.at<double>(0, 2) = jj;
								tempp.at<double>(0, 3) = ii * jj;
								tempp.at<double>(0, 4) = ii * ii;
								tempp.at<double>(0, 5) = jj * jj;
								flat_phase.at<double>(ii - top_pad, jj - left_pad) = sum(tempp.mul(flat_phase_coef))[0];
							}
						}
						ret = phase2cos(flat_phase, temp.re, temp.im);
						if (return_check(ret, "phase2cos()", error_head)) return -1;
						slc = slc * temp;
					}
					
				}
				slc_series.push_back(slc);
				slc_series_filter.push_back(slc);
			}

			//填充相关系数
			if (b_coh_est)
			{
				zeromat = Mat::zeros(flat_phase.rows, flat_phase.cols, CV_64F);
				for (int mm = 0; mm < n_images; mm++)
				{
					zeromat.copyTo(coherence_series[mm]);
				}
			}
			//计算
#pragma omp parallel for schedule(guided)
			for (int ii = (top - top_pad); ii < (bottom - top_pad); ii++)
			{
				ComplexMat coherence_matrix, eigenvector; Mat eigenvalue; int ret;
				for (int jj = (left - left_pad); jj < (right - left_pad); jj++)
				{
					ret = coherence_matrix_estimation(slc_series, coherence_matrix, homogeneous_test_wnd, homogeneous_test_wnd, ii, jj);
					if (ret == 0)
					{
						ret = HermitianEVD(coherence_matrix, eigenvalue, eigenvector);
						if (!eigenvalue.empty() && ret == 0)
						{
							if (eigenvalue.at<double>(1, 0) / (eigenvalue.at<double>(0, 0) + 1e-10) < thresh_c1_to_c2)
							{
								out_mask.at<int>(ii + top_pad, jj + left_pad) = 1;
								for (int kk = 0; kk < n_images; kk++)
								{
									slc_series_filter[kk].re.at<double>(ii, jj) = eigenvector.re.at<double>(kk, 0);
									slc_series_filter[kk].im.at<double>(ii, jj) = eigenvector.im.at<double>(kk, 0);
									if (b_coh_est)
									{
										coherence_series[kk].at<double>(ii, jj) = coherence_matrix(cv::Range(master_indx - 1, master_indx),
											cv::Range(kk, kk + 1)).GetMod().at<double>(0, 0);
									}
									
								}
							}
						}
					}


				}
			}

			//储存
			slc = slc_series_filter[master_indx - 1];
			for (int kk = 0; kk < n_images; kk++)
			{
				coherence_series[kk](cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ret = conversion.write_subarray_to_h5(coherence_files[kk].c_str(), "coherence", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
				//if (kk == master_indx - 1) continue;
				ret = multilook(slc, slc_series_filter[kk], 1, 1, phase);
				if (return_check(ret, "multilook()", error_head)) return -1;
				phase(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ret = conversion.write_subarray_to_h5(phase_files[kk].c_str(), "phase", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
				
			}
			slc_series.clear();
			slc_series_filter.clear();
			
			fprintf(stdout, "估计相位进度：%lf\n", double((i + 1) * block_num_col + j + 1) / double((block_num_col) * (block_num_row)));
		}
	}

	return 0;
}

int Utils::unwrap_region_growing(
	vector<tri_node>& nodes, 
	const vector<tri_edge>& edges, 
	size_t start_edge, 
	double distance_thresh, 
	double quality_thresh
)
{
	if (nodes.size() < 3 ||
		edges.size() < 3 ||
		start_edge < 1 ||
		start_edge > edges.size()
		)
	{
		fprintf(stderr, "unwrap_region_growing(): input check failed!\n\n");
		return -1;
	}
	if (distance_thresh < 1.0) distance_thresh = 1.0;
	if (quality_thresh > 0.9) quality_thresh = 0.9;
	//找到增量积分起始点
	int ix = 0;
	double MC = -1.0;
	size_t num_edges = edges.size();
	//for (int i = 0; i < num_edges; i++)
	//{
	//	if (edges[i].MC > MC)
	//	{
	//		ix = i;
	//		MC = edges[i].MC;
	//	}
	//}
	size_t start = edges[start_edge - 1].end1;

	//采用类似质量图法解缠的算法进行增量积分集成
	edge_index tmp;
	priority_queue<edge_index> que;
	//nodes[start - 1].set_vel(0.0);//起始点形变速率和高程误差设置为0，后续可根据参考点进行校正
	//nodes[start - 1].set_height(0.0);
	nodes[start - 1].set_status(true);
	long* ptr_neigh = NULL;
	int num_neigh, end2, number, row1, col1, row2, col2;
	double distance, phase, MC_total, phase_total, delta_phase;

	nodes[start - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
	for (int i = 0; i < num_neigh; i++)
	{
		end2 = edges[*(ptr_neigh + i) - 1].end1 == start ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
		nodes[start - 1].get_distance(nodes[end2 - 1], &distance);
		if (!nodes[end2 - 1].get_status() &&
			distance <= distance_thresh &&
			edges[*(ptr_neigh + i) - 1].quality > quality_thresh
			)
		{
			tmp.num = *(ptr_neigh + i);
			tmp.quality = -edges[*(ptr_neigh + i) - 1].quality;
			que.push(tmp);
		}
	}

	while (que.size() != 0)
	{
		tmp = que.top();
		que.pop();
		if (nodes[edges[tmp.num - 1].end1 - 1].get_status())
		{
			number = edges[tmp.num - 1].end1;
			end2 = edges[tmp.num - 1].end2;
		}
		else
		{
			number = edges[tmp.num - 1].end2;
			end2 = edges[tmp.num - 1].end1;
		}
		MC_total = 1e-10;
		phase_total = 0.0;
		if (!nodes[end2 - 1].get_status())
		{
			nodes[end2 - 1].get_pos(&row2, &col2);
			nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			for (int i = 0; i < num_neigh; i++)
			{
				number = edges[*(ptr_neigh + i) - 1].end1 == end2 ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
				if (nodes[number - 1].get_status())
				{
					nodes[number - 1].get_phase(&phase);
					nodes[number - 1].get_pos(&row1, &col1);
					if (row1 > row2)
					{
						delta_phase = -edges[*(ptr_neigh + i) - 1].phase_diff;
					}
					if (row1 < row2)
					{
						delta_phase = edges[*(ptr_neigh + i) - 1].phase_diff;
					}
					if (row1 == row2)
					{
						if (col1 > col2)
						{
							delta_phase = -edges[*(ptr_neigh + i) - 1].phase_diff;
						}
						else
						{
							delta_phase = edges[*(ptr_neigh + i) - 1].phase_diff;
						}
					}
					phase_total += (delta_phase + phase) * edges[*(ptr_neigh + i) - 1].quality;
					MC_total += edges[*(ptr_neigh + i) - 1].quality;
				}
			}
			nodes[end2 - 1].set_phase(phase_total / MC_total);
			nodes[end2 - 1].set_status(true);
			nodes[end2 - 1].get_neigh_ptr(&ptr_neigh, &num_neigh);
			number = end2;
			for (int i = 0; i < num_neigh; i++)
			{

				end2 = edges[*(ptr_neigh + i) - 1].end1 == number ? edges[*(ptr_neigh + i) - 1].end2 : edges[*(ptr_neigh + i) - 1].end1;
				nodes[number - 1].get_distance(nodes[end2 - 1], &distance);
				if (!nodes[end2 - 1].get_status() &&
					distance <= distance_thresh &&
					edges[*(ptr_neigh + i) - 1].quality > quality_thresh
					)
				{
					tmp.num = *(ptr_neigh + i);
					tmp.quality = -edges[*(ptr_neigh + i) - 1].quality;
					que.push(tmp);
				}
			}
		}
	}
	return 0;
}

int Utils::unwrap_3D(
	const Mat& mask,
	const vector<Mat>& quality_map,
	vector<Mat>& wrapped_phase_series,
	vector<Mat>& unwrapped_phase_series,
	const char* delaunay_exe_path,
	const char* tmp_file_path,
	double distance_thresh,
	double quality_thresh
)
{
	if (mask.rows < 2 ||
		mask.cols < 2 ||
		mask.type() != CV_32S ||
		wrapped_phase_series.size() < 2||
		delaunay_exe_path == NULL||
		tmp_file_path == NULL ||
		quality_map.size() != wrapped_phase_series.size()
		)
	{
		fprintf(stderr, "unwrap_3D(): input check failed!\n");
		return -1;
	}

	quality_thresh = quality_thresh > 1.0 ? 0.9 : quality_thresh;
	distance_thresh = distance_thresh < 1.0 ? 1.0 : distance_thresh;

	/*----------------------------------*/
	/*            时间维解缠            */
	/*----------------------------------*/

	/*
	* 利用mask生成delaunay三角网络
	*/
	size_t nr = mask.rows;
	size_t nc = mask.cols;
	long num_nodes = cv::countNonZero(mask);
	if (num_nodes < 3)
	{
		fprintf(stderr, "unwrap_3D(): at least 3 nodes are needed!\n");
		return -1;
	}
	int ret;
	int n_images = wrapped_phase_series.size();
	for (size_t i = 0; i < n_images; i++)
	{
		if (quality_map[i].rows != nr || quality_map[i].cols != nc)
		{
			fprintf(stderr, "unwrap_3D(): quality map size mismatch!\n");
			return -1;
		}
		if (wrapped_phase_series[i].rows != nr || wrapped_phase_series[i].cols != nc)
		{
			fprintf(stderr, "unwrap_3D(): wrapped phase size mismatch!\n");
			return -1;
		}
	}
	unwrapped_phase_series.resize(n_images);
	vector<tri_node> nodes; vector<vector<tri_node>> nodes_vec; nodes_vec.resize(n_images);
	vector<tri_edge> edges; vector<vector<tri_edge>> edges_vec; edges_vec.resize(n_images);
	vector<int> node_neighbour;
	string tmp_folder(tmp_file_path);
	string node_file = tmp_folder + "\\triangle.node";
	ret = write_node_file(node_file.c_str(), mask);
	if (return_check(ret, "write_node_file()", error_head)) return -1;
	ret = gen_delaunay(node_file.c_str(), delaunay_exe_path);
	if (return_check(ret, "gen_delaunay()", error_head)) return -1;
	string edge_file = tmp_folder + "\\triangle.1.edge";
	ret = read_edges(edge_file.c_str(), edges, node_neighbour, num_nodes);
	if (return_check(ret, "read_edges()", error_head)) return -1;
	for (int i = 0; i < n_images; i++)
	{
		ret = init_tri_node(nodes, wrapped_phase_series[i], mask, edges, node_neighbour, num_nodes);
		if (return_check(ret, "init_tri_node()", error_head)) return -1;
		ret = init_edge_phase_diff(edges, nodes);
		if (return_check(ret, "init_edge_phase_diff()", error_head)) return -1;
		ret = init_edges_quality(quality_map[i], edges, nodes);
		if (return_check(ret, "init_edges_quality()", error_head)) return -1;
		nodes_vec[i] = nodes;
		edges_vec[i] = edges;
	}

	/*
	* 一维高斯滤波
	*/

	int Gaussian_radius = 2;
	double sigma = 1.0;
	Mat Gaussian_template = Mat::zeros(2 * Gaussian_radius + 1, 1, CV_64F);
	for (int i = 0; i < 2 * Gaussian_radius + 1; i++)
	{
		Gaussian_template.at<double>(i, 0) = exp(-(double(i - Gaussian_radius)) * (double(i - Gaussian_radius)) / (2.0 * sigma * sigma))
			/ (sigma * sqrt(PI * 2.0));
	}
	size_t num_edges = edges.size();
#pragma omp parallel for schedule(guided)
	for (long long i = 0; i < num_edges; i++)
	{
		Mat time_series = Mat::zeros(n_images, 1, CV_64F); Mat temp;
		for (size_t j = 0; j < n_images; j++)
		{
			time_series.at<double>(j, 0) = edges_vec[j][i].phase_diff;
		}
		cv::copyMakeBorder(time_series, time_series, Gaussian_radius, Gaussian_radius, 0, 0, cv::BORDER_REPLICATE);
		for (size_t j = Gaussian_radius; j < n_images + Gaussian_radius; j++)
		{
			time_series.at<double>(j, 0) = cv::sum(time_series(cv::Range(j - Gaussian_radius, j + Gaussian_radius + 1),
				cv::Range(0, 1)).mul(Gaussian_template))[0];
		}
		time_series(cv::Range(Gaussian_radius, n_images + Gaussian_radius), cv::Range(0, 1)).copyTo(temp);

		/*
		* 时间维解缠
		*/
		Mat temp2; temp.copyTo(temp2);
		for (size_t j = 1; j < n_images; j++)
		{
			double tmp = temp.at<double>(j, 0) - temp.at<double>(j - 1, 0);
			tmp = atan2(sin(tmp), cos(tmp));
			temp2.at<double>(j, 0) = temp2.at<double>(j - 1, 0) + tmp;
		}

		for (size_t j = 0; j < n_images; j++)
		{
			edges_vec[j][i].phase_diff = temp2.at<double>(j, 0);
		}
	}

	/*----------------------------------*/
	/*            空间维解缠            */
	/*----------------------------------*/
	
	/*
	* 确定解缠起始点
	*/
	double mean_quality, max_quality = -1.0;
	size_t ix_max;
	for (size_t i = 0; i < num_edges; i++)
	{
		mean_quality = 0.0;
		for (size_t j = 0; j < n_images; j++)
		{
			mean_quality += edges_vec[j][i].quality;
		}
		mean_quality /= (double)n_images;
		if (mean_quality > max_quality)
		{
			max_quality = mean_quality;
			ix_max = i;
		}
	}
//#pragma omp parallel for schedule(guided)
	for (size_t i = 0; i < n_images; i++)
	{
		//vector<tri_node> nodes; vector<tri_edge>edges;
		int row, col; double phase;
		Mat unwrapped_phase;
		wrapped_phase_series[i].copyTo(unwrapped_phase);
		nodes = nodes_vec[i];
		edges = edges_vec[i];
		ret = unwrap_region_growing(nodes, edges, ix_max, distance_thresh, quality_thresh);
		if (return_check(ret, "unwrap_region_growing()", error_head)) return -1;



		/*
		* 从节点中取出解缠相位
		*/

		for (size_t j = 0; j < num_nodes; j++)
		{
			if (nodes[j].get_status())
			{
				nodes[j].get_pos(&row, &col);
				nodes[j].get_phase(&phase);
				unwrapped_phase.at<double>(row, col) = phase;
			}
		}
		unwrapped_phase.copyTo(unwrapped_phase_series[i]);
	}
	
	/*----------------------------------*/
	/*          时间维再解缠            */
	/*----------------------------------*/

	
//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < nr; i++)
//	{
//		double tmp, tmp_old;
//		for (size_t j = 0; j < nc; j++)
//		{
//			if (mask.at<int>(i, j) > 0)
//			{
//				for (size_t k = 1; k < n_images; k++)
//				{
//					if (k == 1) tmp_old = unwrapped_phase_series[0].at<double>(i, j);
//					tmp = unwrapped_phase_series[k].at<double>(i, j) - tmp_old;
//					tmp = atan2(sin(tmp), cos(tmp));
//					tmp_old = unwrapped_phase_series[k].at<double>(i, j);
//					unwrapped_phase_series[k].at<double>(i, j) = unwrapped_phase_series[k - 1].at<double>(i, j) + tmp;
//				}
//			}
//			
//		}
//	}

	return 0;
}

int Utils::unwrap_3D_mcf(
	const Mat& mask, 
	const vector<Mat>& quality_map,
	vector<Mat>& wrapped_phase_series,
	vector<Mat>& unwrapped_phase_series,
	const char* delaunay_exe_path, 
	const char* mcf_exe_path, 
	const char* tmp_file_path,
	double distance_thresh
)
{
	if (mask.rows < 2 ||
		mask.cols < 2 ||
		mask.type() != CV_32S ||
		wrapped_phase_series.size() < 2 ||
		delaunay_exe_path == NULL ||
		tmp_file_path == NULL ||
		mcf_exe_path == NULL ||
		quality_map.size() != wrapped_phase_series.size()
		)
	{
		fprintf(stderr, "unwrap_3D_mcf(): input check failed!\n");
		return -1;
	}

	distance_thresh = distance_thresh < 1.0 ? 1.0 : distance_thresh;

	/*----------------------------------*/
	/*            时间维解缠            */
	/*----------------------------------*/

	/*
	* 利用mask生成delaunay三角网络
	*/
	size_t nr = mask.rows;
	size_t nc = mask.cols;
	long num_nodes = cv::countNonZero(mask);
	if (num_nodes < 3)
	{
		fprintf(stderr, "unwrap_3D_mcf(): at least 3 nodes are needed!\n");
		return -1;
	}
	int ret;
	int n_images = wrapped_phase_series.size();
	for (size_t i = 0; i < n_images; i++)
	{
		if (quality_map[i].rows != nr || quality_map[i].cols != nc)
		{
			fprintf(stderr, "unwrap_3D_mcf(): quality map size mismatch!\n");
			return -1;
		}
		if (wrapped_phase_series[i].rows != nr || wrapped_phase_series[i].cols != nc)
		{
			fprintf(stderr, "unwrap_3D_mcf(): wrapped phase size mismatch!\n");
			return -1;
		}
	}
	Unwrap unwrap;
	unwrapped_phase_series.resize(n_images);
	vector<tri_node> nodes; 
	vector<tri_edge> edges; 
	vector<triangle> tri;
	vector<int> node_neighbour;
	Mat out_mask;
	string tmp_folder(tmp_file_path);
	string node_file = tmp_folder + "\\triangle.node";
	string ele_file = tmp_folder + "\\triangle.1.ele";
	string neigh_file = tmp_folder + "\\triangle.1.neigh";
	string mcf_problem = tmp_folder + "\\mcf_delaunay.net";
	string mcf_solution = tmp_folder + "\\mcf_delaunay.net.sol";
	string edge_file = tmp_folder + "\\triangle.1.edge";
	ret = write_node_file(node_file.c_str(), mask);
	if (return_check(ret, "write_node_file()", error_head)) return -1;
	ret = gen_delaunay(node_file.c_str(), delaunay_exe_path);
	if (return_check(ret, "gen_delaunay()", error_head)) return -1;
	ret = read_edges(edge_file.c_str(), edges, node_neighbour, num_nodes);
	if (return_check(ret, "read_edges()", error_head)) return -1;
	int num_triangle, positive = 0, negative = 0;
	size_t num_edges = edges.size();
	//先将起始点做时间维解缠
	int pos_row, pos_col; double tmp;
	ret = init_tri_node(nodes, wrapped_phase_series[0], mask, edges, node_neighbour, num_nodes);
	if (return_check(ret, "init_tri_node()", error_head)) return -1;
	nodes[199].get_pos(&pos_row, &pos_col);
	Mat time_series(n_images, 1, CV_64F);
	for (int i = 0; i < n_images; i++)
	{
		time_series.at<double>(i, 0) = wrapped_phase_series[i].at<double>(pos_row, pos_col);
	}
	for (int i = 1; i < n_images; i++)
	{
		tmp = time_series.at<double>(i, 0) - time_series.at<double>(i - 1, 0);
		wrapped_phase_series[i].at<double>(pos_row, pos_col) = time_series.at<double>(i - 1, 0) +
			atan2(sin(tmp), cos(tmp));
	}

	for (int i = 0; i < n_images; i++)
	{
		//将edges的gain清零
		for (size_t ii = 0; ii < num_edges; ii++)
		{
			edges[ii].gain = 0.0;
		}
		positive = 0, negative = 0;
		ret = init_tri_node(nodes, wrapped_phase_series[i], mask, edges, node_neighbour, num_nodes);
		if (return_check(ret, "init_tri_node()", error_head)) return -1;
		if (i == 0)
		{
			ret = read_triangle(ele_file.c_str(), neigh_file.c_str(), tri, nodes, edges);
			if (return_check(ret, "read_triangle()", error_head)) return -1;
			num_triangle = tri.size();
		}
		ret = residue(tri, nodes, edges, 10.0);
		if (return_check(ret, "residue()", error_head)) return -1;
		/*
		* 检查残差点数，若无残差点则不使用mcf.exe求解
		*/
		for (int ii = 0; ii < num_triangle; ii++)
		{
			if (tri[ii].residue > 0.7)
			{
				positive++;
			}
			if (tri[ii].residue < -0.7)
			{
				negative++;
			}
		}
		if (positive == 0 && negative == 0)
		{
			
			ret = unwrap.MCF(wrapped_phase_series[i], unwrapped_phase_series[i], out_mask, mask, nodes, edges, 200, false, distance_thresh);
			if (return_check(ret, "unwrap.MCF()", error_head)) return -1;
		}
		else
		{
			ret = write_DIMACS(mcf_problem.c_str(), tri, nodes, edges, quality_map[i]);
			if (return_check(ret, "write_DIMACS()", error_head)) return -1;
			ret = unwrap.mcf_delaunay(mcf_problem.c_str(), mcf_exe_path);
			if (return_check(ret, "unwrap.mcf_delaunay()", error_head)) return -1;
			ret = read_DIMACS(mcf_solution.c_str(), edges, nodes, tri);
			if (return_check(ret, "read_DIMACS()", error_head)) return -1;
			ret = unwrap.MCF(wrapped_phase_series[i], unwrapped_phase_series[i], out_mask, mask, nodes, edges, 200, false, distance_thresh);
			if (return_check(ret, "unwrap.MCF()", error_head)) return -1;
		}

	}
	return 0;
}

int Utils::unwrap_3D_adaptive_tiling(
	const Mat& mask,
	const vector<Mat>& quality_map,
	vector<Mat>& wrapped_phase_series,
	vector<Mat>& unwrapped_phase_series,
	const char* delaunay_exe_path,
	const char* mcf_exe_path,
	const char* tmp_file_path,
	double distance_thresh,
	double quality_thresh
)
{
	if (mask.rows < 2 ||
		mask.cols < 2 ||
		mask.type() != CV_32S ||
		wrapped_phase_series.size() < 2 ||
		delaunay_exe_path == NULL ||
		tmp_file_path == NULL ||
		mcf_exe_path == NULL ||
		quality_map.size() != wrapped_phase_series.size() ||
		distance_thresh < 1.0 ||
		quality_thresh > 1.0
		)
	{
		fprintf(stderr, "unwrap_3D_adaptive_tiling(): input check failed!\n");
		return -1;
	}

	/*----------------------------------*/
	/*           自适应分块             */
	/*----------------------------------*/

	/*
	* 确定搜索范围半径
	*/
	int r = (int)floor(distance_thresh), c, ret;
	Mat c_mat = Mat::zeros(2 * r + 1, 1, CV_32S);
	for (int i = 0; i < 2 * r + 1; i++)
	{
		c_mat.at<int>(i, 0) = (int)floor(sqrt(distance_thresh * distance_thresh - double(r - i) * double(r - i)));
	}

	Mat mask1, mask2; mask.copyTo(mask1);
	int nr = mask.rows; int nc = mask.cols;
	int n_images = wrapped_phase_series.size();
	size_t num_nodes = cv::countNonZero(mask);
	if (num_nodes < 3)
	{
		fprintf(stderr, "unwrap_3D_adaptive_tiling(): at least 3 nodes are needed!\n");
		return -1;
	}
	size_t node_count = 0;
	int block_count = 1, row, col, mask2_start_row = nr, mask2_end_row = -1, mask2_start_col = nc, mask2_end_col = -1;
	bool b_break;
	char str[4096];
	node_index node_ix, node_ix2;
	queue<node_index> que;
	vector<Mat> wrapped_phase, unwrapped_phase, qualitymap;
	wrapped_phase.resize(n_images);
	qualitymap.resize(n_images);
	unwrapped_phase_series.resize(n_images);
	for (int i = 0; i < n_images; i++)
	{
		wrapped_phase_series[i].copyTo(unwrapped_phase_series[i]);
	}
	while (node_count != num_nodes)
	{
		mask2_start_row = nr, mask2_end_row = -1, mask2_start_col = nc, mask2_end_col = -1;
		b_break = false;
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if (mask1.at<int>(i, j) == 1)
				{
					mask2_start_row = mask2_start_row > i ? i : mask2_start_row;
					mask2_end_row = mask2_end_row < i ? i : mask2_end_row;
					mask2_start_col = mask2_start_col > j ? j : mask2_start_col;
					mask2_end_col = mask2_end_col < j ? j : mask2_end_col;

					mask1.at<int>(i, j) += block_count; //标注已选点并将其加入到待处理队列中
					node_ix.row = i; node_ix.col = j;
					que.push(node_ix);
					b_break = true;
					break;
				}
			}
			if (b_break) break;
		}
		
		/*
		* 搜索满足条件的点，并加入队列
		*/

		while (!que.empty())
		{
			node_ix = que.front();
			que.pop();
			node_count++;

			for (int i = 0; i <= 2 * r; i++)
			{
				c = c_mat.at<int>(i, 0);
				for (int j = 0; j <= 2 * c; j++)
				{
					row = node_ix.row + (i - r);
					row = row < 0 ? 0 : row; row = row > nr - 1 ? nr - 1 : row;
					col = node_ix.col + (j - c);
					col = col < 0 ? 0 : col; col = col > nc - 1 ? nc - 1 : col;
					if (mask1.at<int>(row, col) == 1)
					{
						mask2_start_row = mask2_start_row > row ? row : mask2_start_row;
						mask2_end_row = mask2_end_row < row ? row : mask2_end_row;
						mask2_start_col = mask2_start_col > col ? col : mask2_start_col;
						mask2_end_col = mask2_end_col < col ? col : mask2_end_col;

						mask1.at<int>(row, col) += block_count; //标注已选点
						node_ix2.row = row; node_ix2.col = col;
						que.push(node_ix2);
					}
				}
			}
		}

		/*
		* 确定新的掩膜矩阵mask2
		*/

		mask2 = Mat::zeros(mask2_end_row - mask2_start_row + 1, mask2_end_col - mask2_start_col + 1, CV_32S);

		for (int i = mask2_start_row; i <= mask2_end_row; i++)
		{
			for (int j = mask2_start_col; j <= mask2_end_col; j++)
			{
				if (mask1.at<int>(i, j) == block_count + 1)
				{
					mask2.at<int>(i - mask2_start_row, j - mask2_start_col) = 1;
				}
			}
		}
		mask2.convertTo(mask2, CV_64F);
		cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\mask2.bin", mask2);
		mask2.convertTo(mask2, CV_32S);
		if (3 > cv::countNonZero(mask2))
		{
			block_count++;
			continue;
		}

		/*
		* 3D相位解缠
		*/
		
		for (int i = 0; i < n_images; i++)
		{
			wrapped_phase_series[i](cv::Range(mask2_start_row, mask2_end_row + 1), cv::Range(mask2_start_col, mask2_end_col + 1)).copyTo
			(wrapped_phase[i]);
			quality_map[i](cv::Range(mask2_start_row, mask2_end_row + 1), cv::Range(mask2_start_col, mask2_end_col + 1)).copyTo
			(qualitymap[i]);
		}

		//ret = unwrap_3D(mask2, qualitymap, wrapped_phase, unwrapped_phase, delaunay_exe_path, tmp_file_path,
		//	distance_thresh + 1.0, quality_thresh);
		ret = unwrap_3D_mcf(mask2, qualitymap, wrapped_phase, unwrapped_phase, delaunay_exe_path,
			mcf_exe_path, tmp_file_path, distance_thresh);
		if (return_check(ret, "unwrap_3D_mcf()", error_head)) return -1;

		for (int i = 0; i < n_images; i++)
		{
			unwrapped_phase[i].copyTo(unwrapped_phase_series[i](cv::Range(mask2_start_row, mask2_end_row + 1),
				cv::Range(mask2_start_col, mask2_end_col + 1)));

			memset(str, 0, 4096);
			sprintf(str, "H:\\data\\experiment\\test\\Regis2\\unwrapped_phase_%d.jpg", i + 1);
			//cvmat2bin(str, unwrapped_phase[i]);
			savephase(str, "jet", unwrapped_phase[i]);
		}
		block_count++;

	}
	mask1.convertTo(mask1, CV_64F);
	cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\mask1.bin", mask1);
	return 0;
}

















tri_node::tri_node()
{
	this->rows = 0;
	this->cols = 0;
	this->num_neigh_edges = 0;
	this->phase = 0;
	this->b_unwrapped = false;
	this->b_balanced = true;
	this->neigh_edges = NULL;
	this->b_residue = false;
	this->epsilon_height = 0.0;
	this->vel = 0.0;
	//std::cout << "constructor1" << "\n";
}

tri_node::tri_node(const tri_node& node)
{
	this->b_unwrapped = node.b_unwrapped;
	this->b_balanced = node.b_balanced;
	this->cols = node.cols;
	this->b_residue = node.b_residue;
	int num_node = node.num_neigh_edges <= 0 ? 1 : node.num_neigh_edges;
	if (node.neigh_edges == NULL)
	{
		this->neigh_edges = NULL;
	}
	else
	{
		this->neigh_edges = (long*)malloc(sizeof(long) * num_node);
		long* ptr = NULL;
		int num;
		node.get_neigh_ptr(&ptr, &num);
		if (this->neigh_edges != NULL || num > 0 || ptr != NULL)
		{
			std::memcpy(this->neigh_edges, ptr, sizeof(long) * num_node);
		}
	}
	
	this->num_neigh_edges = node.num_neigh_edges;
	this->phase = node.phase;
	this->rows = node.rows;
	this->epsilon_height = node.epsilon_height;
	this->vel = node.vel;
	//std::cout << "constructor2" << "\n";
}

tri_node::tri_node(int row, int col, int num_neigh_edge, double phi)
{
	this->rows = row;
	this->cols = col;
	this->num_neigh_edges = num_neigh_edge;
	this->phase = phi;
	this->b_unwrapped = false;
	this->b_residue = false;
	this->b_balanced = true;
	this->epsilon_height = 0.0;
	this->vel = 0.0;
	if (num_neigh_edge > 0)
	{
		this->neigh_edges = (long*)malloc(sizeof(long) * num_neigh_edge);
	}
	else
	{
		this->neigh_edges = NULL;
	}
	if (this->neigh_edges != NULL)
	{
		for (int i = 0; i < num_neigh_edge; i++)
		{
			*(this->neigh_edges + i) = -1;//初始化邻接边序号都为-1
		}
	}
	
	//std::cout << "constructor3" << "\n";
}

tri_node::~tri_node()
{
	if (this->neigh_edges != NULL)
	{
		free(this->neigh_edges);
		this->neigh_edges = NULL;
	}
	//std::cout << "destructor" << "\n";
}

tri_node tri_node::operator=(const tri_node& src)
{
	if (src.neigh_edges == this->neigh_edges && this->neigh_edges != NULL)//两者相等
	{
		return *this;
	}
	else
	{
		if (this->neigh_edges)
		{
			free(this->neigh_edges);
			this->neigh_edges = NULL;
		}
		if (src.num_neigh_edges > 0)
		{
			this->neigh_edges = (long*)malloc(src.num_neigh_edges * sizeof(long));
			if (this->neigh_edges != NULL && src.neigh_edges != NULL)
			{
				memcpy(this->neigh_edges, src.neigh_edges, src.num_neigh_edges * sizeof(long));
			}
		}
		this->b_balanced = src.b_balanced;
		this->b_residue = src.b_residue;
		this->b_unwrapped = src.b_unwrapped;
		this->cols = src.cols;
		this->rows = src.rows;
		this->num_neigh_edges = src.num_neigh_edges;
		this->phase = src.phase;
		this->epsilon_height = src.epsilon_height;
		this->vel = src.vel;
		return *this;
	}
}

int tri_node::get_phase(double* phi) const
{
	if (phi == NULL)
	{
		fprintf(stderr, "get_phase(): input check failed!\n\n");
		return -1;
	}
	*phi = this->phase;
	return 0;
}

int tri_node::get_pos(int* rows, int* cols) const
{
	if (rows == NULL ||
		cols == NULL)
	{
		fprintf(stderr, "tri_node::get_pos(): input check failed!\n\n");
		return -1;
	}
	*rows = this->rows;
	*cols = this->cols;
	return 0;
}

int tri_node::set_phase(double phi)
{
	this->phase = phi;
	return 0;
}

int tri_node::get_neigh_ptr(long** ptr2ptr, int* num) const
{
	if (ptr2ptr == NULL || num == NULL)
	{
		fprintf(stderr, "get_neigh_ptr(): input check failed!\n\n");
		return -1;
	}
	*ptr2ptr = this->neigh_edges;
	*num = this->num_neigh_edges;
	return 0;
}

int tri_node::set_status(bool b_unwrapped)
{
	this->b_unwrapped = b_unwrapped;
	return 0;
}

int tri_node::set_balance(bool b_balanced)
{
	this->b_balanced = b_balanced;
	return 0;
}

int tri_node::print_neighbour() const
{
	if (this->num_neigh_edges <= 0)
	{
		fprintf(stdout, "no neighbour edges!\n");
		return 0;
	}
	for (int i = 0; i < this->num_neigh_edges; i++)
	{
		fprintf(stdout, "%ld ", *(this->neigh_edges + i));
	}
	fprintf(stdout, "\n");
	return 0;
}

int tri_node::get_num_neigh(int* num_neigh) const
{
	*num_neigh = this->num_neigh_edges;
	return 0;
}

int tri_node::get_distance(tri_node node, double* distance) const
{
	*distance = sqrt(((double)node.rows - (double)this->rows) * ((double)node.rows - (double)this->rows) +
		((double)node.cols - (double)this->cols) * ((double)node.cols - (double)this->cols));
	return 0;
}

bool tri_node::get_status() const
{
	return this->b_unwrapped;
}

bool tri_node::get_balance() const
{
	return this->b_balanced;
}

bool tri_node::is_residue_node() const
{
	return this->b_residue;
}

int tri_node::set_residue(bool b_res)
{
	this->b_residue = b_res;
	return 0;
}

double tri_node::get_vel() const
{
	return this->vel;
}

double tri_node::get_height() const
{
	return this->epsilon_height;
}

int tri_node::set_vel(double vel)
{
	this->vel = vel;
	return 0;
}

int tri_node::set_height(double height)
{
	this->epsilon_height = height;
	return 0;
}






