#include<complex.h>
#include "stdafx.h"
#include"..\include\Utils.h"
#include<direct.h>
#include<SensAPI.h>
#include<urlmon.h>
#include<Windows.h>
#include<tchar.h>
#include <atlconv.h>
#include"gdal_priv.h"
#include"gdal.h"
#include"../include/FormatConversion.h"
#include"Eigen/Dense"

#pragma comment(lib,"URlmon")
#pragma comment(lib, "Sensapi.lib")

#ifdef _DEBUG
#pragma comment(lib,"ComplexMat_d.lib")
#pragma comment(lib, "FormatConversion_d.lib")
#else
#pragma comment(lib,"ComplexMat.lib")
#pragma comment(lib, "FormatConversion.lib")
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

}

Utils::~Utils()
{
}

int Utils::createVandermondeMatrix(Mat& inArray, Mat& vandermondeMatrix, int degree)
{
	if (inArray.cols != 1 || inArray.rows < 1 || degree < 1)
	{
		fprintf(stderr, "createVandermondeMatrix(): input check failed!\n");
		return -1;
	}
	vandermondeMatrix.create(inArray.rows, degree + 1, CV_64F);
	if (inArray.type() != CV_64F) inArray.convertTo(inArray, CV_64F);
	for (int i = 0; i < inArray.rows; i++)
	{
		for (int j = 0; j < degree + 1; j++)
		{
			vandermondeMatrix.at<double>(i, j) = pow(inArray.at<double>(i, 0), (double)j);
		}
	}
	return 0;
}

int Utils::ployFit(Mat& a, Mat& B, Mat& x)
{
	Mat A, b;
	a.copyTo(A);
	B.copyTo(b);
	if (A.rows != b.rows || A.cols > A.rows || A.empty())
	{
		fprintf(stderr, "ployFit(): input check failed!\n");
		return -1;
	}
	if (A.type() != CV_64F) A.convertTo(A, CV_64F);
	if (b.type() != CV_64F) b.convertTo(b, CV_64F);
	Mat A_t;
	cv::transpose(A, A_t);
	A = A_t * A;
	b = A_t * b;
	if (!cv::solve(A, b, x, cv::DECOMP_LU))
	{
		fprintf(stderr, "ployFit(): matrix defficiency!\n");
		return -1;
	}
	return 0;
}

int Utils::polyVal(Mat& coefficient, double x, double* val)
{
	if (!val || coefficient.rows < 1 || coefficient.cols != 1) return -1;
	double sum = 0.0;
	if (coefficient.type() != CV_64F) coefficient.convertTo(coefficient, CV_64F);
	for (int i = 0; i < coefficient.rows; i++)
	{
		sum += coefficient.at<double>(i, 0) * pow(x, (double)i);
	}
	*val = sum;
	return 0;
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
		Slave.GetCols() != Master.GetCols() ||
		Master.type() != Slave.type() ||
		(Master.type() != CV_64F && Master.type() != CV_32F))
	{
		fprintf(stderr, "generate_phase(): input check failed!\n\n");
		return -1;
	}
	int rows = Master.GetRows();
	int cols = Master.GetCols();
	phase.create(rows, cols, CV_64F);
	if (Master.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double real = Master.re.at<double>(i, j) * Slave.re.at<double>(i, j) + Master.im.at<double>(i, j) * Slave.im.at<double>(i, j);
				double imag = Slave.re.at<double>(i, j) * Master.im.at<double>(i, j) - Master.re.at<double>(i, j) * Slave.im.at<double>(i, j);
				phase.at<double>(i, j) = atan2(imag, real);
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double real = Master.re.at<float>(i, j) * Slave.re.at<float>(i, j) + Master.im.at<float>(i, j) * Slave.im.at<float>(i, j);
				double imag = Slave.re.at<float>(i, j) * Master.im.at<float>(i, j) - Master.re.at<float>(i, j) * Slave.im.at<float>(i, j);
				phase.at<double>(i, j) = atan2(imag, real);
			}
		}
	}
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
	if (/*!b_balanced*/true)
	{
		Nodes_num = residue.rows * residue.cols + 1;
		Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1) +
			2 * 2 * residue.cols + 2 * 2 * (residue.rows - 2);
	}
	//else
	//{
	//	Nodes_num = residue.rows * residue.cols;
	//	Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1);
	//}
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
	fprintf(fp, "n %ld %lf\n", node_index, -sum);
	long earth_node_index = node_index;

	/*
	* 写入每个有向弧的费用（流费用）
	*/
	long lower_bound = 0;
	long upper_bound = 5;
	double mean_coherence1, mean_coherence2, mean_coherence3, mean_coherence4;
	fprintf(fp, "c Arc descriptor lines (from, to, minflow, maxflow, cost)\n");
	/*接地节点的有向弧流费用*/
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

int Utils::write_DIMACS(const char* DIMACS_problem_file, const Mat& residue, Mat& mask, const Mat& cost, double thresh)
{
	if (residue.cols < 2 ||
		residue.rows < 2 ||
		cost.cols < 2 ||
		cost.rows < 2 ||
		residue.type() != CV_64F ||
		cost.type() != CV_64F ||
		mask.type() != CV_32S ||
		mask.rows != cost.rows ||
		mask.cols != cost.cols ||
		(cost.rows - residue.rows) != 1 ||
		(cost.cols - residue.cols) != 1 ||
		thresh < 0.0)
	{
		fprintf(stderr, "write_DIMACS(): input check failed!\n\n");
		return -1;
	}
	long nr = residue.rows;
	long nc = residue.cols;
	long positive, negative, total, Arcs_num = 0, Nodes_num, feasible_node_num;
	Mat new_mask; mask.copyTo(new_mask);
	new_mask = 1 - new_mask;
	Mat residue_mask = Mat::zeros(nr, nc, CV_32S);
	//根据输入掩膜数据和残差点数据更新掩膜
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(residue.at<double>(i, j)) > thresh)
			{
				residue_mask.at<int>(i, j) = 1;
				new_mask.at<int>(i, j) = 1;
				new_mask.at<int>(i + 1, j) = 1;
				new_mask.at<int>(i + 1, j + 1) = 1;
				new_mask.at<int>(i, j + 1) = 1;
			}
		}
	}
	mask = 1 - new_mask;
	//根据更新的掩膜计算可行的网络节点和流数量
	for (int i = 0; i < nr + 1; i++)
	{
		for (int j = 0; j < nc + 1; j++)
		{
			if (new_mask.at<int>(i, j) == 1)
			{
				int ii, jj;
				ii = i - 1; ii = ii < 0 ? 0 : ii;
				jj = j; jj = jj > nc - 1 ? nc - 1 : jj;
				residue_mask.at<int>(ii, jj) = 1;

				ii = i - 1; ii = ii < 0 ? 0 : ii;
				jj = j - 1; jj = jj < 0 ? 0 : jj;
				residue_mask.at<int>(ii, jj) = 1;

				ii = i; ii = ii > nr - 1 ? nr - 1 : ii;
				jj = j; jj = jj > nc - 1 ? nc - 1 : jj;
				residue_mask.at<int>(ii, jj) = 1;

				ii = i; ii = ii > nr - 1 ? nr - 1 : ii;
				jj = j - 1; jj = jj < 0 ? 0 : jj;
				residue_mask.at<int>(ii, jj) = 1;
			}
		}
	}
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (residue_mask.at<int>(i, j) == 1)
			{
				int ii, jj;
				ii = i; jj = j - 1;
				if (jj >= 0 && residue_mask.at<int>(ii, jj) == 1) Arcs_num++;
				ii = i; jj = j + 1;
				if (jj < nc && residue_mask.at<int>(ii, jj) == 1) Arcs_num++;
				ii = i - 1; jj = j;
				if (ii >= 0 && residue_mask.at<int>(ii, jj) == 1) Arcs_num++;
				ii = i + 1; jj = j;
				if (ii < nr && residue_mask.at<int>(ii, jj) == 1) Arcs_num++;
			}
		}
	}
	//计算可行节点掩膜边缘点数
	int edge_node_num = 0;
	for (int j = 0; j < nc; j++)
	{
		if (residue_mask.at<int>(0, j) == 1) edge_node_num++;
		if (residue_mask.at<int>(nr - 1, j) == 1) edge_node_num++;
	}
	for (int i = 1; i < nr - 1; i++)
	{
		if (residue_mask.at<int>(i, 0) == 1) edge_node_num++;
		if (residue_mask.at<int>(i, nc - 1) == 1) edge_node_num++;
	}
	long i, j;
	long node_index = 1;
	long node_index2;
	double sum = 0.0;
	//统计正负残差点数
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
	if (/*!b_balanced*/1)
	{
		Nodes_num = residue.rows * residue.cols + 1;
		Arcs_num += 2 * 2 * residue.cols + 2 * 2 * (residue.rows - 2);
		/*Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1) +
			2 * 2 * residue.cols + 2 * 2 * (residue.rows - 2);*/
	}
	else
	{
		Nodes_num = residue.rows * residue.cols;
		Arcs_num = Arcs_num;
		//Arcs_num = 2 * (residue.rows - 1) * residue.cols + 2 * residue.rows * (residue.cols - 1);
	}
	ofstream fout;
	FILE* fp = NULL;
	fopen_s(&fp, DIMACS_problem_file, "wt");
	if (!fp)
	{
		fprintf(stderr, "write_DIMACS(): cant't open file %s\n\n", DIMACS_problem_file);
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
	if (/*!b_balanced*/1)
	{
		fprintf(fp, "n %ld %lf\n", node_index, -sum);
	}


	long earth_node_index = node_index;

	/*
	* 写入每个有向弧的费用（流费用）
	*/
	long lower_bound = 0;
	long upper_bound = 5;
	double mean_cost;
	fprintf(fp, "c Arc descriptor lines (from, to, minflow, maxflow, cost)\n");
	/*接地节点的有向弧流费用*/
	if (/*!b_balanced*/1)
	{
		//top
		for (i = 0; i < nc; i++)
		{
			node_index = i + 1;
			mean_cost = cost.at<double>(0, i);
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_cost,
				earth_node_index, node_index, lower_bound, upper_bound, mean_cost);
		}
		//bottom
		for (i = 0; i < nc; i++)
		{
			node_index = nc * (nr - 1) + i + 1;
			mean_cost = cost.at<double>(nr - 1, i);
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_cost,
				earth_node_index, node_index, lower_bound, upper_bound, mean_cost);
		}
		//left
		for (i = 1; i < nr - 1; i++)
		{
			node_index = nc * i + 1;
			mean_cost = cost.at<double>(i, 0);
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_cost,
				earth_node_index, node_index, lower_bound, upper_bound, mean_cost);
		}
		//right
		for (i = 1; i < nr - 1; i++)
		{
			node_index = nc * (i + 1);
			mean_cost = cost.at<double>(i, nc - 1);
			fprintf(fp, "a %ld %ld %ld %ld %lf\na %ld %ld %ld %ld %lf\n",
				node_index, earth_node_index, lower_bound, upper_bound, mean_cost,
				earth_node_index, node_index, lower_bound, upper_bound, mean_cost);
		}
	}

	/*非接地节点的有向弧流费用*/
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++)
		{
			if (residue_mask.at<int>(i, j) == 1)
			{
				node_index = i * nc + j + 1;
				int ii, jj;
				
				ii = i; jj = j - 1;
				if (jj >= 0 && residue_mask.at<int>(ii, jj) == 1)
				{
					node_index2 = ii * nc + jj + 1;
					mean_cost = mean(cost(Range(i, i + 1), Range(j, j + 1))).val[0];
					fprintf(fp, "a %ld %ld %ld %ld %lf\n",
						node_index, node_index2, lower_bound, upper_bound, mean_cost);
				}

				ii = i; jj = j + 1;
				if (jj < nc && residue_mask.at<int>(ii, jj) == 1)
				{
					node_index2 = ii * nc + jj + 1;
					mean_cost = mean(cost(Range(i, i + 1), Range(j, j + 1))).val[0];
					fprintf(fp, "a %ld %ld %ld %ld %lf\n",
						node_index, node_index2, lower_bound, upper_bound, mean_cost);
				}
				ii = i - 1; jj = j;
				if (ii >= 0 && residue_mask.at<int>(ii, jj) == 1)
				{
					node_index2 = ii * nc + jj + 1;
					mean_cost = mean(cost(Range(i, i + 1), Range(j, j + 1))).val[0];
					fprintf(fp, "a %ld %ld %ld %ld %lf\n",
						node_index, node_index2, lower_bound, upper_bound, mean_cost);
				}
				ii = i + 1; jj = j;
				if (ii < nr && residue_mask.at<int>(ii, jj) == 1)
				{
					node_index2 = ii * nc + jj + 1;
					mean_cost = mean(cost(Range(i, i + 1), Range(j, j + 1))).val[0];
					fprintf(fp, "a %ld %ld %ld %ld %lf\n",
						node_index, node_index2, lower_bound, upper_bound, mean_cost);
				}
			}
		}
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
		vec1.type() != vec2.type() ||
		vec1.channels() != 1 ||
		vec1.cols != vec2.cols ||
		vec1.rows != vec2.rows ||
		vec2.channels() != 1 ||
		(vec2.type() != CV_64F && vec2.type() != CV_32F)
		)
	{
		fprintf(stderr, "cross(): input check failed!\n\n");
		return -1;
	}
	if (vec2.type() == CV_64F)
	{
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
	}
	else
	{
		Mat out_tmp = Mat::zeros(vec1.rows, vec1.cols, CV_32F);
		int rows = vec1.rows;
		for (int i = 0; i < rows; i++)
		{
			out_tmp.at<float>(i, 0) = vec1.at<float>(i, 1) * vec2.at<float>(i, 2) -
				vec1.at<float>(i, 2) * vec2.at<float>(i, 1);//a(2) * b(3) - a(3) * b(2)

			out_tmp.at<float>(i, 1) = vec1.at<float>(i, 2) * vec2.at<float>(i, 0) -
				vec1.at<float>(i, 0) * vec2.at<float>(i, 2);//a(3) * b(1) - a(1) * b(3)

			out_tmp.at<float>(i, 2) = vec1.at<float>(i, 0) * vec2.at<float>(i, 1) -
				vec1.at<float>(i, 1) * vec2.at<float>(i, 0);//a(1) * b(2) - a(2) * b(1)
		}
		out_tmp.copyTo(out);
	}
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
	if (rows < 1 || cols < 1 || (Src.type() != CV_64F && Src.type() != CV_32F))
	{
		fprintf(stderr, "wrap(): input check failed!\n\n");
		return -1;
	}
	Mat tmp = Mat::zeros(rows, cols, CV_64F);
	if (Src.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				tmp.at<double>(i, j) = atan2(sin(Src.at<double>(i, j)), cos(Src.at<double>(i, j)));
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				tmp.at<double>(i, j) = atan2(sin(Src.at<float>(i, j)), cos(Src.at<float>(i, j)));
			}
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

int Utils::gen_mask_pdv(Mat& phase_derivatives_variance, Mat& mask, int wndsize, double thresh)
{
	if (phase_derivatives_variance.rows < 2 ||
		phase_derivatives_variance.cols < 2 ||
		phase_derivatives_variance.channels() != 1 ||
		phase_derivatives_variance.type() != CV_64F ||
		wndsize < 0 ||
		wndsize > phase_derivatives_variance.rows ||
		thresh < 0.0 ||
		thresh > 1.0
		)
	{
		fprintf(stderr, "gen_mask(): input check failed!\n\n");
		return -1;
	}
	int nr = phase_derivatives_variance.rows;
	int nc = phase_derivatives_variance.cols;
	Mat temp_coh;
	phase_derivatives_variance.copyTo(temp_coh);
	int radius = (wndsize + 1) / 2;
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
			if (mean < thresh && phase_derivatives_variance.at<double>(i, j) < thresh) mask.at<int>(i, j) = 1;
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
		master_image.type() != slave_image.type() ||
		(slave_image.type() != CV_64F && slave_image.type() != CV_32F) ||
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
	if (master_image.type() == CV_64F)
	{
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
	}
	else
	{
		Mat Coherence(na_new, nr_new, CV_32F, Scalar::all(0));
#pragma omp parallel for schedule(guided)
		for (int i = win_a + 1; i <= na - win_a; i++)
		{
			for (int j = win_r + 1; j <= nr - win_r; j++)
			{
				Mat planes_master[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F) };
				Mat planes_slave[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F) };
				Mat planes[] = { Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F), Mat::zeros(2 * win_a + 1, 2 * win_r + 1, CV_32F) };
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
				Coherence.at<float>(i - 1 - win_a, j - 1 - win_r) = up / (down + 0.0000001);
			}
		}
		copyMakeBorder(Coherence, Coherence, win_a, win_a, win_r, win_r, BORDER_REFLECT);
		Coherence.copyTo(coherence);
	}
	
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
		(master.type() != CV_64F && master.type() != CV_32F) ||
		slave.type() != master.type() ||
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
	tmp.re = master.re.mul(slave.re) + master.im.mul(slave.im);
	tmp.im = slave.re.mul(master.im) - master.re.mul(slave.im);
	int nr = tmp.GetRows();
	int nc = tmp.GetCols();
	int radius_rg = multilook_rg / 2;
	int radius_az = multilook_az / 2;
	phase.create(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		int left, right, bottom, top; double real, imag;
		for (int j = 0; j < nc; j++)
		{
			left = j - radius_rg; left = left < 0 ? 0 : left;
			right = left + multilook_rg; right = right > nc - 1 ? nc - 1 : right;
			top = i - radius_az; top = top < 0 ? 0 : top;
			bottom = top + multilook_az; bottom = bottom > nr - 1 ? nr - 1 : bottom;
			real = cv::mean(tmp.re(Range(top, bottom + 1),Range(left, right + 1)))[0];
			imag = cv::mean(tmp.im(Range(top, bottom + 1), Range(left, right + 1)))[0];
			phase.at<double>(i, j) = atan2(imag, real);
		}
	}
	return 0;
}

int Utils::Multilook(
	const ComplexMat& master, 
	const ComplexMat& slave,
	int multilook_rg, 
	int multilook_az,
	Mat& phase
)
{
	if (master.GetRows() != slave.GetRows() ||
		master.GetCols() != slave.GetCols() ||
		(master.type() != CV_64F && master.type() != CV_32F) ||
		(slave.type() != CV_64F && slave.type() != CV_32F) ||
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
	tmp.re = master.re.mul(slave.re) + master.im.mul(slave.im);
	tmp.im = slave.re.mul(master.im) - master.re.mul(slave.im);
	//ret = master.Mul(slave, tmp, true);
	//if (return_check(ret, "Master.Mul(*, *, *)", error_head)) return -1;
	int nr = tmp.GetRows();
	int nc = tmp.GetCols();
	int nr_new = nr / multilook_az;
	int nc_new = nc / multilook_rg;

	//Mat real = Mat::zeros(nr_new, nc_new, CV_32F);
	//Mat imag = Mat::zeros(nr_new, nc_new, CV_32F);
	phase.create(nr_new, nc_new, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr_new; i++)
	{
		int left, right, bottom, top; double real, imag;
		top = i * multilook_az; top = top < 0 ? 0 : top;
		bottom = top + multilook_az; bottom = bottom > nr ? nr : bottom;
		for (int j = 0; j < nc_new; j++)
		{
			left = j * multilook_rg; left = left < 0 ? 0 : left;
			right = left + multilook_rg; right = right > nc ? nc : right;
			real = cv::mean(tmp.re(Range(top, bottom), Range(left, right)))[0];
			imag = cv::mean(tmp.im(Range(top, bottom), Range(left, right)))[0];
			phase.at<double>(i, j) = atan2(imag, real);
		}
	}
	//tmp.SetRe(real);
	//tmp.SetIm(imag);
	//tmp.GetPhase().copyTo(phase);
	return 0;
}

int Utils::multilook(const Mat& phase, Mat& outPhase, int multi_rg, int multi_az)
{
	if (multi_rg <= 1 && multi_az <= 1)
	{
		phase.copyTo(outPhase);
		return 0;
	}
	ComplexMat slc;
	int ret;
	phase.copyTo(outPhase);
	outPhase.convertTo(outPhase, CV_32F);//节省内存
	ret = phase2cos(outPhase, slc.re, slc.im);
	if (return_check(ret, "phase2cos()", error_head)) return -1;
	int nr = slc.GetRows();
	int nc = slc.GetCols();
	int nr_new = nr / multi_az;
	int nc_new = nc / multi_rg;
	outPhase.create(nr_new, nc_new, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr_new; i++)
	{
		int left, right, bottom, top; double real, imag;
		top = i * multi_az; top = top < 0 ? 0 : top;
		bottom = top + multi_az; bottom = bottom > nr ? nr : bottom;
		for (int j = 0; j < nc_new; j++)
		{
			left = j * multi_rg; left = left < 0 ? 0 : left;
			right = left + multi_rg; right = right > nc ? nc : right;
			real = cv::mean(slc.re(Range(top, bottom), Range(left, right)))[0];
			imag = cv::mean(slc.im(Range(top, bottom), Range(left, right)))[0];
			outPhase.at<double>(i, j) = atan2(imag, real);
		}
	}
	return 0;
}

int Utils::multilook_SAR(const Mat& amplitude, Mat& outAmplitude, int multilook_rg, int multilook_az)
{
	if (amplitude.empty() ||
		amplitude.rows < multilook_rg ||
		amplitude.cols < multilook_az ||
		multilook_rg < 1 || multilook_az < 1||
		(amplitude.type() != CV_64F && amplitude.type() != CV_32F)
		)
	{
		fprintf(stderr, "multilook_SAR(): input check failed!\n");
		return -1;
	}
	int nr = amplitude.rows;
	int nc = amplitude.cols;
	int nr_new = (int)((double)nr / (double)multilook_az);
	int nc_new = (int)((double)nc / (double)multilook_rg);
	Mat tmp;
	if (amplitude.type() == CV_64F)
	{
		tmp = Mat::zeros(nr_new, nc_new, CV_64F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr_new; i++)
		{
			int left, right, bottom, top;
			top = i * multilook_az; top = top < 0 ? 0 : top;
			bottom = top + multilook_az; bottom = bottom > nr ? nr : bottom;
			for (int j = 0; j < nc_new; j++)
			{
				left = j * multilook_rg; left = left < 0 ? 0 : left;
				right = left + multilook_rg; right = right > nc ? nc : right;
				tmp.at<double>(i, j) = cv::mean(amplitude(Range(top, bottom), Range(left, right)))[0];
			}
		}
	}
	else
	{
		tmp = Mat::zeros(nr_new, nc_new, CV_32F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr_new; i++)
		{
			int left, right, bottom, top;
			top = i * multilook_az; top = top < 0 ? 0 : top;
			bottom = top + multilook_az; bottom = bottom > nr ? nr : bottom;
			for (int j = 0; j < nc_new; j++)
			{
				left = j * multilook_rg; left = left < 0 ? 0 : left;
				right = left + multilook_rg; right = right > nc ? nc : right;
				tmp.at<float>(i, j) = cv::mean(amplitude(Range(top, bottom), Range(left, right)))[0];
			}
		}
	}

	tmp.copyTo(outAmplitude);
	return 0;
}

int Utils::phase2cos(const Mat& phase, Mat& cos, Mat& sin)
{
	if (phase.rows < 1 ||
		phase.cols < 1 ||
		(phase.type() != CV_64F && phase.type() != CV_32F) ||
		phase.channels() != 1)
	{
		fprintf(stderr, "phase2cos(): input check failed!\n\n");
		return -1;
	}
	int nr = phase.rows;
	int nc = phase.cols;
	if (phase.type() == CV_32F)
	{
		cos.create(nr, nc, CV_32F);
		sin.create(nr, nc, CV_32F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				cos.at<float>(i, j) = std::cos(phase.at<float>(i, j));
				sin.at<float>(i, j) = std::sin(phase.at<float>(i, j));
			}
		}
	}
	else
	{
		cos.create(nr, nc, CV_64F);
		sin.create(nr, nc, CV_64F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				cos.at<double>(i, j) = std::cos(phase.at<double>(i, j));
				sin.at<double>(i, j) = std::sin(phase.at<double>(i, j));
			}
		}
	}
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

	const double epsilon = 0.000000000000001;
	const double pi = 3.14159265358979323846;
	const double d2r = pi / 180;
	const double r2d = 180 / pi;

	const double a = 6378137.0;		//椭球长半轴
	const double f_inverse = 298.257223563;			//扁率倒数
	const double b = a - a / f_inverse;
	//const double b = 6356752.314245;			//椭球短半轴

	const double e = sqrt(a * a - b * b) / a;

	double x = xyz.at<double>(0, 0);
	double y = xyz.at<double>(0, 1);
	double z = xyz.at<double>(0, 2);

	double tmpX = x;
	double temY = y;
	double temZ = z;

	double curB = 0;
	double N = 0;
	double calB = atan2(temZ, sqrt(tmpX * tmpX + temY * temY));

	int counter = 0;
	while (abs(curB - calB) * r2d > epsilon && counter < 25)
	{
		curB = calB;
		N = a / sqrt(1 - e * e * sin(curB) * sin(curB));
		calB = atan2(temZ + N * e * e * sin(curB), sqrt(tmpX * tmpX + temY * temY));
		counter++;
	}

	x = atan2(temY, tmpX) * r2d;
	y = curB * r2d;
	z = temZ / sin(curB) - N * (1 - e * e);

	llh = Mat::zeros(1, 3, CV_64F);
	llh.at<double>(0, 0) = y;
	llh.at<double>(0, 1) = x;
	llh.at<double>(0, 2) = z;
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
	const double epsilon = 0.000000000000001;
	const double pi = 3.14159265358979323846;
	const double d2r = pi / 180;
	const double r2d = 180 / pi;

	const double a = 6378137.0;		//椭球长半轴
	const double f_inverse = 298.257223563;			//扁率倒数
	const double b = a - a / f_inverse;
	const double e = sqrt(a * a - b * b) / a;

	double y = llh.at<double>(0, 0);
	double x = llh.at<double>(0, 1);
	double z = llh.at<double>(0, 2);

	double L = x * d2r;
	double B = y * d2r;
	double H = z;
	double N = a / sqrt(1 - e * e * sin(B) * sin(B));
	x = (N + H) * cos(B) * cos(L);
	y = (N + H) * cos(B) * sin(L);
	z = (N * (1 - e * e) + H) * sin(B);
	Mat tmp;
	tmp.create(1, 3, CV_64F);
	tmp.at<double>(0, 0) = x;
	tmp.at<double>(0, 1) = y;
	tmp.at<double>(0, 2) = z;
	tmp.copyTo(xyz);
	return 0;
}

int Utils::ell2xyz(double lon, double lat, double elevation, Position& xyz)
{
	if (fabs(lon) > 180.0 || fabs(lat) > 90.0)
	{
		fprintf(stderr, "ell2xyz(): input check failed!\n");
		return -1;
	}
	const double epsilon = 0.000000000000001;
	const double pi = 3.14159265358979323846;
	const double d2r = pi / 180;
	const double r2d = 180 / pi;
	const double a = 6378137.0;		//椭球长半轴
	const double f_inverse = 298.257223563;			//扁率倒数
	const double b = a - a / f_inverse;
	const double e = sqrt(a * a - b * b) / a;
	double y = lat;
	double x = lon;
	double z = elevation;
	double L = x * d2r;
	double B = y * d2r;
	double H = z;
	double N = a / sqrt(1 - e * e * sin(B) * sin(B));
	x = (N + H) * cos(B) * cos(L);
	y = (N + H) * cos(B) * sin(L);
	z = (N * (1 - e * e) + H) * sin(B);
	xyz.x = x;
	xyz.y = y;
	xyz.z = z;
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
	size_t imagesize = nr * nc;
	if (imagesize > 10000 * 10000)
	{
		int mul_times = (int)sqrt(double(imagesize / 10000.0 / 10000.0));
		mul_times = mul_times < 1 ? 1 : mul_times;
		multilook_SAR(mod, mod, mul_times, mul_times);
	}
	double max, min, std;
	this->std(mod, &std);
	nr = mod.rows;
	nc = mod.cols;
	double mean = cv::mean(mod)[0];
	if (SLC.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if ((mod.at<double>(i, j) - mean) >= 2.0 * std) mod.at<double>(i, j) = mean + 2.0 * std;
				if ((mod.at<double>(i, j) - mean) < -2.0 * std) mod.at<double>(i, j) = mean - 2.0 * std;
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
				if ((mod.at<float>(i, j) - mean) >= 2.0 * std) mod.at<float>(i, j) = mean + 2.0 * std;
				if ((mod.at<float>(i, j) - mean) < -2.0 * std) mod.at<float>(i, j) = mean - 2.0 * std;
			}
		}
	}
	cv::minMaxLoc(mod, &min, &max);
	if (fabs(max - min) < 0.00000001)
	{
		fprintf(stderr, "SLC image intensity is the same for every pixel\n\n");
		return -1;
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

int Utils::SAR_image_quantify(const char* filename, double db, ComplexMat& SLC)
{
	if (filename == NULL ||
		db < 0 ||
		SLC.GetRows() < 1 ||
		SLC.GetCols() < 1 /*||
		SLC.type() != CV_64F*/)
	{
		fprintf(stderr, "SAR_image_quantify(): input check failed!\n\n");
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
		min = min < -1.0 ? -1.0 : min;
		max = min + db;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if (mod.at<double>(i, j) >= max) mod.at<double>(i, j) = max;
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
		min = min < -1.0 ? -1.0 : min;
		max = min + db;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if (mod.at<float>(i, j) >= max) mod.at<float>(i, j) = max;
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

int Utils::saveAmplitude(const char* filename, Mat& amplitude)
{
	if (!filename || amplitude.empty() || (amplitude.type() != CV_64F && amplitude.type() != CV_32F))
	{
		fprintf(stderr, "saveAmplitude(): input check failed!\n");
		return -1;
	}
	double min, max;
	int nr = amplitude.rows;
	int nc = amplitude.cols;
	
	double mean = cv::mean(amplitude)[0];
	double std;
	this->std(amplitude, &std);
	if (amplitude.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				if ((amplitude.at<double>(i, j) - mean) >= 3.0 * std) amplitude.at<double>(i, j) = mean + 3.0 * std;
				if ((amplitude.at<double>(i, j) - mean) <= -3.0 * std) amplitude.at<double>(i, j) = mean - 3.0 * std;
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
				if ((amplitude.at<float>(i, j) - mean) >= 3.0 * std) amplitude.at<float>(i, j) = mean + 3.0 * std;
				if ((amplitude.at<float>(i, j) - mean) <= -3.0 * std) amplitude.at<float>(i, j) = mean - 3.0 * std;
			}
		}
	}
	cv::minMaxLoc(amplitude, &min, &max);
	if (fabs(max - min) < 0.00000001)
	{
		fprintf(stderr, "SLC image intensity is the same for every pixel\n\n");
		return -1;
	}
	amplitude = (amplitude - min) / (max - min) * 255.0;
	amplitude.convertTo(amplitude, CV_8U);
	bool ret = cv::imwrite(filename, amplitude);
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

int Utils::save_coherence(const char* filename, const char* colormap, Mat coherence)
{
	if (filename == NULL ||
		colormap == NULL ||
		coherence.rows < 1 ||
		coherence.cols < 1 ||
		(coherence.type() != CV_32F && coherence.type() != CV_64F) ||
		coherence.channels() != 1)
	{
		fprintf(stderr, "save_coherence(): input check failed!\n\n");
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
	coherence.copyTo(tmp);
	if (coherence.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < tmp.rows; i++)
		{
			for (int j = 0; j < tmp.cols; j++)
			{
				if (tmp.at<double>(i, j) < 0.0) tmp.at<double>(i, j) = 0.0;
				if (tmp.at<double>(i, j) > 1.0) tmp.at<double>(i, j) = 1.0;
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < tmp.rows; i++)
		{
			for (int j = 0; j < tmp.cols; j++)
			{
				if (tmp.at<float>(i, j) < 0.0) tmp.at<float>(i, j) = 0.0;
				if (tmp.at<float>(i, j) > 1.0) tmp.at<float>(i, j) = 1.0;
			}
		}
	}
	min = 0.0; max = 1.0;
	tmp = (tmp - min) / (max - min) * 255.0;
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

int Utils::savephase_black(const char* filename, const char* colormap, Mat& phase, Mat& mask)
{
	if (filename == NULL ||
		colormap == NULL ||
		phase.rows < 1 ||
		phase.cols < 1 ||
		phase.type() != CV_64F ||
		phase.channels() != 1 ||
		mask.size() != phase.size() ||
		mask.type() != CV_32S
		)
	{
		fprintf(stderr, "savephase_black(): input check failed!\n\n");
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
		fprintf(stderr, "savephase_black(): phase value is the same for every pixel\n\n");
		return -1;
	}
	tmp = (tmp - min) / (max - min) * 255.0;
	tmp.convertTo(tmp, CV_8UC3);
	if (!gray)
	{
		cv::applyColorMap(tmp, tmp, type);
	}
	for (int i = 0; i < tmp.rows; i++)
	{
		for (int j = 0; j < tmp.cols; j++)
		{
			if (mask.at<int>(i, j) == 0)
			{
				tmp.at<Vec<uchar, 3>>(i, j)[0] = 0;
				tmp.at<Vec<uchar, 3>>(i, j)[1] = 0;
				tmp.at<Vec<uchar, 3>>(i, j)[2] = 0;
			}
		}
	}
	bool ret = cv::imwrite(filename, tmp);
	if (!ret)
	{
		fprintf(stderr, "cv::imwrite(): can't write to %s\n\n", filename);
		return -1;
	}
	return 0;
}

int Utils::savephase_white(const char* filename, const char* colormap, Mat& phase, Mat& mask)
{
	if (filename == NULL ||
		colormap == NULL ||
		phase.rows < 1 ||
		phase.cols < 1 ||
		phase.type() != CV_64F ||
		phase.channels() != 1 ||
		mask.size() != phase.size() ||
		mask.type() != CV_32S
		)
	{
		fprintf(stderr, "savephase_white(): input check failed!\n\n");
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
		fprintf(stderr, "savephase_white() : phase value is the same for every pixel\n\n");
		return -1;
	}
	tmp = (tmp - min) / (max - min) * 255.0;
	tmp.convertTo(tmp, CV_8UC3);
	if (!gray)
	{
		cv::applyColorMap(tmp, tmp, type);
	}
	for (int i = 0; i < tmp.rows; i++)
	{
		for (int j = 0; j < tmp.cols; j++)
		{
			if (mask.at<int>(i, j) == 0)
			{
				tmp.at<Vec<uchar, 3>>(i, j)[0] = 255;
				tmp.at<Vec<uchar, 3>>(i, j)[1] = 255;
				tmp.at<Vec<uchar, 3>>(i, j)[2] = 255;
			}
		}
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
	node_neighbours.resize(num_nodes);
	for (int i = 0; i < num_nodes; i++)
	{
		node_neighbours[i] = 0;
	}
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
		char delaunay_job_name[512]; delaunay_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(delaunay_job_name, "DELAUNAY_%lld", tt);
		string delaunay_job_name_string(delaunay_job_name);
		HANDLE hd = CreateJobObjectA(NULL, delaunay_job_name_string.c_str());
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
		//input.type() != CV_64F ||
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

//int Utils::stack_coregistration(
//	vector<ComplexMat>& SAR_images,
//	Mat& offset,
//	int Master_index,
//	int coh_method,
//	int interp_times,
//	int blocksize
//)
//{
//	if (SAR_images.size() < 2||
//		Master_index < 1 ||
//		Master_index > SAR_images.size() ||
//		coh_method < 0 ||
//		coh_method > 1 ||
//		interp_times < 1||
//		blocksize < 4
//		)
//	{
//		fprintf(stderr, "stack_coregistration(): input check failed!\n\n");
//		return -1;
//	}
//	vector<ComplexMat> SAR_images_out;
//	SAR_images_out.resize(SAR_images.size());
//	Registration coregis;
//	int nr = SAR_images[Master_index - 1].GetRows();
//	int nc = SAR_images[Master_index - 1].GetCols();
//	int n_images = SAR_images.size();
//	Mat move_row = Mat::zeros(1, n_images, CV_32S);
//	Mat move_col = Mat::zeros(1, n_images, CV_32S);
//	Mat Slave_indx = Mat::zeros(1, n_images - 1, CV_32S);
//	Mat offset_topleft = Mat::zeros(n_images, 2, CV_32S);
//	int count = 0;
//	for (int i = 0; i < n_images; i++)
//	{
//		if (i != Master_index - 1)
//		{
//			Slave_indx.at<int>(0, count) = i;
//			count++;
//		}
//	}
//	////////////粗配准///////////////
//	ComplexMat master;
//	master = SAR_images[Master_index - 1];
//	volatile bool parallel_flag = true;
//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < n_images - 1; i++)
//	{
//		if (!parallel_flag) continue;
//		ComplexMat slave;
//		int ret, move_r, move_c;
//		int slave_img = Slave_indx.at<int>(0, i);
//		slave = SAR_images[slave_img];
//		ret = coregis.real_coherent(master, slave, &move_r, &move_c);
//		if (ret < 0)
//		{
//			parallel_flag = false;
//			continue;
//		}
//		move_row.at<int>(0, slave_img) = move_r;
//		move_col.at<int>(0, slave_img) = move_c;
//	}
//	if (parallel_check(parallel_flag, "stack_coregistration()", parallel_error_head)) return -1;
//	////////////粗配准公共部分裁剪///////////////
//	int move_r_min, move_r_max, move_c_min, move_c_max;
//	move_r_min = 100000000;
//	move_r_max = -100000000;
//	move_c_min = 100000000;
//	move_c_max = -100000000;
//
//	int cut_rows;
//	for (int i = 0; i < n_images; i++)
//	{
//		move_r_min = move_r_min < move_row.at<int>(0, i) ? move_r_min : move_row.at<int>(0, i);
//		move_c_min = move_c_min < move_col.at<int>(0, i) ? move_c_min : move_col.at<int>(0, i);
//		move_r_max = move_r_max > move_row.at<int>(0, i) ? move_r_max : move_row.at<int>(0, i);
//		move_c_max = move_c_max > move_col.at<int>(0, i) ? move_c_max : move_col.at<int>(0, i);
//	}
//	if (abs(move_r_max) >= nr ||
//		abs(move_r_min) >= nr ||
//		abs(move_c_min) >= nc ||
//		abs(move_c_max) >= nc ||
//		abs(move_r_max - move_r_min) >= nr||
//		abs(move_c_max - move_c_min) >= nc
//		)
//	{
//		fprintf(stderr, "stack_coregistration(): SAR images have no common area!\n\n");
//		return -1;
//	}
//	if (move_r_min >= 0)
//	{
//		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(0, nr - move_r_max), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = 0;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = move_row.at<int>(0, slave_ix);
//			int row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
//			int col_start = 0;
//			int col_end = SAR_images[slave_ix].GetCols();
//			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix);
//		}
//	}
//	if (move_r_max <= 0 && move_r_min < 0)
//	{
//		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(-move_r_min, nr), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = move_row.at<int>(0, slave_ix) - move_r_min;
//			int row_end = nr + move_row.at<int>(0, slave_ix);
//			int col_start = 0;
//			int col_end = SAR_images[slave_ix].GetCols();
//			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
//		}
//	}
//	if (move_r_max > 0 && move_r_min < 0)
//	{
//		SAR_images_out[Master_index - 1] = SAR_images[Master_index - 1](Range(-move_r_min, nr - move_r_max), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = move_row.at<int>(0, slave_ix) - move_r_min;
//			int row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
//			int col_start = 0;
//			int col_end = SAR_images[slave_ix].GetCols();
//			SAR_images_out[slave_ix] = SAR_images[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
//		}
//	}
//
//
//	if (move_c_min >= 0)
//	{
//		cut_rows = SAR_images_out[Master_index - 1].GetRows();
//		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(0, nc - move_c_max));
//		offset_topleft.at<int>(Master_index - 1, 1) = 0;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = 0;
//			int row_end = SAR_images_out[slave_ix].GetRows();
//			int col_start = move_col.at<int>(0, slave_ix);
//			int col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
//			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix);
//		}
//	}
//	if (move_c_max <= 0 && move_c_min < 0)
//	{
//		cut_rows = SAR_images_out[Master_index - 1].GetRows();
//		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(-move_c_min, nc));
//		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = 0;
//			int row_end = SAR_images_out[slave_ix].GetRows();
//			int col_start = move_col.at<int>(0, slave_ix) - move_c_min;
//			int col_end = nc + move_col.at<int>(0, slave_ix);
//			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
//		}
//	}
//	if (move_c_max > 0 && move_c_min < 0)
//	{
//		cut_rows = SAR_images_out[Master_index - 1].GetRows();
//		SAR_images_out[Master_index - 1] = SAR_images_out[Master_index - 1](Range(0, cut_rows), Range(-move_c_min, nc - move_c_max));
//		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < n_images - 1; i++)
//		{
//			int slave_ix = Slave_indx.at<int>(0, i);
//			int row_start = 0;
//			int row_end = SAR_images_out[slave_ix].GetRows();
//			int col_start = move_col.at<int>(0, slave_ix) - move_c_min;
//			int col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
//			SAR_images_out[slave_ix] = SAR_images_out[slave_ix](Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
//		}
//	}
//	offset_topleft.copyTo(offset);
//
//	////////////////////////精配准////////////////////////////
//	SAR_images[Master_index - 1] = SAR_images_out[Master_index - 1];
//	master = SAR_images_out[Master_index - 1];
//	//这里并行加速，因为子函数已经进行了加速
////#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < n_images - 1; i++)
//	{
//		ComplexMat slave;
//		int slave_ix = Slave_indx.at<int>(0, i);
//		slave = SAR_images_out[slave_ix];
//		int ret;
//		ret = coregis.coregistration_subpixel(master, slave, blocksize, interp_times);
//		if (return_check(ret, "registration_subpixel(*, *, *, *)", error_head)) return -1;
//		SAR_images[slave_ix] = slave;
//		//cout << i + 1 << "/" << n_images - 1 << "\n";
//	}
//	return 0;
//}
//
//int Utils::stack_coregistration(
//	vector<string>& SAR_images,
//	vector<string>& SAR_images_out,
//	Mat& offset,
//	int Master_index,
//	int interp_times,
//	int blocksize
//)
//{
//	if (SAR_images.size() < 2 ||
//		Master_index < 1 ||
//		Master_index > SAR_images.size() ||
//		SAR_images_out.size() != SAR_images.size() ||
//		interp_times < 1 ||
//		blocksize < 16
//		)
//	{
//		fprintf(stderr, "stack_coregistration(): input check failed!\n\n");
//		return -1;
//	}
//
//	Registration coregis; FormatConversion conversion;
//	int n_images = SAR_images.size(), ret;
//	Mat move_row = Mat::zeros(1, n_images, CV_32S);
//	Mat move_col = Mat::zeros(1, n_images, CV_32S);
//	Mat Slave_indx = Mat::zeros(1, n_images - 1, CV_32S);
//	Mat offset_topleft = Mat::zeros(n_images, 2, CV_32S);
//	int count = 0;
//	//创建配准后输出的h5文件
//	for (int i = 0; i < n_images; i++)
//	{
//		ret = conversion.creat_new_h5(SAR_images_out[i].c_str());
//		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
//	}
//	for (int i = 0; i < n_images; i++)
//	{
//		if (i != Master_index - 1)
//		{
//			Slave_indx.at<int>(0, count) = i;
//			count++;
//		}
//	}
//	//检查并确保SAR图像尺寸大小一致
//	int min_row = 10000000, min_col = 10000000;
//	Mat azimuth_len, range_len;
//	for (int i = 0; i < n_images; i++)
//	{
//		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "azimuth_len", azimuth_len);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "range_len", range_len);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		min_row = azimuth_len.at<int>(0, 0) > min_row ? min_row : azimuth_len.at<int>(0, 0);
//		min_col = range_len.at<int>(0, 0) > min_col ? min_col : range_len.at<int>(0, 0);
//	}
//	if (min_col < 1 || min_row < 1)
//	{
//		fprintf(stderr, "invalide SAR images size\n");
//		return -1;
//	}
//	
//	////////////粗配准///////////////
//	ComplexMat master, slave, master_out, slave_out;
//	ret = conversion.read_slc_from_h5(SAR_images[Master_index - 1].c_str(), master);
//	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//	master = master(cv::Range(0, min_row), cv::Range(0, min_col));
//	if (master.type() != CV_64F) master.convertTo(master, CV_64F);
//	int nr = master.GetRows();
//	int nc = master.GetCols();
//	for (int i = 0; i < n_images - 1; i++)
//	{
//		
//		int move_r, move_c;
//		int slave_img = Slave_indx.at<int>(0, i);
//		ret = conversion.read_slc_from_h5(SAR_images[slave_img].c_str(), slave);
//		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//		slave = slave(cv::Range(0, min_row), cv::Range(0, min_col));
//		if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
//		ret = coregis.real_coherent(master, slave, &move_r, &move_c);
//		if (return_check(ret, "real_coherent()", error_head)) return -1;
//		move_row.at<int>(0, slave_img) = move_r;
//		move_col.at<int>(0, slave_img) = move_c;
//	}
//	////////////粗配准公共部分裁剪///////////////
//	int move_r_min, move_r_max, move_c_min, move_c_max;
//	move_r_min = 100000000;
//	move_r_max = -100000000;
//	move_c_min = 100000000;
//	move_c_max = -100000000;
//
//	int cut_rows;
//	for (int i = 0; i < n_images; i++)
//	{
//		move_r_min = move_r_min < move_row.at<int>(0, i) ? move_r_min : move_row.at<int>(0, i);
//		move_c_min = move_c_min < move_col.at<int>(0, i) ? move_c_min : move_col.at<int>(0, i);
//		move_r_max = move_r_max > move_row.at<int>(0, i) ? move_r_max : move_row.at<int>(0, i);
//		move_c_max = move_c_max > move_col.at<int>(0, i) ? move_c_max : move_col.at<int>(0, i);
//	}
//	if (abs(move_r_max) >= nr ||
//		abs(move_r_min) >= nr ||
//		abs(move_c_min) >= nc ||
//		abs(move_c_max) >= nc ||
//		abs(move_r_max - move_r_min) >= nr ||
//		abs(move_c_max - move_c_min) >= nc
//		)
//	{
//		fprintf(stderr, "stack_coregistration(): SAR images have no common area!\n\n");
//		return -1;
//	}
//
//	//主图像裁剪
//	if (move_r_min >= 0)
//	{
//		master = master(Range(0, nr - move_r_max), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = 0;
//	}
//	if (move_r_max <= 0 && move_r_min < 0)
//	{
//		master = master(Range(-move_r_min, nr), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
//	}
//	if (move_r_max > 0 && move_r_min < 0)
//	{
//		master = master(Range(-move_r_min, nr - move_r_max), Range(0, nc));
//		offset_topleft.at<int>(Master_index - 1, 0) = -move_r_min;
//	}
//	if (move_c_min >= 0)
//	{
//		cut_rows = master.GetRows();
//		master = master(Range(0, cut_rows), Range(0, nc - move_c_max));
//		offset_topleft.at<int>(Master_index - 1, 1) = 0;
//	}
//	if (move_c_max <= 0 && move_c_min < 0)
//	{
//		cut_rows = master.GetRows();
//		master = master(Range(0, cut_rows), Range(-move_c_min, nc));
//		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
//	}
//	if (move_c_max > 0 && move_c_min < 0)
//	{
//		cut_rows = master.GetRows();
//		master = master(Range(0, cut_rows), Range(-move_c_min, nc - move_c_max));
//		offset_topleft.at<int>(Master_index - 1, 1) = -move_c_min;
//	}
//	ret = conversion.write_slc_to_h5(SAR_images_out[Master_index - 1].c_str(), master);
//	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
//	azimuth_len.at<int>(0, 0) = master.GetRows();
//	range_len.at<int>(0, 0) = master.GetCols();
//	ret = conversion.write_array_to_h5(SAR_images_out[Master_index - 1].c_str(), "azimuth_len", azimuth_len);
//	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//	ret = conversion.write_array_to_h5(SAR_images_out[Master_index - 1].c_str(), "range_len", range_len);
//	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//	//辅图像裁剪
//	for (int i = 0; i < n_images - 1; i++)
//	{
//		int slave_ix, row_start, row_end, col_start, col_end;
//		slave_ix = Slave_indx.at<int>(0, i);
//		ret = conversion.read_slc_from_h5(SAR_images[slave_ix].c_str(), slave);
//		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//		if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
//		if (move_r_min >= 0)
//		{
//			
//			row_start = move_row.at<int>(0, slave_ix);
//			row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
//			col_start = 0;
//			col_end = slave.GetCols();
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix);
//		}
//		if (move_r_max <= 0 && move_r_min < 0)
//		{
//			row_start = move_row.at<int>(0, slave_ix) - move_r_min;
//			row_end = nr + move_row.at<int>(0, slave_ix);
//			col_start = 0;
//			col_end = slave.GetCols();
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
//		}
//		if (move_r_max > 0 && move_r_min < 0)
//		{
//			row_start = move_row.at<int>(0, slave_ix) - move_r_min;
//			row_end = nr + move_row.at<int>(0, slave_ix) - move_r_max;
//			col_start = 0;
//			col_end = slave.GetCols();
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 0) = move_row.at<int>(0, slave_ix) - move_r_min;
//		}
//
//
//		if (move_c_min >= 0)
//		{
//			row_start = 0;
//			row_end = slave.GetRows();
//			col_start = move_col.at<int>(0, slave_ix);
//			col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix);
//		}
//		if (move_c_max <= 0 && move_c_min < 0)
//		{
//			row_start = 0;
//			row_end = slave.GetRows();
//			col_start = move_col.at<int>(0, slave_ix) - move_c_min;
//			col_end = nc + move_col.at<int>(0, slave_ix);
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
//		}
//		if (move_c_max > 0 && move_c_min < 0)
//		{
//			row_start = 0;
//			row_end = slave.GetRows();
//			col_start = move_col.at<int>(0, slave_ix) - move_c_min;
//			col_end = nc + move_col.at<int>(0, slave_ix) - move_c_max;
//			slave = slave(Range(row_start, row_end), Range(col_start, col_end));
//			offset_topleft.at<int>(slave_ix, 1) = move_col.at<int>(0, slave_ix) - move_c_min;
//		}
//
//		//精配准
//		ret = coregis.coregistration_subpixel(master, slave, blocksize, interp_times);
//		if (return_check(ret, "coregistration_subpixel()", error_head)) return -1;
//		//写出
//		ret = conversion.write_slc_to_h5(SAR_images_out[slave_ix].c_str(), slave);
//		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
//		azimuth_len.at<int>(0, 0) = slave.GetRows();
//		range_len.at<int>(0, 0) = slave.GetCols();
//		ret = conversion.write_array_to_h5(SAR_images_out[slave_ix].c_str(), "azimuth_len", azimuth_len);
//		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		ret = conversion.write_array_to_h5(SAR_images_out[slave_ix].c_str(), "range_len", range_len);
//		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//	}
//	offset_topleft.copyTo(offset);
//
//
//	return 0;
//}
//
//int Utils::stack_coregistration(
//	vector<string>& SAR_images, 
//	vector<string>& SAR_images_out,
//	int Master_index,
//	int interp_times,
//	int blocksize
//)
//{
//	if (SAR_images.size() < 2 ||
//		Master_index < 1 ||
//		Master_index > SAR_images.size() ||
//		SAR_images_out.size() != SAR_images.size() ||
//		interp_times < 1 ||
//		blocksize < 16
//		)
//	{
//		fprintf(stderr, "stack_coregistration(): input check failed!\n\n");
//		return -1;
//	}
//	//获取各图像的尺寸，并创建输出h5文件
//	FormatConversion conversion;
//	int ret, type;
//	int n_images = SAR_images.size();
//	Mat images_rows, images_cols, tmp;
//	images_rows = Mat::zeros(n_images, 1, CV_32S); images_cols = Mat::zeros(n_images, 1, CV_32S);
//	for (int i = 0; i < n_images; i++)
//	{
//		ret = conversion.creat_new_h5(SAR_images_out[i].c_str());
//		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
//
//		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "range_len", tmp);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		images_cols.at<int>(i, 0) = tmp.at<int>(0, 0);
//
//		ret = conversion.read_array_from_h5(SAR_images[i].c_str(), "azimuth_len", tmp);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		images_rows.at<int>(i, 0) = tmp.at<int>(0, 0);
//	}
//	//分块读取数据并求取偏移量
//	Utils util; Registration regis;
//	int rows = images_rows.at<int>(Master_index - 1, 0); int cols = images_cols.at<int>(Master_index - 1, 0);
//	int m = rows / blocksize;
//	int n = cols / blocksize;
//	if (m * n < 10)
//	{
//		fprintf(stderr, "stack_coregistration(): try smaller blocksize!\n");
//		return -1;
//	}
//	Mat offset_r = Mat::zeros(m, n, CV_64F); Mat offset_c = Mat::zeros(m, n, CV_64F);
//	Mat offset_coord_row = Mat::zeros(m, n, CV_64F);
//	Mat offset_coord_col = Mat::zeros(m, n, CV_64F);
//	//子块中心坐标
//	for (int i = 0; i < m; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			offset_coord_row.at<double>(i, j) = ((double)blocksize) / 2 * (double)(2 * i + 1);
//			offset_coord_col.at<double>(i, j) = ((double)blocksize) / 2 * (double)(2 * j + 1);
//		}
//	}
//	//根据输入图像尺寸大小判断是否分块读取（超过20000×20000则分块读取，否则一次性读取）
//	ComplexMat master_w, slave_w;
//	bool b_block = true; bool master_read = false;
//	if (rows * cols < 20000 * 20000) b_block = false;
//	for (int ii = 0; ii < n_images; ii++)
//	{
//		if (ii == Master_index - 1) continue;
//		if (!b_block)//不分块读取
//		{
//			if (!master_read)
//			{
//				ret = conversion.read_slc_from_h5(SAR_images[Master_index - 1].c_str(), master_w);
//				if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//				master_read = true;
//				type = master_w.type();
//				ret = conversion.write_slc_to_h5(SAR_images_out[Master_index - 1].c_str(), master_w);//写主图像
//			}
//			
//			ret = conversion.read_slc_from_h5(SAR_images[ii].c_str(), slave_w);
//			if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//			if (type != slave_w.type())
//			{
//				fprintf(stderr, "stack_coregistration(): images type mismatch!\n");
//				return -1;
//			}
//			if (type != CV_16S && type != CV_64F)
//			{
//				fprintf(stderr, "stack_coregistration(): data type not supported!\n");
//				return -1;
//			}
//		}
//		//分块读取并计算偏移量
//		int mm, nn;
//		mm = images_rows.at<int>(ii, 0) / blocksize;
//		nn = images_cols.at<int>(ii, 0) / blocksize;
//		if (!b_block)
//		{
//# pragma omp parallel for schedule(guided)
//			
//			for (int j = 0; j < m; j++)
//			{
//				int offset_row, offset_col, move_r, move_c;
//				ComplexMat master, slave, master_interp, slave_interp;
//				for (int k = 0; k < n; k++)
//				{
//					offset_row = j * blocksize; offset_col = k * blocksize;
//					if ((j + 1) * blocksize < images_rows.at<int>(ii, 0) && (k + 1) * blocksize < images_cols.at<int>(ii, 0))
//					{
//						master = master_w(cv::Range(offset_row, offset_row + blocksize), cv::Range(offset_col, offset_col + blocksize));
//						slave = slave_w(cv::Range(offset_row, offset_row + blocksize), cv::Range(offset_col, offset_col + blocksize));
//
//
//
//						//计算偏移量
//						if (master.type() != CV_64F) master.convertTo(master, CV_64F);
//						if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
//						move_r = 0; move_c = 0;
//						ret = regis.interp_paddingzero(master, master_interp, interp_times);
//						//if (return_check(ret, "interp_paddingzero()", error_head)) return -1;
//						ret = regis.interp_paddingzero(slave, slave_interp, interp_times);
//						//if (return_check(ret, "interp_paddingzero()", error_head)) return -1;
//						ret = regis.real_coherent(master_interp, slave_interp, &move_r, &move_c);
//						//if (return_check(ret, "real_coherent()", error_head)) return -1;
//						offset_r.at<double>(j, k) = double(move_r) / double(interp_times);
//						offset_c.at<double>(j, k) = double(move_c) / double(interp_times);
//					}
//
//				}
//			}
//			
//		}
//		else
//		{
//			int offset_row, offset_col, move_r, move_c;
//			ComplexMat master, slave, master_interp, slave_interp;
//			for (int j = 0; j < m; j++)
//			{
//				for (int k = 0; k < n; k++)
//				{
//					offset_row = j * blocksize; offset_col = k * blocksize;
//					if ((j + 1) * blocksize < images_rows.at<int>(ii, 0) && (k + 1) * blocksize < images_cols.at<int>(ii, 0))
//					{
//						//mm = j + 1; nn = k + 1;//记录实际的子块行列数
//						ret = conversion.read_subarray_from_h5(SAR_images[Master_index - 1].c_str(), "s_im", offset_row, offset_col, blocksize, blocksize, master.im);
//						if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//						ret = conversion.read_subarray_from_h5(SAR_images[Master_index - 1].c_str(), "s_re", offset_row, offset_col, blocksize, blocksize, master.re);
//						if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//						ret = conversion.read_subarray_from_h5(SAR_images[ii].c_str(), "s_im", offset_row, offset_col, blocksize, blocksize, slave.im);
//						if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//						ret = conversion.read_subarray_from_h5(SAR_images[ii].c_str(), "s_re", offset_row, offset_col, blocksize, blocksize, slave.re);
//						if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//
//						//计算偏移量
//						if (master.type() != CV_64F) master.convertTo(master, CV_64F);
//						if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
//
//						ret = regis.interp_paddingzero(master, master_interp, interp_times);
//						if (return_check(ret, "interp_paddingzero()", error_head)) return -1;
//						ret = regis.interp_paddingzero(slave, slave_interp, interp_times);
//						if (return_check(ret, "interp_paddingzero()", error_head)) return -1;
//						ret = regis.real_coherent(master_interp, slave_interp, &move_r, &move_c);
//						if (return_check(ret, "real_coherent()", error_head)) return -1;
//						offset_r.at<double>(j, k) = double(move_r) / double(interp_times);
//						offset_c.at<double>(j, k) = double(move_c) / double(interp_times);
//					}
//
//				}
//			}
//		}
//		
//
//		//剔除outliers
//		m = mm; n = nn;//更新实际子块行列数
//		Mat sentinel = Mat::zeros(m, n, CV_64F);
//		int ix, iy, count = 0, c = 0; double delta, thresh = 2.0;
//		for (int i = 0; i < m; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				count = 0;
//				//上
//				ix = j;
//				iy = i - 1; iy = iy < 0 ? 0 : iy;
//				delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
//				delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
//				if (fabs(delta) >= thresh) count++;
//				//下
//				ix = j;
//				iy = i + 1; iy = iy > m - 1 ? m - 1 : iy;
//				delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
//				delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
//				if (fabs(delta) >= thresh) count++;
//				//左
//				ix = j - 1; ix = ix < 0 ? 0 : ix;
//				iy = i;
//				delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
//				delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
//				if (fabs(delta) >= thresh) count++;
//				//右
//				ix = j + 1; ix = ix > n - 1 ? n - 1 : ix;
//				iy = i;
//				delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
//				delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
//				if (fabs(delta) >= thresh) count++;
//
//				if (count > 2) { sentinel.at<double>(i, j) = 1.0; c++; }
//			}
//		}
//		Mat offset_c_0, offset_r_0, offset_coord_row_0, offset_coord_col_0;
//		offset_c_0 = Mat::zeros(m * n - c, 1, CV_64F);
//		offset_r_0 = Mat::zeros(m * n - c, 1, CV_64F);
//		offset_coord_row_0 = Mat::zeros(m * n - c, 1, CV_64F);
//		offset_coord_col_0 = Mat::zeros(m * n - c, 1, CV_64F);
//		count = 0;
//		for (int i = 0; i < m; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				if (sentinel.at<double>(i, j) < 0.5)
//				{
//					offset_r_0.at<double>(count, 0) = offset_r.at<double>(i, j);
//					offset_c_0.at<double>(count, 0) = offset_c.at<double>(i, j);
//					offset_coord_row_0.at<double>(count, 0) = offset_coord_row.at<double>(i, j);
//					offset_coord_col_0.at<double>(count, 0) = offset_coord_col.at<double>(i, j);
//					count++;
//				}
//			}
//		}
//
//
//		m = 1; n = count;
//		if (count < 10)
//		{
//			fprintf(stderr, "stack_coregistration(): insufficient valide sub blocks!\n");
//			return -1;
//		}
//		//偏移量拟合（坐标做归一化处理）
//		//拟合公式为 offser_row / offser_col = a0 + a1 * x + a2 * y;
//		double offset_x = (double)cols / 2;
//		double offset_y = (double)rows / 2;
//		double scale_x = (double)cols;
//		double scale_y = (double)rows;
//		offset_coord_row_0 -= offset_y;
//		offset_coord_col_0 -= offset_x;
//		offset_coord_row_0 /= scale_y;
//		offset_coord_col_0 /= scale_x;
//		Mat A = Mat::ones(m * n, 3, CV_64F);
//		Mat temp, A_t;
//		offset_coord_col_0.copyTo(A(Range(0, m * n), Range(1, 2)));
//
//		offset_coord_row_0.copyTo(A(Range(0, m * n), Range(2, 3)));
//
//	
//		cv::transpose(A, A_t);
//
//		Mat b_r, b_c, coef_r, coef_c, error_r, error_c, b_t, a, a_t;
//
//		A.copyTo(a);
//		cv::transpose(a, a_t);
//		offset_r_0.copyTo(b_r);
//		b_r = A_t * b_r;
//
//		offset_c_0.copyTo(b_c);
//		b_c = A_t * b_c;
//
//		A = A_t * A;
//
//		double rms1 = -1.0; double rms2 = -1.0;
//		Mat eye = Mat::zeros(m * n, m * n, CV_64F);
//		for (int i = 0; i < m * n; i++)
//		{
//			eye.at<double>(i, i) = 1.0;
//		}
//		if (cv::invert(A, error_r, cv::DECOMP_LU) > 0)
//		{
//			cv::transpose(offset_r_0, b_t);
//			error_r = b_t * (eye - a * error_r * a_t) * offset_r_0;
//			rms1 = sqrt(error_r.at<double>(0, 0) / double(m * n));
//		}
//		if (cv::invert(A, error_c, cv::DECOMP_LU) > 0)
//		{
//			cv::transpose(offset_c_0, b_t);
//			error_c = b_t * (eye - a * error_c * a_t) * offset_c_0;
//			rms2 = sqrt(error_c.at<double>(0, 0) / double(m * n));
//		}
//		if (!cv::solve(A, b_r, coef_r, cv::DECOMP_NORMAL))
//		{
//			fprintf(stderr, "stack_coregistration(): matrix defficiency!\n");
//			return -1;
//		}
//		if (!cv::solve(A, b_c, coef_c, cv::DECOMP_NORMAL))
//		{
//			fprintf(stderr, "stack_coregistration(): matrix defficiency!\n");
//			return -1;
//		}
//
//		/*---------------------------------------*/
//	    /*    双线性插值获取重采样后的辅图像     */
//	    /*---------------------------------------*/
//		ComplexMat slave1;
//		ret = conversion.read_slc_from_h5(SAR_images[ii].c_str(), slave1);
//		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
//		//if (slave1.type() != CV_16S) slave1.convertTo(slave1, CV_16S);
//		int rows_slave, cols_slave;
//		rows_slave = slave1.GetRows(); cols_slave = slave1.GetCols();
//		type = slave1.type();
//		ComplexMat slave_tmp; slave_tmp.re = Mat::zeros(rows, cols, type); slave_tmp.im = Mat::zeros(rows, cols, type);
//#pragma omp parallel for schedule(guided)
//		for (int i = 0; i < rows; i++)
//		{
//			double x, y, iiii, jjjj; Mat tmp(1, 3, CV_64F); Mat result;
//			int mm0, nn0, mm1, nn1;
//			double offset_rows, offset_cols, upper, lower;
//			for (int j = 0; j < cols; j++)
//			{
//				jjjj = (double)j;
//				iiii = (double)i;
//				x = (jjjj - offset_x) / scale_x;
//				y = (iiii - offset_y) / scale_y;
//				tmp.at<double>(0, 0) = 1.0;
//				tmp.at<double>(0, 1) = x;
//				tmp.at<double>(0, 2) = y;
//				//tmp.at<double>(0, 3) = x * y;
//				//tmp.at<double>(0, 4) = x * x;
//				//tmp.at<double>(0, 5) = y * y;
//				result = tmp * coef_r;
//				offset_rows = result.at<double>(0, 0);
//				result = tmp * coef_c;
//				offset_cols = result.at<double>(0, 0);
//
//				iiii += offset_rows;
//				jjjj += offset_cols;
//
//				mm0 = (int)floor(iiii); nn0 = (int)floor(jjjj);
//				if (mm0 < 0 || nn0 < 0 || mm0 > rows_slave - 1 || nn0 > cols_slave - 1)
//				{
//					if (type == CV_64F)
//					{
//						slave_tmp.re.at<double>(i, j) = 0;
//						slave_tmp.im.at<double>(i, j) = 0;
//					}
//					else
//					{
//						slave_tmp.re.at<short>(i, j) = 0;
//						slave_tmp.im.at<short>(i, j) = 0;
//					}
//				}
//				else
//				{
//					mm1 = mm0 + 1; nn1 = nn0 + 1;
//					mm1 = mm1 >= rows_slave - 1 ? rows_slave - 1 : mm1;
//					nn1 = nn1 >= cols_slave - 1 ? cols_slave - 1 : nn1;
//					if (type == CV_16S)
//					{
//						//实部插值
//						upper = (double)slave1.re.at<short>(mm0, nn0) + double(slave1.re.at<short>(mm0, nn1) - slave1.re.at<short>(mm0, nn0)) * (jjjj - (double)nn0);
//						lower = (double)slave1.re.at<short>(mm1, nn0) + double(slave1.re.at<short>(mm1, nn1) - slave1.re.at<short>(mm1, nn0)) * (jjjj - (double)nn0);
//						slave_tmp.re.at<short>(i, j) = upper + double(lower - upper) * (iiii - (double)mm0);
//						//虚部插值
//						upper = (double)slave1.im.at<short>(mm0, nn0) + double(slave1.im.at<short>(mm0, nn1) - slave1.im.at<short>(mm0, nn0)) * (jjjj - (double)nn0);
//						lower = (double)slave1.im.at<short>(mm1, nn0) + double(slave1.im.at<short>(mm1, nn1) - slave1.im.at<short>(mm1, nn0)) * (jjjj - (double)nn0);
//						slave_tmp.im.at<short>(i, j) = upper + double(lower - upper) * (iiii - (double)mm0);
//					}
//					else
//					{
//						//实部插值
//						upper = slave1.re.at<double>(mm0, nn0) + (slave1.re.at<double>(mm0, nn1) - slave1.re.at<double>(mm0, nn0)) * (jjjj - (double)nn0);
//						lower = slave1.re.at<double>(mm1, nn0) + (slave1.re.at<double>(mm1, nn1) - slave1.re.at<double>(mm1, nn0)) * (jjjj - (double)nn0);
//						slave_tmp.re.at<double>(i, j) = upper + (lower - upper) * (iiii - (double)mm0);
//						//虚部插值
//						upper = slave1.im.at<double>(mm0, nn0) + (slave1.im.at<double>(mm0, nn1) - slave1.im.at<double>(mm0, nn0)) * (jjjj - (double)nn0);
//						lower = slave1.im.at<double>(mm1, nn0) + (slave1.im.at<double>(mm1, nn1) - slave1.im.at<double>(mm1, nn0)) * (jjjj - (double)nn0);
//						slave_tmp.im.at<double>(i, j) = upper + (lower - upper) * (iiii - (double)mm0);
//					}
//					
//				}
//
//			}
//		}
//
//		ret = conversion.write_slc_to_h5(SAR_images_out[ii].c_str(), slave_tmp);
//
//
//	}
//	return 0;
//}

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

int Utils::hist(Mat& input, double lowerbound, double upperbound, double interval, Mat& out_x, Mat& out_y)
{
	if (input.empty() ||
		input.channels() != 1 ||
		lowerbound >= upperbound ||
		interval >= (upperbound - lowerbound)
		)
	{
		fprintf(stderr, "hist(): input check failed!\n");
		return -1;
	}
	int num_bins = (upperbound - lowerbound) / interval + 1;
	double cmp = lowerbound + interval;
	out_x.create(1, num_bins, CV_64F);
	out_y.create(1, num_bins, CV_64F);
	for (int i = 0; i < num_bins; i++)
	{
		out_x.at<double>(i) = lowerbound + interval * (0.5 + i);
		out_y.at<double>(i) = 0.0;
	}
	Mat tmp;
	input.copyTo(tmp);
	if (tmp.rows != 1) tmp = tmp.reshape(0, 1);
	if (tmp.type() != CV_64F) tmp.convertTo(tmp, CV_64F);
	int num_element = tmp.cols; int count = 0; int bin_count = 0; int bin_num_count = 0;
	cv::sort(tmp, tmp, SORT_ASCENDING + SORT_EVERY_ROW);
	while (count < num_element)
	{
		if (tmp.at<double>(count) >= cmp)
		{
			out_y.at<double>(bin_count) = bin_num_count;
			bin_num_count = 0;
			bin_count++;
			cmp += interval;
		}
		bin_num_count++;
		count++;
	}
	out_y.at<double>(num_bins - 1) = out_y.at<double>(num_bins - 2);

	return 0;
}

int Utils::gaussian_curve_fit(Mat& input_x, Mat& input_y, double* mu, double* sigma_square, double* scale)
{
	if (input_x.empty() ||
		input_y.size() != input_x.size() ||
		input_x.type() != CV_64F ||
		input_y.type() != CV_64F ||
		!mu || !sigma_square || !scale
		)
	{
		fprintf(stderr, "gaussian_curve_fit(): input check failed!\n");
		return -1;
	}
	if (input_x.rows != 1)
	{
		input_x = input_x.reshape(0, 1);
		input_y = input_y.reshape(0, 1);
	}
	int cols = input_x.cols;
	Mat b(cols, 1, CV_64F), A(cols, 3, CV_64F); A = 1.0;
	for (int i = 0; i < cols; i++)
	{
		b.at<double>(i) = log(input_y.at<double>(i));
		A.at<double>(i, 1) = input_x.at<double>(i);
		A.at<double>(i, 2) = input_x.at<double>(i) * input_x.at<double>(i);
	}
	Mat A_t;
	cv::transpose(A, A_t);
	A = A_t * A;
	b = A_t * b;
	Mat x;
	if (!cv::solve(A, b, x, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "gaussian_curve_fit(): can't solve LS problem!\n");
		return -1;
	}
	*sigma_square = -1.0 / x.at<double>(2) / 2.0;
	*mu = *sigma_square * x.at<double>(1);
	*scale = exp(x.at<double>(0) + (*mu) * (*mu) / 2.0 / *sigma_square);
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
	int j_col = (int)cols / 2;
	row.create(rows, 1, CV_64F); col.create(rows, 1, CV_64F);
	for (int i = 0; i < rows; i++)
	{
		row.at<double>(i, 0) = i + offset_row;
		col.at<double>(i, 0) = j_col + offset_col;
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
		tmp.at<double>(0, 0) = lat.at<double>(i, 0);
		tmp.at<double>(0, 1) = lon.at<double>(i, 0);
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
		tmp.at<double>(0, 0) = lat.at<double>(i, 0);
		tmp.at<double>(0, 1) = lon.at<double>(i, 0);
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
		tmp.at<double>(0, 0) = lat.at<double>(i, 0);
		tmp.at<double>(0, 1) = lon.at<double>(i, 0);
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

int Utils::baseline_estimation(
	const Mat& stateVec1, 
	const Mat& stateVec2, 
	double lon_center,
	double lat_center, 
	int offset_row, 
	int offset_col, 
	int scene_height,
	int scene_width, 
	double time_interval,
	double time_interval2,
	double* B_effect,
	double* B_parallel
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.type() != CV_64F ||
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
	* 图像1成像点位置计算
	*/

	Mat sate1 = Mat::zeros(rows, 3, CV_64F);
	Mat sate2 = Mat::zeros(rows, 3, CV_64F);
	Mat satev1 = Mat::zeros(rows, 3, CV_64F);
	Mat satev2 = Mat::zeros(rows, 3, CV_64F);
	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat_center;
		tmp.at<double>(0, 1) = lon_center;
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
		tmp.at<double>(0, 0) = lat_center;
		tmp.at<double>(0, 1) = lon_center;
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
	Mat R, B, tmp, xyz, effect_dir; double r;
	tmp = Mat::zeros(1, 3, CV_64F);
	tmp.at<double>(0, 0) = lat_center;
	tmp.at<double>(0, 1) = lon_center;
	tmp.at<double>(0, 2) = 0;
	ell2xyz(tmp, xyz);

	R = xyz - sate1(Range(0, 1), Range(0, 3));
	r = sqrt(sum(R.mul(R))[0]);
	R = R / r;
	B = sate2(Range(0, 1), Range(0, 3)) - sate1(Range(0, 1), Range(0, 3));
	*B_parallel = sum(R.mul(B))[0];//平行基线

	tmp = satev1(Range(0, 1), Range(0, 3));
	cross(tmp, R, effect_dir);
	r = sqrt(sum(effect_dir.mul(effect_dir))[0]);
	effect_dir = effect_dir / r;
	r = sqrt(sum(xyz.mul(xyz))[0]);
	xyz = xyz / r;
	r = sum(xyz.mul(effect_dir))[0];
	effect_dir = r < 0.0 ? -effect_dir : effect_dir;
	*B_effect = sum(effect_dir.mul(B))[0];//平行基线


	return 0;
}

int Utils::baseline_estimation(
	const Mat& stateVec1,
	const Mat& stateVec2,
	double lon_center,
	double lat_center,
	double dem_center,
	int offset_row,
	int offset_col,
	int scene_height, 
	int scene_width, 
	double time_interval,
	double time_interval2,
	double* B_effect, 
	double* B_parallel
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.type() != CV_64F ||
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
	* 图像1成像点位置计算
	*/

	Mat sate1 = Mat::zeros(rows, 3, CV_64F);
	Mat sate2 = Mat::zeros(rows, 3, CV_64F);
	Mat satev1 = Mat::zeros(rows, 3, CV_64F);
	Mat satev2 = Mat::zeros(rows, 3, CV_64F);
	for (int i = 0; i < 1; i++)
	{
		Mat tmp(1, 3, CV_64F); Mat xyz;
		tmp.at<double>(0, 0) = lat_center;
		tmp.at<double>(0, 1) = lon_center;
		tmp.at<double>(0, 2) = dem_center;
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
		tmp.at<double>(0, 0) = lat_center;
		tmp.at<double>(0, 1) = lon_center;
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
	Mat R, B, tmp, xyz, effect_dir; double r;
	tmp = Mat::zeros(1, 3, CV_64F);
	tmp.at<double>(0, 0) = lat_center;
	tmp.at<double>(0, 1) = lon_center;
	tmp.at<double>(0, 2) = 0;
	ell2xyz(tmp, xyz);

	R = xyz - sate1(Range(0, 1), Range(0, 3));
	r = sqrt(sum(R.mul(R))[0]);
	R = R / r;
	B = sate2(Range(0, 1), Range(0, 3)) - sate1(Range(0, 1), Range(0, 3));
	*B_parallel = sum(R.mul(B))[0];//平行基线

	tmp = satev1(Range(0, 1), Range(0, 3));
	cross(tmp, R, effect_dir);
	r = sqrt(sum(effect_dir.mul(effect_dir))[0]);
	effect_dir = effect_dir / r;
	r = sqrt(sum(xyz.mul(xyz))[0]);
	xyz = xyz / r;
	r = sum(xyz.mul(effect_dir))[0];
	effect_dir = r < 0.0 ? -effect_dir : effect_dir;
	*B_effect = sum(effect_dir.mul(B))[0];//平行基线
	return 0;
}

int Utils::homogeneous_selection_and_phase_linking(
	vector<ComplexMat>& slc_stack,
	vector<ComplexMat>& slc_stack_filtered,
	int test_wndsize, int est_wndsize
)
{
	if (slc_stack.size() < 2 ||
		test_wndsize < 5 ||
		est_wndsize < 3 ||
		test_wndsize < est_wndsize ||
		test_wndsize % 2 != 1 ||
		est_wndsize % 2 != 1
		)
	{
		fprintf(stderr, "homogeneous_selection_and_phase_linking(): input check failed!\n");
		return -1;
	}
	if (slc_stack[0].type() != CV_32F)
	{
		fprintf(stderr, "homogeneous_selection_and_phase_linking(): input check failed!\n");
		return -1;
	}
	int n_images = slc_stack.size();
	//检查各图像数据尺寸是否一致
	int nr = slc_stack[0].GetRows();
	int nc = slc_stack[0].GetCols();
	for (int i = 0; i < n_images; i++)
	{
		if (slc_stack[i].GetCols() != nc || slc_stack[i].GetRows() != nr)
		{
			fprintf(stderr, "homogeneous_selection_and_phase_linking(): stack imagesize mismatch!\n");
			return -1;
		}
	}
	//将待处理数据拷贝至结果数据
	for (int i = 0; i < n_images; i++)
	{
		slc_stack_filtered.push_back(slc_stack[i]);
	}
	int radius_test = (test_wndsize - 1) / 2;
	int radius_estimation = (est_wndsize - 1) / 2;
	//扩展数据边缘
	 for (int i = 0; i < n_images; i++)
	{
		 cv::copyMakeBorder(slc_stack[i].re, slc_stack[i].re, radius_test + radius_estimation,
			 radius_test + radius_estimation, radius_test + radius_estimation, radius_test + radius_estimation,
			 BORDER_REFLECT);
		 cv::copyMakeBorder(slc_stack[i].im, slc_stack[i].im, radius_test + radius_estimation,
			 radius_test + radius_estimation, radius_test + radius_estimation, radius_test + radius_estimation,
			 BORDER_REFLECT);
	}
	

	double Look = 5.5;
	double thresh = 13.3;
	double rho = 0.8409;
	//预先计算第一块相关矩阵
	vector<ComplexMat> pre_covariance;
	ComplexMat temp((2 * radius_test + 1) * (nc + 2 * radius_test), 4);
	int x = n_images * (n_images - 1) / 2;
	for (int i = 0; i < x; i++)
	{
		pre_covariance.push_back(temp);
	}
#pragma omp parallel for schedule(guided)
	for (int mm = 0; mm < n_images; mm++)
	{
		for (int nn = mm + 1; nn < n_images; nn++)
		{
			int count = (n_images - 1 + n_images - mm) * mm / 2 + nn - mm - 1;
			ComplexMat Cov_t;
			for (int i = radius_estimation; i < 2 * radius_test + 1 + radius_estimation; i++)
			{
				for (int j = radius_estimation; j < nc + 2 * radius_test + radius_estimation; j++)
				{
					ComplexMat Cov(2, 2); ComplexMat vec(2, 1); ComplexMat vec_h, tmp_vec;
					for (int ii = i - radius_estimation; ii < i + 1 + radius_estimation; ii++)
					{
						for (int jj = j - radius_estimation; jj < j + 1 + radius_estimation; jj++)
						{
							vec.re.at<double>(0, 0) = slc_stack[mm].re.at<float>(ii, jj);
							vec.im.at<double>(0, 0) = slc_stack[mm].im.at<float>(ii, jj);
							vec.re.at<double>(1, 0) = slc_stack[nn].re.at<float>(ii, jj);
							vec.im.at<double>(1, 0) = slc_stack[nn].im.at<float>(ii, jj);
							vec_h = vec.transpose();
							vec.mul(vec_h, tmp_vec);
							Cov = Cov + tmp_vec;
						}
					}
					Cov.reshape(1, 4, Cov_t);
					int row = (i - radius_estimation) * (nc + 2 * radius_test) + j - radius_estimation;
					pre_covariance[count].SetValue(cv::Range(row, row + 1), cv::Range(0, 4), Cov_t);
				}
			}
		}
	}

	Mat coherence = Mat::zeros(nr, nc, CV_64F);
	Mat phase = Mat::zeros(nr, nc, CV_64F);
	Mat gamma = Mat::zeros(nr, nc, CV_64F);
	Mat homo_num = Mat::zeros(nr, nc, CV_32S);
	Mat homo_index = Mat::zeros(nr * nc, test_wndsize * test_wndsize, CV_16S);
	int total_count = 0;
	//预先对第一行元素的进行同质点选取和phase-linking
#pragma omp parallel for schedule(guided)
	for (int j = radius_test + radius_estimation; j < nc + radius_test + radius_estimation; j++)
	{
		ComplexMat Covariance, eigenvector;
		Mat Covariance_x1 = Mat::zeros(n_images, n_images, CV_64F);
		Mat Covariance_x2 = Mat::zeros(n_images, n_images, CV_64F);
		Mat eigenvalue;
		Covariance.re = Mat::eye(n_images, n_images, CV_64F);
		Covariance.im = Mat::eye(n_images, n_images, CV_64F);
		int count = 0;
		int i = radius_estimation + radius_test;
		for (int mm = 0; mm < n_images; mm++)
		{
			for (int nn = mm + 1; nn < n_images; nn++)
			{
				//选取中心点相关矩阵
				ComplexMat Cov_center, Cov_other, vec(2, 1), tmp;
				double a_re, a_im, b_re, b_im;
				int row = (i - radius_estimation) * (nc + 2 * radius_test) + j - radius_estimation;
				Cov_center = pre_covariance[count](cv::Range(row, row + 1), cv::Range(0, 4));
				Cov_center.reshape(2, 2, Cov_center);
				double det_center = Cov_center.determinant().real();
				a_re = slc_stack[mm].re.at<float>(i, j);
				a_im = slc_stack[mm].im.at<float>(i, j);
				b_re = slc_stack[nn].re.at<float>(i, j);
				b_im = slc_stack[nn].im.at<float>(i, j);
				Covariance.re.at<double>(mm, nn) += a_re * b_re + a_im * b_im;
				Covariance.im.at<double>(mm, nn) += a_im * b_re - a_re * b_im;
				Covariance_x1.at<double>(mm, nn) += a_re * a_re + a_im * a_im;
				Covariance_x2.at<double>(mm, nn) += b_re * b_re + b_im * b_im;
				Mat mask = Mat::zeros(test_wndsize, test_wndsize, CV_16S);
				mask.at<short>(radius_test, radius_test) = 1;
				//homo_num.at<int>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) += 1;
				//测试窗口内同质像元
				for (int ii = i - radius_test; ii < i + 1 + radius_test; ii++)
				{
					for (int jj = j - radius_test; jj < j + 1 + radius_test; jj++)
					{
						int row = (ii - radius_estimation) * (nc + 2 * radius_test) + jj - radius_estimation;
						Cov_other = pre_covariance[count](cv::Range(row, row + 1), cv::Range(0, 4));
						Cov_other.reshape(2, 2, Cov_other);
						double det_other = Cov_other.determinant().real();
						double det_sum = (Cov_other + Cov_center).determinant().real();
						double lnQ = Look * (2 * 2 * log(2) + log(det_center) + log(det_other) - 2.0 * log(det_sum));
						if (-2.0 * rho * lnQ <= thresh)
						{
							a_re = slc_stack[mm].re.at<float>(ii, jj);
							a_im = slc_stack[mm].im.at<float>(ii, jj);
							b_re = slc_stack[nn].re.at<float>(ii, jj);
							b_im = slc_stack[nn].im.at<float>(ii, jj);
							Covariance.re.at<double>(mm, nn) += a_re * b_re + a_im * b_im;
							Covariance.im.at<double>(mm, nn) += a_im * b_re - a_re * b_im;
							Covariance_x1.at<double>(mm, nn) += a_re * a_re + a_im * a_im;
							Covariance_x2.at<double>(mm, nn) += b_re * b_re + b_im * b_im;
							int xx = ii - i + radius_test;
							int yy = jj - j + radius_test;
							mask.at<short>(xx, yy) = 1;
							//homo_num.at<int>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) += 1;
						}
					}
				}
				int ix = i - radius_estimation - radius_test;
				int jx = j - radius_estimation - radius_test;
				cv::transpose(mask, mask);
				mask = mask.reshape(1, 1);
				mask.copyTo(homo_index(cv::Range(ix* nc + jx, ix* nc + jx + 1), cv::Range(0, test_wndsize* test_wndsize)));
				count++;
			}
		}
		for (int mm = 0; mm < n_images; mm++)
		{
			for (int nn = mm + 1; nn < n_images; nn++)
			{
				Covariance.re.at<double>(mm, nn) = Covariance.re.at<double>(mm, nn) / sqrt(Covariance_x1.at<double>(mm, nn) *
					Covariance_x2.at<double>(mm, nn));
				Covariance.im.at<double>(mm, nn) = Covariance.im.at<double>(mm, nn) / sqrt(Covariance_x1.at<double>(mm, nn) *
					Covariance_x2.at<double>(mm, nn));
				Covariance.re.at<double>(nn, mm) = Covariance.re.at<double>(mm, nn);
				Covariance.im.at<double>(nn, mm) = -Covariance.im.at<double>(mm, nn);
			}
		}
		//保存相位和相关系数（调试用）
		double re = Covariance.re.at<double>(0, 1);
		double im = Covariance.im.at<double>(0, 1);
		phase.at<double>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
			atan2(im, re);
		coherence.at<double>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
			sqrt(im * im + re * re);
		//phase-linking
		//int ret = HermitianEVD(Covariance, eigenvalue, eigenvector);
		//if (ret == 0)
		//{
		//	//计算goodness of fit
		//	/*for (int mm = 0; mm < n_images; mm++)
		//	{
		//		for (int nn = mm + 1; nn < n_images; nn++)
		//		{
		//			int rows = i - radius_estimation - radius_test;
		//			int cols = j - radius_estimation - radius_test;
		//			double re2, im2, re3, im3, phi_ik, phi_ik2;
		//			re = Covariance.re.at<double>(mm, nn);
		//			im = Covariance.im.at<double>(mm, nn);
		//			phi_ik = atan2(im, re);

		//			re = eigenvector.re.at<double>(mm, 0);;
		//			im = eigenvector.im.at<double>(mm, 0);
		//			re2 = eigenvector.re.at<double>(nn, 0);
		//			im2 = eigenvector.im.at<double>(nn, 0);

		//			re3 = re * re2 + im * im2;
		//			im3 = im * re2 - re * im2;
		//			phi_ik2 = atan2(im3, re3);
		//			gamma.at<double>(rows, cols) += cos(phi_ik - phi_ik2);
		//		}
		//	}*/
		//	for (int kk = 0; kk < n_images; kk++)
		//	{
		//		slc_stack_filtered[kk].re.at<float>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
		//			eigenvector.re.at<double>(kk, 0);
		//		slc_stack_filtered[kk].im.at<float>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
		//			eigenvector.im.at<double>(kk, 0);
		//	}
		//}
		
	}
	total_count++;
	printf("\r估计进度：%lf%%", double(total_count) / double(nr) * 100.0);
	fflush(stdout);
	//对剩余行进行循环处理


	for (int i = radius_estimation + radius_test + 1; i < nr + radius_estimation + radius_test; i++)
	{
		//更新相关矩阵块
#pragma omp parallel for schedule(guided)
		for (int mm = 0; mm < n_images; mm++)
		{
			for (int nn = mm + 1; nn < n_images; nn++)
			{
				int count = (n_images - 1 + n_images - mm) * mm / 2 + nn - mm - 1;
				//删除已经处理的第一行并将剩余部分上移一行
				pre_covariance[count].SetValue(cv::Range(0, (2 * radius_test)* (nc + 2 * radius_test)), cv::Range(0, 4),
					pre_covariance[count](cv::Range(nc + 2 * radius_test, (2 * radius_test + 1)* (nc + 2 * radius_test)),
						cv::Range(0, 4)));
				//重新计算新的一行并加入到相关矩阵块
				int i_last_row = i + radius_test;//新的一行行号
				for (int j = radius_estimation; j < nc + 2 * radius_test + radius_estimation; j++)
				{
					ComplexMat Cov(2, 2); ComplexMat Cov_t;
					for (int ii = i_last_row - radius_estimation; ii < i_last_row + 1 + radius_estimation; ii++)
					{
						for (int jj = j - radius_estimation; jj < j + 1 + radius_estimation; jj++)
						{
							ComplexMat vec(2, 1); ComplexMat vec_h, tmp_vec;
							vec.re.at<double>(0, 0) = slc_stack[mm].re.at<float>(ii, jj);
							vec.im.at<double>(0, 0) = slc_stack[mm].im.at<float>(ii, jj);
							vec.re.at<double>(1, 0) = slc_stack[nn].re.at<float>(ii, jj);
							vec.im.at<double>(1, 0) = slc_stack[nn].im.at<float>(ii, jj);
							vec_h = vec.transpose();
							vec.mul(vec_h, tmp_vec);
							Cov = Cov + tmp_vec;
						}
					}
					Cov.reshape(1, 4, Cov_t); 
					int row = (2 * radius_test) * (nc + 2 * radius_test) + j - radius_estimation;
					pre_covariance[count].SetValue(cv::Range(row, row + 1), cv::Range(0, 4), Cov_t);
				}
			}
		}


		//处理新的一行
#pragma omp parallel for schedule(guided)
		for (int j = radius_test + radius_estimation; j < nc + radius_test + radius_estimation; j++)
		{
			ComplexMat Covariance, eigenvector;
			Mat Covariance_x1 = Mat::zeros(n_images, n_images, CV_64F);
			Mat Covariance_x2 = Mat::zeros(n_images, n_images, CV_64F);
			Mat eigenvalue;
			Covariance.re = Mat::eye(n_images, n_images, CV_64F);
			Covariance.im = Mat::eye(n_images, n_images, CV_64F);
			int count = 0;
			//int i = radius_estimation + radius_test;
			for (int mm = 0; mm < n_images; mm++)
			{
				for (int nn = mm + 1; nn < n_images; nn++)
				{
					//选取中心点相关矩阵
					ComplexMat Cov_center, Cov_other, vec(2, 1), tmp;
					double a_re, a_im, b_re, b_im;
					/*注意此处总是处理中间行，因为相关矩阵块已经更新了*/
					int row = radius_test * (nc + 2 * radius_test) + j - radius_estimation;
					Cov_center = pre_covariance[count](cv::Range(row, row + 1), cv::Range(0, 4));
					Cov_center.reshape(2, 2, Cov_center);
					double det_center = Cov_center.determinant().real();
					a_re = slc_stack[mm].re.at<float>(i, j);
					a_im = slc_stack[mm].im.at<float>(i, j);
					b_re = slc_stack[nn].re.at<float>(i, j);
					b_im = slc_stack[nn].im.at<float>(i, j);
					Covariance.re.at<double>(mm, nn) += a_re * b_re + a_im * b_im;
					Covariance.im.at<double>(mm, nn) += a_im * b_re - a_re * b_im;
					Covariance_x1.at<double>(mm, nn) += a_re * a_re + a_im * a_im;
					Covariance_x2.at<double>(mm, nn) += b_re * b_re + b_im * b_im;
					Mat mask = Mat::zeros(test_wndsize, test_wndsize, CV_16S);
					mask.at<short>(radius_test, radius_test) = 1;
					//homo_num.at<int>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) += 1;
					//测试窗口内同质像元
					for (int ii = i - radius_test; ii < i + 1 + radius_test; ii++)
					{
						for (int jj = j - radius_test; jj < j + 1 + radius_test; jj++)
						{
							int row = (ii - i + radius_test) * (nc + 2 * radius_test) + jj - radius_estimation;
							Cov_other = pre_covariance[count](cv::Range(row, row + 1), cv::Range(0, 4));
							Cov_other.reshape(2, 2, Cov_other);
							double det_other = Cov_other.determinant().real();
							double det_sum = (Cov_other + Cov_center).determinant().real();
							double lnQ = Look * (2 * 2 * log(2) + log(det_center) + log(det_other) - 2.0 * log(det_sum));
							if (-2.0 * rho * lnQ <= thresh)
							{
								a_re = slc_stack[mm].re.at<float>(ii, jj);
								a_im = slc_stack[mm].im.at<float>(ii, jj);
								b_re = slc_stack[nn].re.at<float>(ii, jj);
								b_im = slc_stack[nn].im.at<float>(ii, jj);
								Covariance.re.at<double>(mm, nn) += a_re * b_re + a_im * b_im;
								Covariance.im.at<double>(mm, nn) += a_im * b_re - a_re * b_im;
								Covariance_x1.at<double>(mm, nn) += a_re * a_re + a_im * a_im;
								Covariance_x2.at<double>(mm, nn) += b_re * b_re + b_im * b_im;
								int xx = ii - i + radius_test;
								int yy = jj - j + radius_test;
								mask.at<short>(xx, yy) = 1;
								//homo_num.at<int>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) += 1;
							}
						}
					}
					int ix = i - radius_estimation - radius_test;
					int jx = j - radius_estimation - radius_test;
					cv::transpose(mask, mask);
					mask = mask.reshape(1, 1);
					mask.copyTo(homo_index(cv::Range(ix * nc + jx, ix * nc + jx + 1), cv::Range(0, test_wndsize * test_wndsize)));
					count++;
				}
			}
			for (int mm = 0; mm < n_images; mm++)
			{
				for (int nn = mm + 1; nn < n_images; nn++)
				{
					Covariance.re.at<double>(mm, nn) = Covariance.re.at<double>(mm, nn) / sqrt(Covariance_x1.at<double>(mm, nn) *
						Covariance_x2.at<double>(mm, nn));
					Covariance.im.at<double>(mm, nn) = Covariance.im.at<double>(mm, nn) / sqrt(Covariance_x1.at<double>(mm, nn) *
						Covariance_x2.at<double>(mm, nn));
					Covariance.re.at<double>(nn, mm) = Covariance.re.at<double>(mm, nn);
					Covariance.im.at<double>(nn, mm) = -Covariance.im.at<double>(mm, nn);
				}
			}

			//保存相位和相关系数（调试用）
			double re = Covariance.re.at<double>(0, 1);
			double im = Covariance.im.at<double>(0, 1);
			phase.at<double>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
				atan2(im, re);
			coherence.at<double>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
				sqrt(im * im + re * re);
			//phase-linking
			//int ret = HermitianEVD(Covariance, eigenvalue, eigenvector);
			//if (0)
			//{
			//	//计算goodness of fit
			//	/*for (int mm = 0; mm < n_images; mm++)
			//	{
			//		for (int nn = mm + 1; nn < n_images; nn++)
			//		{
			//			int rows = i - radius_estimation - radius_test;
			//			int cols = j - radius_estimation - radius_test;
			//			double re2, im2, re3, im3, phi_ik, phi_ik2;
			//			re = Covariance.re.at<double>(mm, nn);
			//			im = Covariance.im.at<double>(mm, nn);
			//			phi_ik = atan2(im, re);

			//			re = eigenvector.re.at<double>(mm, 0);;
			//			im = eigenvector.im.at<double>(mm, 0);
			//			re2 = eigenvector.re.at<double>(nn, 0);
			//			im2 = eigenvector.im.at<double>(nn, 0);

			//			re3 = re * re2 + im * im2;
			//			im3 = im * re2 - re * im2;
			//			phi_ik2 = atan2(im3, re3);
			//			gamma.at<double>(rows, cols) += cos(phi_ik - phi_ik2);
			//		}
			//	}*/
			//	for (int kk = 0; kk < n_images; kk++)
			//	{
			//		slc_stack_filtered[kk].re.at<float>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
			//			eigenvector.re.at<double>(kk, 0);
			//		slc_stack_filtered[kk].im.at<float>(i - radius_estimation - radius_test, j - radius_estimation - radius_test) =
			//			eigenvector.im.at<double>(kk, 0);
			//	}
			//}
		}
		total_count++;
		printf("\r估计进度：%lf%%", double(total_count) / double(nr) * 100.0);
		fflush(stdout);
	}
	//gamma = gamma / ((n_images - 1) * n_images / 2);
	//homo_index.convertTo(homo_index, CV_64F);
	//cvmat2bin("E:\\working_dir\\papers\\homogeneous_selection\\gamma.bin", gamma);
	FormatConversion conversion;
	conversion.creat_new_h5("E:\\working_dir\\papers\\homogeneous_selection\\Sentinel\\sceneA\\singlepol\\homo_index.h5");
	conversion.write_array_to_h5("E:\\working_dir\\papers\\homogeneous_selection\\Sentinel\\sceneA\\singlepol\\homo_index.h5",
		"homo_index", homo_index);
	cvmat2bin("E:\\working_dir\\papers\\homogeneous_selection\\Sentinel\\sceneA\\singlepol\\phase.bin", phase);
	cvmat2bin("E:\\working_dir\\papers\\homogeneous_selection\\Sentinel\\sceneA\\singlepol\\coherence.bin", coherence);
	return 0;
}

int Utils::SKP_decomposition(
	ComplexMat& inputMat,
	int nr1,
	int nc1,
	int nr2,
	int nc2,
	vector<ComplexMat>& outputMat1,
	vector<ComplexMat>& outputMat2
)
{
	assert(!inputMat.isempty() &&
		inputMat.type() == CV_64F &&
		(nr1 * nr2 == inputMat.GetRows()) &&
		(nc1 * nc2 == inputMat.GetCols()));
	
	//完成复矩阵reshape
	int nr_new = nr1 * nc1;
	int nc_new = nr2 * nc2;
	ComplexMat R(nr1 * nc1, nr2 * nc2);
	Mat tmp_re, tmp_im;
	int count = 0;
	for (int j = 0; j < nc1; j++)
	{
		for (int i = 0; i < nr1; i++)
		{
			inputMat.re(cv::Range(i * nr2, (i + 1) * nr2), cv::Range(j * nc2, (j + 1) * nc2)).copyTo(tmp_re);
			cv::transpose(tmp_re, tmp_re);
			tmp_re = tmp_re.reshape(1, 1);
			tmp_re.copyTo(R.re(cv::Range(count, count + 1), cv::Range(0, nr2 * nc2)));

			inputMat.im(cv::Range(i * nr2, (i + 1) * nr2), cv::Range(j * nc2, (j + 1) * nc2)).copyTo(tmp_im);
			cv::transpose(tmp_im, tmp_im);
			tmp_im = tmp_im.reshape(1, 1);
			tmp_im.copyTo(R.im(cv::Range(count, count + 1), cv::Range(0, nr2 * nc2)));

			count++;
		}
	}
	//cout << R.re << endl;
	//完成奇异值分解
	Eigen::MatrixXcd x(nr_new, nc_new);
	complex<double> d;
	for (int i = 0; i < nr_new; i++)
	{
		for (int j = 0; j < nc_new; j++)
		{
			d.real(R.re.at<double>(i, j));
			d.imag(R.im.at<double>(i, j));
			x(i, j) = d;
		}
	}
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
	svd.compute(x);
	Eigen::MatrixXcd U = svd.matrixU();
	Eigen::MatrixXcd V = svd.matrixV();
	Eigen::VectorXd singular_value = svd.singularValues();
	int M = singular_value.size();
	
	for (int i = 0; i < M; i++)
	{
		ComplexMat U_matrix(nr_new, 1), V_matrix(nc_new, 1);
		for (int j = 0; j < nr_new; j++)
		{
			U_matrix.re.at<double>(j, 0) = U(j, i).real();
			U_matrix.im.at<double>(j, 0) = U(j, i).imag();
		}
		for (int j = 0; j < nc_new; j++)
		{
			V_matrix.re.at<double>(j, 0) = V(j, i).real();
			V_matrix.im.at<double>(j, 0) = V(j, i).imag();
		}
		U_matrix = U_matrix * sqrt(singular_value(i, 0));
		V_matrix = V_matrix * sqrt(singular_value(i, 0));
		U_matrix.reshape(nr1, nc1, U_matrix);
		V_matrix.reshape(nr2, nc2, V_matrix);
		U_matrix = U_matrix.transpose(false);
		V_matrix = V_matrix.transpose();
		outputMat1.push_back(U_matrix);
		outputMat2.push_back(V_matrix);
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

int Utils::homogeneous_test(
	const vector<ComplexMat>& slc_series,
	int windsize_az,
	int windsize_rg,
	Mat& homo_num,
	Mat& homo_index,
	double alpha,
	const char* method
)
{
	if (slc_series.size() < 5 ||
		windsize_az < 3 ||
		windsize_az % 2 == 0 ||
		windsize_rg < 3 ||
		windsize_rg % 2 == 0
		)
	{
		fprintf(stderr, "homogeneous_test(): input check failed!\n");
		return -1;
	}
	int nr = slc_series[0].GetRows();
	int nc = slc_series[0].GetCols();
	int n_images = slc_series.size();
	int type = slc_series[0].type();
	for (int i = 0; i < n_images; i++)
	{
		int nr_tmp = slc_series[i].GetRows();
		int nc_tmp = slc_series[i].GetCols();
		int type_tmp = slc_series[i].type();
		if (type != type_tmp || nr_tmp != nr || nc_tmp != nc)
		{
			fprintf(stderr, "homogeneous_test(): images stack type or size mismatch!\n");
			return -1;
		}
	}
	homo_num.create(nr, nc, CV_32S); homo_num = 0;
	homo_index.create(nr * nc, windsize_az * windsize_rg, CV_8U); homo_index = 0;
	int radius_rg = (windsize_rg - 1) / 2;
	int radius_az = (windsize_az - 1) / 2;
	//统计同质检验
	if (type == CV_64F)
	{
		for (int m = 0; m < nr; m++)
		{
#pragma omp parallel for schedule(guided)
			for (int n = 0; n < nc; n++)
			{
				ComplexMat pix1(n_images, 1); ComplexMat pix2(n_images, 1);
				Mat pix1_amp, pix2_amp;
				Mat mask = Mat::zeros(windsize_az, windsize_rg, CV_32S);
				int ref_row, ref_col;
				ref_row = m, ref_col = n;
				mask.at<int>(radius_az, radius_rg) = 1;
				int count = 1;

				for (int i = 0; i < n_images; i++)
				{
					pix1.re.at<double>(i, 0) = slc_series[i].re.at<double>(ref_row, ref_col);
					pix1.im.at<double>(i, 0) = slc_series[i].im.at<double>(ref_row, ref_col);
				}
				pix1_amp = pix1.GetMod();
				for (int i = -radius_az; i <= radius_az; i++)
				{

					for (int j = -radius_rg; j <= radius_rg; j++)
					{
						if ((ref_row + i) < 0 || (ref_row + i) >= nr || (ref_col + j) < 0 || (ref_col + j) >= nc || (i == 0 && j == 0)) continue;
						for (int k = 0; k < n_images; k++)
						{
							pix2.re.at<double>(k, 0) = slc_series[k].re.at<double>(ref_row + i, ref_col + j);
							pix2.im.at<double>(k, 0) = slc_series[k].im.at<double>(ref_row + i, ref_col + j);
						}
						pix2_amp = pix2.GetMod();
						int b_homo = -1;
						homogeneous_test(pix1_amp, pix2_amp, &b_homo, alpha);
						if (b_homo == 0) { mask.at<int>(i + radius_az, j + radius_rg) = 1; count++; }
					}
				}
				homo_num.at<int>(m, n) = count;
				for (int i = 0; i < windsize_az; i++)
				{
					for (int j = 0; j < windsize_rg; j++)
					{
						homo_index.at<uchar>(m * nc + n, i * windsize_rg + j) = mask.at<int>(i, j);
					}
				}
			}
			printf("\r估计进度：%lf%%", double(m + 1) / double(nr) * 100.0);
			fflush(stdout);
		}
	}
	else if (type == CV_32F)
	{
		for (int m = 0; m < nr; m++)
		{
#pragma omp parallel for schedule(guided)
			for (int n = 0; n < nc; n++)
			{
				ComplexMat pix1(n_images, 1); ComplexMat pix2(n_images, 1);
				Mat pix1_amp, pix2_amp;
				Mat mask = Mat::zeros(windsize_az, windsize_rg, CV_32S);
				int ref_row, ref_col;
				ref_row = m, ref_col = n;
				mask.at<int>(radius_az, radius_rg) = 1;
				int count = 1;

				for (int i = 0; i < n_images; i++)
				{
					pix1.re.at<double>(i, 0) = slc_series[i].re.at<float>(ref_row, ref_col);
					pix1.im.at<double>(i, 0) = slc_series[i].im.at<float>(ref_row, ref_col);
				}
				pix1_amp = pix1.GetMod();
				for (int i = -radius_az; i <= radius_az; i++)
				{

					for (int j = -radius_rg; j <= radius_rg; j++)
					{
						if ((ref_row + i) < 0 || (ref_row + i) >= nr || (ref_col + j) < 0 || (ref_col + j) >= nc || (i == 0 && j == 0)) continue;
						for (int k = 0; k < n_images; k++)
						{
							pix2.re.at<double>(k, 0) = slc_series[k].re.at<float>(ref_row + i, ref_col + j);
							pix2.im.at<double>(k, 0) = slc_series[k].im.at<float>(ref_row + i, ref_col + j);
						}
						pix2_amp = pix2.GetMod();
						int b_homo = -1;
						homogeneous_test(pix1_amp, pix2_amp, &b_homo, alpha);
						if (b_homo == 0) { mask.at<int>(i + radius_az, j + radius_rg) = 1; count++; }
					}
				}
				homo_num.at<int>(m, n) = count;
				for (int i = 0; i < windsize_az; i++)
				{
					for (int j = 0; j < windsize_rg; j++)
					{
						homo_index.at<uchar>(m * nc + n, i * windsize_rg + j) = mask.at<int>(i, j);
					}
				}
			}
			printf("\r估计进度：%lf%%", double(m + 1) / double(nr) * 100.0);
			fflush(stdout);
		}
	}
	else if (type == CV_16S)
	{
		for (int m = 0; m < nr; m++)
		{
#pragma omp parallel for schedule(guided)
			for (int n = 0; n < nc; n++)
			{
				ComplexMat pix1(n_images, 1); ComplexMat pix2(n_images, 1);
				Mat pix1_amp, pix2_amp;
				Mat mask = Mat::zeros(windsize_az, windsize_rg, CV_32S);
				int ref_row, ref_col;
				ref_row = m, ref_col = n;
				mask.at<int>(radius_az, radius_rg) = 1;
				int count = 1;

				for (int i = 0; i < n_images; i++)
				{
					pix1.re.at<double>(i, 0) = slc_series[i].re.at<short>(ref_row, ref_col);
					pix1.im.at<double>(i, 0) = slc_series[i].im.at<short>(ref_row, ref_col);
				}
				pix1_amp = pix1.GetMod();
				for (int i = -radius_az; i <= radius_az; i++)
				{

					for (int j = -radius_rg; j <= radius_rg; j++)
					{
						if ((ref_row + i) < 0 || (ref_row + i) >= nr || (ref_col + j) < 0 || (ref_col + j) >= nc || (i == 0 && j == 0)) continue;
						for (int k = 0; k < n_images; k++)
						{
							pix2.re.at<double>(k, 0) = slc_series[k].re.at<short>(ref_row + i, ref_col + j);
							pix2.im.at<double>(k, 0) = slc_series[k].im.at<short>(ref_row + i, ref_col + j);
						}
						pix2_amp = pix2.GetMod();
						int b_homo = -1;
						homogeneous_test(pix1_amp, pix2_amp, &b_homo, alpha);
						if (b_homo == 0) { mask.at<int>(i + radius_az, j + radius_rg) = 1; count++; }
					}
				}
				homo_num.at<int>(m, n) = count;
				for (int i = 0; i < windsize_az; i++)
				{
					for (int j = 0; j < windsize_rg; j++)
					{
						homo_index.at<uchar>(m * nc + n, i * windsize_rg + j) = mask.at<int>(i, j);
					}
				}
			}
			printf("\r估计进度：%lf%%", double(m + 1) / double(nr) * 100.0);
			fflush(stdout);
		}
	}
	else
	{
		fprintf(stderr, "homogeneous_test(): data type not supported!\n");
		return -1;
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
	if (slc_series.size() < 3||
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
		//估计相关矩阵
		int count2 = rows * cols; int count;
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
					A.at<double>(k, count) = slc_series[k].re.at<double>(i + top, j + left);
					B.at<double>(k, count) = slc_series[k].im.at<double>(i + top, j + left);
					count++;
					s += A.at<double>(k, count - 1) * A.at<double>(k, count - 1)
						+ B.at<double>(k, count - 1) * B.at<double>(k, count - 1);
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
	return 0;
}

//int Utils::MB_phase_estimation(
//	vector<string> coregis_slc_files,
//	vector<string> phase_files, 
//	vector<string> coherence_files,
//	int master_indx, 
//	int blocksize_row, 
//	int blocksize_col, 
//	Mat& out_mask,
//	bool b_coh_est,
//	int homogeneous_test_wnd,
//	double thresh_c1_to_c2,
//	bool b_flat,
//	bool b_normalize
//)
//{
//	if (coregis_slc_files.size() < 2 ||
//		phase_files.size() != coregis_slc_files.size() ||
//		coherence_files.size() != phase_files.size() ||
//		master_indx < 1 ||
//		master_indx > phase_files.size() ||
//		blocksize_row < 100 ||
//		blocksize_col < 100 ||
//		thresh_c1_to_c2 < 0.0 ||
//		thresh_c1_to_c2 > 1.0 ||
//		homogeneous_test_wnd % 2 == 0
//		)
//	{
//		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
//		return -1;
//	}
//	int homotest_radius = (homogeneous_test_wnd - 1) / 2;
//	if (blocksize_row <= homotest_radius || blocksize_col <= homotest_radius)
//	{
//		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
//		return -1;
//	}
//	int nr, nc, ret, n_images; Mat tmp;
//	n_images = coregis_slc_files.size();
//	FormatConversion conversion; Deflat flat;
//	ret = conversion.read_array_from_h5(coregis_slc_files[0].c_str(), "azimuth_len", tmp);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	nr = tmp.at<int>(0, 0);
//	ret = conversion.read_array_from_h5(coregis_slc_files[0].c_str(), "range_len", tmp);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	nc = tmp.at<int>(0, 0);
//	Mat mask = Mat::zeros(nr, nc, CV_32S); mask.copyTo(out_mask); mask.release();
//	//检查输入SAR图像尺寸是否相同
//	for (int i = 1; i < n_images; i++)
//	{
//		ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "azimuth_len", tmp);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		if (tmp.at<int>(0, 0) != nr)
//		{
//			fprintf(stderr, "%s images size mismatch!\n", coregis_slc_files[i].c_str());
//			return -1;
//		}
//		ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "range_len", tmp);
//		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//		if (tmp.at<int>(0, 0) != nc)
//		{
//			fprintf(stderr, "%s images size mismatch!\n", coregis_slc_files[i].c_str());
//			return -1;
//		}
//	}
//
//	//去平地相位
//	Mat stateVec1, prf, lon_coef, lat_coef, 
//		carrier_frequency, stateVec2, prf2, phase,
//		phase_deflat, flat_phase_coef, azimuth_len, range_len;
//	phase = Mat::zeros(nr, nc, CV_64F);
//	azimuth_len = Mat::zeros(1, 1, CV_32S); range_len = Mat::zeros(1, 1, CV_32S);
//	int offset_row, offset_col;
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "offset_row", tmp);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	offset_row = tmp.at<int>(0, 0);
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "offset_col", tmp);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	offset_col = tmp.at<int>(0, 0);
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "state_vec", stateVec1);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "prf", prf);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "lon_coefficient", lon_coef);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "lat_coefficient", lat_coef);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	ret = conversion.read_array_from_h5(coregis_slc_files[master_indx - 1].c_str(), "carrier_frequency", carrier_frequency);
//	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//	for (int i = 0; i < n_images; i++)
//	{
//		//预先填充相关系数
//		if (b_coh_est)
//		{
//			ret = conversion.creat_new_h5(coherence_files[i].c_str());
//			if (return_check(ret, "creat_new_h5()", error_head)) return -1;
//			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "coherence", phase);
//			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//			azimuth_len.at<int>(0, 0) = phase.rows; range_len.at<int>(0, 0) = phase.cols;
//			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "azimuth_len", azimuth_len);
//			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//			ret = conversion.write_array_to_h5(coherence_files[i].c_str(), "range_len", range_len);
//			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		}
//		//预先填充干涉相位
//		ret = conversion.creat_new_h5(phase_files[i].c_str());
//		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
//		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "phase", phase);
//		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		azimuth_len.at<int>(0, 0) = phase.rows; range_len.at<int>(0, 0) = phase.cols;
//		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "azimuth_len", azimuth_len);
//		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		ret = conversion.write_array_to_h5(phase_files[i].c_str(), "range_len", range_len);
//		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		if (b_flat)
//		{
//			ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "state_vec", stateVec2);
//			if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//			ret = conversion.read_array_from_h5(coregis_slc_files[i].c_str(), "prf", prf2);
//			if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//			ret = flat.deflat(stateVec1, stateVec2, lon_coef, lat_coef, phase, offset_row, offset_col, 0, 1 / prf.at<double>(0, 0),
//				1 / prf2.at<double>(0, 0), 1, 3e8 / carrier_frequency.at<double>(0, 0), phase, flat_phase_coef);
//			ret = conversion.write_array_to_h5(phase_files[i].c_str(), "flat_phase_coefficient", flat_phase_coef);
//			if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
//		}
//		fprintf(stdout, "去平地进度：%d/%d\n", i, n_images - 1);
//	}
//
//	//分块读取、计算和储存
//
//	int left, right, top, bottom, block_num_row, block_num_col, left_pad, right_pad, top_pad, bottom_pad;
//	vector<ComplexMat> slc_series, slc_series_filter;
//	vector<Mat> coherence_series; coherence_series.resize(n_images);
//	ComplexMat slc, slc2, temp;
//	Mat flat_phase, ph, zeromat;
//	if (nr % blocksize_row == 0) block_num_row = nr / blocksize_row;
//	else block_num_row = int(floor((double)nr / (double)blocksize_row)) + 1;
//	if (nc % blocksize_col == 0) block_num_col = nc / blocksize_col;
//	else block_num_col = int(floor((double)nc / (double)blocksize_col)) + 1;
//	for (int i = 0; i < block_num_row; i++)
//	{
//		for (int j = 0; j < block_num_col; j++)
//		{
//			top = i * blocksize_row;
//			top_pad = top - homotest_radius; top_pad = top_pad < 0 ? 0 : top_pad;
//			bottom = top + blocksize_row; bottom = bottom > nr ? nr : bottom;
//			bottom_pad = bottom + homotest_radius; bottom_pad = bottom_pad > nr ? nr : bottom_pad;
//			left = j * blocksize_col;
//			left_pad = left - homotest_radius; left_pad = left_pad < 0 ? 0 : left_pad;
//			right = left + blocksize_col; right = right > nc ? nc : right;
//			right_pad = right + homotest_radius; right_pad = right_pad > nc ? nc : right_pad;
//
//			//读取数据
//			for (int k = 0; k < n_images; k++)
//			{
//				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_re",
//					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.re);
//				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_im",
//					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.im);
//				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
//				if (b_flat)
//				{
//					if (k != master_indx - 1)
//					{
//						flat_phase.create(bottom_pad - top_pad, right_pad - left_pad, CV_64F);
//						ret = conversion.read_array_from_h5(phase_files[k].c_str(), "flat_phase_coefficient", flat_phase_coef);
//						if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
//#pragma omp parallel for schedule(guided)
//						for (int ii = top_pad; ii < bottom_pad; ii++)
//						{
//							Mat tempp(1, 6, CV_64F);
//							for (int jj = left_pad; jj < right_pad; jj++)
//							{
//								tempp.at<double>(0, 0) = 1.0;
//								tempp.at<double>(0, 1) = ii;
//								tempp.at<double>(0, 2) = jj;
//								tempp.at<double>(0, 3) = ii * jj;
//								tempp.at<double>(0, 4) = ii * ii;
//								tempp.at<double>(0, 5) = jj * jj;
//								flat_phase.at<double>(ii - top_pad, jj - left_pad) = sum(tempp.mul(flat_phase_coef))[0];
//							}
//						}
//						ret = phase2cos(flat_phase, temp.re, temp.im);
//						if (return_check(ret, "phase2cos()", error_head)) return -1;
//						slc = slc * temp;
//					}
//					
//				}
//				slc_series.push_back(slc);
//				slc_series_filter.push_back(slc);
//			}
//
//			//填充相关系数
//			if (b_coh_est)
//			{
//				zeromat = Mat::zeros(flat_phase.rows, flat_phase.cols, CV_64F);
//				for (int mm = 0; mm < n_images; mm++)
//				{
//					zeromat.copyTo(coherence_series[mm]);
//				}
//			}
//			//计算
//#pragma omp parallel for schedule(guided)
//			for (int ii = (top - top_pad); ii < (bottom - top_pad); ii++)
//			{
//				ComplexMat coherence_matrix, eigenvector; Mat eigenvalue; int ret;
//				for (int jj = (left - left_pad); jj < (right - left_pad); jj++)
//				{
//					ret = coherence_matrix_estimation(slc_series, coherence_matrix, homogeneous_test_wnd, homogeneous_test_wnd, ii, jj);
//					if (ret == 0)
//					{
//						ret = HermitianEVD(coherence_matrix, eigenvalue, eigenvector);
//						if (!eigenvalue.empty() && ret == 0)
//						{
//							if (eigenvalue.at<double>(1, 0) / (eigenvalue.at<double>(0, 0) + 1e-10) < thresh_c1_to_c2)
//							{
//								out_mask.at<int>(ii + top_pad, jj + left_pad) = 1;
//								for (int kk = 0; kk < n_images; kk++)
//								{
//									slc_series_filter[kk].re.at<double>(ii, jj) = eigenvector.re.at<double>(kk, 0);
//									slc_series_filter[kk].im.at<double>(ii, jj) = eigenvector.im.at<double>(kk, 0);
//									if (b_coh_est)
//									{
//										coherence_series[kk].at<double>(ii, jj) = coherence_matrix(cv::Range(master_indx - 1, master_indx),
//											cv::Range(kk, kk + 1)).GetMod().at<double>(0, 0);
//									}
//									
//								}
//							}
//						}
//					}
//
//
//				}
//			}
//
//			//储存
//			slc = slc_series_filter[master_indx - 1];
//			for (int kk = 0; kk < n_images; kk++)
//			{
//				coherence_series[kk](cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
//				ret = conversion.write_subarray_to_h5(coherence_files[kk].c_str(), "coherence", ph, top, left, bottom - top, right - left);
//				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
//				//if (kk == master_indx - 1) continue;
//				ret = multilook(slc, slc_series_filter[kk], 1, 1, phase);
//				if (return_check(ret, "multilook()", error_head)) return -1;
//				phase(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
//				ret = conversion.write_subarray_to_h5(phase_files[kk].c_str(), "phase", ph, top, left, bottom - top, right - left);
//				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
//				
//			}
//			slc_series.clear();
//			slc_series_filter.clear();
//			
//			fprintf(stdout, "估计相位进度：%lf\n", double((i + 1) * block_num_col + j + 1) / double((block_num_col) * (block_num_row)));
//		}
//	}
//
//	return 0;
//}

int Utils::spatialTemporalBaselineEstimation(
	vector<string>& SLCH5Files,
	int reference,
	Mat& temporal,
	Mat& spatial
)
{
	if (SLCH5Files.size() < 2 ||
		reference < 1 || reference > SLCH5Files.size()
		)
	{
		fprintf(stderr, "spatialTemporalBaselineEstimation(): input check failed!\n");
		return -1;
	}
	Utils util; FormatConversion conversion;
	int ret, offset_row, offset_col, num_images, sceneHeight, sceneWidth;
	double prf1, prf2, acquisitionTime1, acquisitionTime2, B_temporal, B_spatial_para, B_spatial_effect;
	double topleft_lon, topright_lon, bottomleft_lon, bottomright_lon,
		topleft_lat, topright_lat, bottomleft_lat, bottomright_lat;
	Mat lon_coef, lat_coef, statevec1, statevec2;
	string start;
	num_images = SLCH5Files.size();
	temporal.create(1, num_images, CV_64F); spatial.create(1, num_images, CV_64F);
	temporal.at<double>(0, reference - 1) = 0; spatial.at<double>(0, reference - 1) = 0;
	ret = conversion.read_array_from_h5(SLCH5Files[reference - 1].c_str(), "lon_coefficient", lon_coef);
	ret = conversion.read_array_from_h5(SLCH5Files[reference - 1].c_str(), "lat_coefficient", lat_coef);
	ret = conversion.read_array_from_h5(SLCH5Files[reference - 1].c_str(), "state_vec", statevec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "prf", &prf1);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(SLCH5Files[reference - 1].c_str(), "acquisition_start_time", start);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start.c_str(), &acquisitionTime1);
	ret = conversion.read_int_from_h5(SLCH5Files[reference - 1].c_str(), "offset_row", &offset_row);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(SLCH5Files[reference - 1].c_str(), "offset_col", &offset_col);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(SLCH5Files[reference - 1].c_str(), "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(SLCH5Files[reference - 1].c_str(), "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	int ret2 = 0;
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "topLeftLon", &topleft_lon);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "topLeftLat", &topleft_lat);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "topRightLon", &topright_lon);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "topRightLat", &topright_lat);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "bottomLeftLon", &bottomleft_lon);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "bottomLeftLat", &bottomleft_lat);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "bottomRightLon", &bottomright_lon);
	ret2 += conversion.read_double_from_h5(SLCH5Files[reference - 1].c_str(), "bottomRightLat", &bottomright_lat);
	double lon_center, lat_center;
	if (ret2 == 0)
	{
		lon_center = (topleft_lon + topright_lon + bottomleft_lon + bottomright_lon) / 4.0;
		lat_center = (topleft_lat + topright_lat + bottomleft_lat + bottomright_lat) / 4.0;
	}

	for (int i = 0; i < num_images; i++)
	{
		if (i == reference - 1) continue;
		ret = conversion.read_array_from_h5(SLCH5Files[i].c_str(), "state_vec", statevec2);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		ret = conversion.read_str_from_h5(SLCH5Files[i].c_str(), "acquisition_start_time", start);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start.c_str(), &acquisitionTime2);
		ret = conversion.read_double_from_h5(SLCH5Files[i].c_str(), "prf", &prf2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		if (ret2 == 0)
		{
			util.baseline_estimation(statevec1, statevec2, lon_center, lat_center,
				offset_row, offset_col, sceneHeight, sceneWidth, 1.0 / prf1, 1.0 / prf2, &B_spatial_effect, &B_spatial_para);
			if (return_check(ret, "baseline_estimation()", error_head)) return -1;
		}
		else
		{
			ret = util.baseline_estimation(statevec1, statevec2, lon_coef, lat_coef, offset_row, offset_col,
				sceneHeight, sceneWidth, 1.0 / prf1, 1.0 / prf2, &B_spatial_effect, &B_spatial_para);
			if (return_check(ret, "baseline_estimation()", error_head)) return -1;
		}
		B_temporal = (acquisitionTime2 - acquisitionTime1) / 60 / 60 / 24;
		temporal.at<double>(0, i) = B_temporal;
		spatial.at<double>(0, i) = B_spatial_effect;
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

int Utils::computeImageGeoBoundry(
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
	if (lon_coefficient.rows != 1 ||
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

int Utils::computeImageGeoBoundry(
	double topleft_lon,
	double topleft_lat,
	double topright_lon, 
	double topright_lat,
	double bottomleft_lon, 
	double bottomleft_lat, 
	double bottomright_lon,
	double bottomright_lat, 
	double* lonMax,
	double* latMax,
	double* lonMin,
	double* latMin
)
{
	Mat lon, lat;
	lon.create(1, 4, CV_64F);
	lat.create(1, 4, CV_64F);

	lon.at<double>(0, 0) = topleft_lon;
	lon.at<double>(0, 1) = topright_lon;
	lon.at<double>(0, 2) = bottomleft_lon;
	lon.at<double>(0, 3) = bottomright_lon;

	lat.at<double>(0, 0) = topleft_lat;
	lat.at<double>(0, 1) = topright_lat;
	lat.at<double>(0, 2) = bottomleft_lat;
	lat.at<double>(0, 3) = bottomright_lat;

	cv::minMaxLoc(lon, lonMin, lonMax);
	cv::minMaxLoc(lat, latMin, latMax);
	double extra = 5.0 / 6000;
	*lonMin = *lonMin - extra * 30;
	*lonMax = *lonMax + extra * 30;
	*latMin = *latMin - extra * 30;
	*latMax = *latMax + extra * 30;
	return 0;
}

int Utils::getSRTMDEM(
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
	double latSpacing = 5.0 / 6000.0;
	double lonSpacing = 5.0 / 6000.0;
	if (!filepath || !lonUL || !latUL) return -1;
	//this->DEMPath = filepath;
	if (GetFileAttributesA(filepath) == -1)
	{
		if (_mkdir(filepath) != 0) return -1;
	}
	string DEMPath = filepath;
	vector<string> srtmFileName;
	vector<bool> bAlreadyExist;
	int ret = getSRTMFileName(lonMin, lonMax, latMin, latMax, srtmFileName);
	if (ret < 0)//不在SRTM数据范围内（-60°,60°）,则以0填充
	{
		int rows = (latMax - latMin) / latSpacing;
		int cols = (lonMax - lonMin) / lonSpacing;
		Mat temp = Mat::zeros(rows, cols, CV_16S);
		temp.copyTo(DEM_out);
		*lonUL = lonMin;
		*latUL = latMax;
		return 0;
	}
	//if (return_check(ret, "getSRTMFileName()", error_head)) return -1;
	//判断文件是否已经存在
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string tmp = DEMPath + "\\" + srtmFileName[i];
		std::replace(tmp.begin(), tmp.end(), '/', '\\');
		if (-1 != GetFileAttributesA(tmp.c_str()))bAlreadyExist.push_back(true);
		else bAlreadyExist.push_back(false);
	}
	//不存在则下载
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		if (!bAlreadyExist[i])
		{
			ret = downloadSRTM(srtmFileName[i].c_str(), DEMPath.c_str());
			if (ret < 0)//未下载到DEM数据,则以0填充
			{
				int rows = (latMax - latMin) / latSpacing;
				int cols = (lonMax - lonMin) / lonSpacing;
				Mat temp = Mat::zeros(rows, cols, CV_16S);
				temp.copyTo(DEM_out);
				*lonUL = lonMin;
				*latUL = latMax;
				return 0;
			}
		}
	}
	//解压文件
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string folderName = srtmFileName[i];
		folderName = folderName.substr(0, folderName.length() - 4);
		string path = DEMPath + string("\\") + folderName;
		std::replace(path.begin(), path.end(), '/', '\\');
		if (-1 != GetFileAttributesA(path.c_str())) continue;
		string srcFile = DEMPath + "\\" + srtmFileName[i];
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

		startRow = (latUpperLeft - latMax) / latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		startCol = (lonMin - lonUpperLeft) / lonSpacing;
		startCol = startCol < 1 ? 1 : startCol;
		startCol = startCol > total_cols ? total_cols : startCol;
		endCol = (lonMax - lonUpperLeft) / lonSpacing;
		endCol = endCol < 1 ? 1 : endCol;
		endCol = endCol > total_cols ? total_cols : endCol;

		string folderName = srtmFileName[0];
		folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
		string path = DEMPath + string("\\") + folderName;
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

			startRow = (latUpperLeft - latMax) / latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;


			Mat outDEM, outDEM2;

			if (yy < yy2)
			{
				string folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				string path = DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[1];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = DEMPath + string("\\") + folderName;
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
				string path = DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = DEMPath + string("\\") + folderName;
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
				startRow = (latUpperLeft - latMax) / latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMax - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMin - lonUpperLeft + 360.0) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;

				Mat outDEM, outDEM2;

				if (xx > xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = DEMPath + string("\\") + folderName;
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
					string path = DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = DEMPath + string("\\") + folderName;
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

				startRow = (latUpperLeft - latMax) / latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMin - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMax - lonUpperLeft) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;


				Mat outDEM, outDEM2;

				if (xx < xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = DEMPath + string("\\") + folderName;
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
					string path = DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = DEMPath + string("\\") + folderName;
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
			string path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);


			startRow = (latUpperLeft - latMax) / latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMax - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin - lonUpperLeft + 360.0) / lonSpacing;
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
			string path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
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
			path = DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = DigitalElevationModel::geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);

			startRow = (latUpperLeft - latMax) / latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		}

	}
	else return -1;
	return 0;
}

int Utils::getCopernicusDEM(
	const char* filepath,
	Mat& DEM_out,
	double* lonUL,
	double* latUL,
	double* lon_spacing,
	double* lat_spacing,
	double lonMin,
	double lonMax,
	double latMin,
	double latMax
)
{
	double latSpacing = 1.0 / 3600.0;
	double lonSpacing = 1.0 / 3600.0;
	if (!filepath || !lonUL || !latUL) return -1;
	if (GetFileAttributesA(filepath) == -1)
	{
		if (_mkdir(filepath) != 0) return -1;
	}
	string DEMPath = filepath;
	vector<string> CopernicusDEMFileName;
	vector<bool> bAlreadyExist;
	int ret = getCopernicusDEMFileName(lonMin, lonMax, latMin, latMax, CopernicusDEMFileName);
	if (ret < 0)//以0填充
	{
		int rows = (latMax - latMin) / latSpacing;
		int cols = (lonMax - lonMin) / lonSpacing;
		if (fabs(lonMax - lonMin) > 180.0)
		{
			cols = (- lonMax + lonMin + 360.0) / lonSpacing;
		}
		Mat temp = Mat::zeros(rows, cols, CV_32F);
		temp.copyTo(DEM_out);
		*lonUL = fabs(lonMax - lonMin) < 180.0 ? lonMin : lonMax;
		*latUL = latMax;
		*lon_spacing = lonSpacing;
		*lat_spacing = latSpacing;
		return 0;
	}
	//判断文件是否已经存在
	for (int i = 0; i < CopernicusDEMFileName.size(); i++)
	{
		string tmp = DEMPath + "\\" + CopernicusDEMFileName[i];
		std::replace(tmp.begin(), tmp.end(), '/', '\\');
		if (-1 != GetFileAttributesA(tmp.c_str()))bAlreadyExist.push_back(true);
		else bAlreadyExist.push_back(false);
	}
	//不存在则下载
	for (int i = 0; i < CopernicusDEMFileName.size(); i++)
	{
		if (!bAlreadyExist[i])
		{
			ret = downloadCopernicusDEM(CopernicusDEMFileName[i].c_str(), DEMPath.c_str());
			if (ret < 0)//未下载到DEM数据,则以0填充
			{
				int rows = (latMax - latMin) / latSpacing;
				int cols = (lonMax - lonMin) / lonSpacing;
				if (fabs(lonMax - lonMin) > 180.0)
				{
					cols = (-lonMax + lonMin + 360.0) / lonSpacing;
				}
				Mat temp = Mat::zeros(rows, cols, CV_32F);
				temp.copyTo(DEM_out);
				*lonUL = fabs(lonMax - lonMin) < 180.0 ? lonMin : lonMax;
				*latUL = latMax;
				*lon_spacing = lonSpacing;
				*lat_spacing = latSpacing;
				return 0;
			}
		}
	}



	int startRow, startCol, endRow, endCol;
	double lonUpperLeft, lonLowerRight, latUpperLeft, latLowerRight;
	int total_rows, total_cols;

	//DEM在一个SRTM方格内
	if (CopernicusDEMFileName.size() == 1)
	{
		int xx, yy;
		
		if (CopernicusDEMFileName[0].substr(22, 1) == string("N") && CopernicusDEMFileName[0].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx, &yy);
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("N") && CopernicusDEMFileName[0].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx, &yy);
			yy = -yy;
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("S") && CopernicusDEMFileName[0].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx, &yy);
			xx = -xx;
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("S") && CopernicusDEMFileName[0].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx, &yy);
			yy = -yy;
			xx = -xx;
		}
		else
		{
			return -1;
		}
		latUpperLeft = xx + 1.0;
		latLowerRight = xx;
		lonUpperLeft = yy;
		lonLowerRight = yy + 1.0;


		string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
		Mat outDEM;
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);
		latSpacing = 1.0 / double(outDEM.rows);
		lonSpacing = 1.0 / double(outDEM.cols);

		total_rows = outDEM.rows;
		total_cols = outDEM.cols;

		startRow = (latUpperLeft - latMax) / latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		if (fabs(lonMax - lonMin) > 180.0)
		{
			startCol = (lonMax - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin + 360.0 - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}
		else
		{
			startCol = (lonMin - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}
		

		outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
		*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
		*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		*lon_spacing = lonSpacing;
		*lat_spacing = latSpacing;
	}
	//DEM在2个方格内
	else if (CopernicusDEMFileName.size() == 2)
	{
		int xx, yy, xx2, yy2;

		if (CopernicusDEMFileName[0].substr(22, 1) == string("N") && CopernicusDEMFileName[0].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx, &yy);
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("N") && CopernicusDEMFileName[0].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx, &yy);
			yy = -yy;
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("S") && CopernicusDEMFileName[0].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx, &yy);
			xx = -xx;
		}
		else if (CopernicusDEMFileName[0].substr(22, 1) == string("S") && CopernicusDEMFileName[0].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[0].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx, &yy);
			yy = -yy;
			xx = -xx;
		}
		else
		{
			return -1;
		}

		if (CopernicusDEMFileName[1].substr(22, 1) == string("N") && CopernicusDEMFileName[1].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[1].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx2, &yy2);
		}
		else if (CopernicusDEMFileName[1].substr(22, 1) == string("N") && CopernicusDEMFileName[1].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[1].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx2, &yy2);
			yy2 = -yy2;
		}
		else if (CopernicusDEMFileName[1].substr(22, 1) == string("S") && CopernicusDEMFileName[1].substr(29, 1) == string("E"))
		{
			sscanf(CopernicusDEMFileName[1].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx2, &yy2);
			xx2 = -xx2;
		}
		else if (CopernicusDEMFileName[1].substr(22, 1) == string("S") && CopernicusDEMFileName[1].substr(29, 1) == string("W"))
		{
			sscanf(CopernicusDEMFileName[1].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx2, &yy2);
			yy2 = -yy2;
			xx2 = -xx2;
		}
		else
		{
			return -1;
		}

		//同一列
		if (yy == yy2)
		{
			Mat outDEM, outDEM2;

			string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);


			path = DEMPath + string("\\") + CopernicusDEMFileName[1];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);
			latSpacing = 1.0 / double(outDEM2.rows);
			lonSpacing = 1.0 / double(outDEM2.cols);

			if (outDEM.size() != outDEM2.size())
			{
				resize(outDEM, outDEM, outDEM2.size());
			}

			cv::vconcat(outDEM, outDEM2, outDEM);
			total_rows = outDEM.rows;
			total_cols = outDEM.cols;

			latUpperLeft = xx + 1;
			latLowerRight = xx - 1;
			lonUpperLeft = yy;
			lonLowerRight = yy + 1;

			startRow = (latUpperLeft - latMax) / latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			if (fabs(lonMax - lonMin) > 180.0)
			{
				startCol = (lonMax - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMin + 360.0 - lonUpperLeft) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;
			}
			else
			{
				startCol = (lonMin - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMax - lonUpperLeft) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;
			}

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
			*lon_spacing = lonSpacing;
			*lat_spacing = latSpacing;
		}
		//同一行
		else if (xx == xx2)
		{
			Mat outDEM, outDEM2;

			string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);
			

			path = DEMPath + string("\\") + CopernicusDEMFileName[1];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);
			latSpacing = 1.0 / double(outDEM2.rows);
			lonSpacing = 1.0 / double(outDEM2.cols);

			if (outDEM.size() != outDEM2.size())
			{
				resize(outDEM, outDEM, outDEM2.size());
			}

			cv::hconcat(outDEM, outDEM2, outDEM);
			total_rows = outDEM.rows;
			total_cols = outDEM.cols;

			lonUpperLeft = yy;
			lonLowerRight = yy + 2;
			latUpperLeft = xx + 1;
			latLowerRight = xx;

			startRow = (latUpperLeft - latMax) / latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			if (fabs(lonMax - lonMin) > 180.0)
			{
				startCol = (lonMax - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMin + 360.0 - lonUpperLeft) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;
			}
			else
			{
				startCol = (lonMin - lonUpperLeft) / lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMax - lonUpperLeft) / lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;
			}


			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
			*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
			*latUL = latUpperLeft - (startRow - 1) * latSpacing;
			*lon_spacing = lonSpacing;
			*lat_spacing = latSpacing;
		}
		else
		{
			return -1;
		}



	}
	//DEM在4个方格内
	else if (CopernicusDEMFileName.size() == 4)
	{
		int xx[4], yy[4], temp;

		for (int i = 0; i < 4; i++)
		{
			if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
				xx[i] = -xx[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
				xx[i] = -xx[i];
			}
			else
			{
				return -1;
			}
		}

		

		Mat outDEM, outDEM2, outDEM3;

		string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);


		path = DEMPath + string("\\") + CopernicusDEMFileName[1];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);


		cv::hconcat(outDEM, outDEM2, outDEM);

		path = DEMPath + string("\\") + CopernicusDEMFileName[2];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);
		latSpacing = 1.0 / double(outDEM2.rows);
		lonSpacing = 1.0 / double(outDEM2.cols);

		path = DEMPath + string("\\") + CopernicusDEMFileName[3];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);


		cv::hconcat(outDEM2, outDEM3, outDEM2);

		if (outDEM.size() != outDEM2.size())
		{
			resize(outDEM, outDEM, outDEM2.size());
		}

		cv::vconcat(outDEM, outDEM2, outDEM);

		total_rows = outDEM.rows;
		total_cols = outDEM.cols;



		latUpperLeft = xx[0] + 1;
		latLowerRight = xx[0] - 1;
		lonUpperLeft = yy[0];
		lonLowerRight = yy[0] + 2;

		startRow = (latUpperLeft - latMax) / latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		if (fabs(lonMax - lonMin) > 180.0)
		{
			startCol = (lonMax - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin + 360.0 - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}
		else
		{
			startCol = (lonMin - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}

		outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
		*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
		*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		*lon_spacing = lonSpacing;
		*lat_spacing = latSpacing;

	}
	//DEM在6个方格内	
	else if (CopernicusDEMFileName.size() == 6)
	{
		int xx[6], yy[6], temp;

		for (int i = 0; i < 6; i++)
		{
			if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
				xx[i] = -xx[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
				xx[i] = -xx[i];
			}
			else
			{
				return -1;
			}
		}
		Mat outDEM, outDEM2, outDEM3, outDEM4;
		//三行两列
		if (xx[0] != xx[2] && xx[0] == xx[1])
		{
			string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);


			path = DEMPath + string("\\") + CopernicusDEMFileName[1];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);

			cv::hconcat(outDEM, outDEM2, outDEM);

			path = DEMPath + string("\\") + CopernicusDEMFileName[2];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);

			
			path = DEMPath + string("\\") + CopernicusDEMFileName[3];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);
			
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			path = DEMPath + string("\\") + CopernicusDEMFileName[4];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);

			path = DEMPath + string("\\") + CopernicusDEMFileName[5];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM4);
			latSpacing = 1.0 / double(outDEM3.rows);
			lonSpacing = 1.0 / double(outDEM3.cols);

			cv::hconcat(outDEM3, outDEM4, outDEM3);

			

			if (outDEM.size() != outDEM3.size())
			{
				resize(outDEM, outDEM, outDEM3.size());
			}
			if (outDEM2.size() != outDEM3.size())
			{
				resize(outDEM2, outDEM2, outDEM3.size());
			}

			cv::vconcat(outDEM, outDEM2, outDEM);
			cv::vconcat(outDEM, outDEM3, outDEM);

			total_rows = outDEM.rows;
			total_cols = outDEM.cols;

			latUpperLeft = xx[0] + 1;
			latLowerRight = xx[0] - 2;
			lonUpperLeft = yy[0];
			lonLowerRight = yy[0] + 2;
		}
		//两行三列
		else if(xx[0] == xx[2] && xx[0] != xx[3])
		{
			string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);

			path = DEMPath + string("\\") + CopernicusDEMFileName[1];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);

			path = DEMPath + string("\\") + CopernicusDEMFileName[2];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);

			cv::hconcat(outDEM, outDEM2, outDEM);
			cv::hconcat(outDEM, outDEM3, outDEM);

			path = DEMPath + string("\\") + CopernicusDEMFileName[3];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);


			path = DEMPath + string("\\") + CopernicusDEMFileName[4];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);

			path = DEMPath + string("\\") + CopernicusDEMFileName[5];
			std::replace(path.begin(), path.end(), '/', '\\');
			ret = CopernicusDEM_geotiffread(path.c_str(), outDEM4);
			latSpacing = 1.0 / double(outDEM4.rows);
			lonSpacing = 1.0 / double(outDEM4.cols);

			cv::hconcat(outDEM2, outDEM3, outDEM2);
			cv::hconcat(outDEM2, outDEM4, outDEM2);




			if (outDEM.size() != outDEM2.size())
			{
				resize(outDEM, outDEM, outDEM2.size());
			}

			cv::vconcat(outDEM, outDEM2, outDEM);

			total_rows = outDEM.rows;
			total_cols = outDEM.cols;

			latUpperLeft = xx[0] + 1;
			latLowerRight = xx[0] - 1;
			lonUpperLeft = yy[0];
			lonLowerRight = yy[0] + 3;
		}




		

		startRow = (latUpperLeft - latMax) / latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		if (fabs(lonMax - lonMin) > 180.0)
		{
			startCol = (lonMax - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin + 360.0 - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}
		else
		{
			startCol = (lonMin - lonUpperLeft) / lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;
		}

		outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(DEM_out);
		*lonUL = lonUpperLeft + (startCol - 1) * lonSpacing;
		*latUL = latUpperLeft - (startRow - 1) * latSpacing;
		*lon_spacing = lonSpacing;
		*lat_spacing = latSpacing;
	}
	//DEM在9个方格内
	else if (CopernicusDEMFileName.size() == 9)
	{
		int xx[9], yy[9], temp;

		for (int i = 0; i < 9; i++)
		{
			if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("N") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_N%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("E"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_E%d_00_DEM.tif", &xx[i], &yy[i]);
				xx[i] = -xx[i];
			}
			else if (CopernicusDEMFileName[i].substr(22, 1) == string("S") && CopernicusDEMFileName[i].substr(29, 1) == string("W"))
			{
				sscanf(CopernicusDEMFileName[i].c_str(), "Copernicus_DSM_COG_10_S%d_00_W%d_00_DEM.tif", &xx[i], &yy[i]);
				yy[i] = -yy[i];
				xx[i] = -xx[i];
			}
			else
			{
				return -1;
			}
		}
		Mat outDEM, outDEM2, outDEM3, outDEM4, outDEM5;
		//三行三列
		string path = DEMPath + string("\\") + CopernicusDEMFileName[0];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM);


		path = DEMPath + string("\\") + CopernicusDEMFileName[1];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);

		

		path = DEMPath + string("\\") + CopernicusDEMFileName[2];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);

		cv::hconcat(outDEM, outDEM2, outDEM);
		cv::hconcat(outDEM, outDEM3, outDEM);


		path = DEMPath + string("\\") + CopernicusDEMFileName[3];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM2);


		path = DEMPath + string("\\") + CopernicusDEMFileName[4];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);



		path = DEMPath + string("\\") + CopernicusDEMFileName[5];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM4);

		cv::hconcat(outDEM2, outDEM3, outDEM2);
		cv::hconcat(outDEM2, outDEM4, outDEM2);

		path = DEMPath + string("\\") + CopernicusDEMFileName[6];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM3);


		path = DEMPath + string("\\") + CopernicusDEMFileName[7];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM4);



		path = DEMPath + string("\\") + CopernicusDEMFileName[8];
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = CopernicusDEM_geotiffread(path.c_str(), outDEM5);
		latSpacing = 1.0 / double(outDEM5.rows);
		lonSpacing = 1.0 / double(outDEM5.cols);

		cv::hconcat(outDEM3, outDEM4, outDEM3);
		cv::hconcat(outDEM3, outDEM5, outDEM3);



		if (outDEM.size() != outDEM3.size())
		{
			resize(outDEM, outDEM, outDEM3.size());
		}
		if (outDEM2.size() != outDEM3.size())
		{
			resize(outDEM2, outDEM2, outDEM3.size());
		}

		cv::vconcat(outDEM, outDEM2, outDEM);
		cv::vconcat(outDEM, outDEM3, outDEM);

		total_rows = outDEM.rows;
		total_cols = outDEM.cols;

		latUpperLeft = xx[0] + 1;
		latLowerRight = xx[0] - 2;
		lonUpperLeft = yy[0];
		lonLowerRight = yy[0] + 3;
	}
	else return -1;
	return 0;
}

int Utils::getSRTMFileName(double lonMin, double lonMax, double latMin, double latMax, vector<string>& name)
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

int Utils::getCopernicusDEMFileName(double lonMin, double lonMax, double latMin, double latMax, vector<string>& name)
{
	if (fabs(lonMin) > 180.0 ||
		fabs(lonMax) > 180.0 ||
		fabs(latMin) >= 90.0 ||
		fabs(latMax) >= 90.0
		)
	{
		fprintf(stderr, "getCopernicusDEMFileName(): input check failed!\n");
		return -1;
	}
	name.clear();
	char tmp[512];
	int startLat, endLat, startLon, startLon2, endLon, endLon2;
	startLat = ceil(latMax);
	endLat = floor(latMin);
	startLon = floor(lonMin);
	startLon2 = startLon;
	endLon = ceil(lonMax);
	endLon2 = endLon;
	if (fabs(lonMax - lonMin) > 180.0)//crossing the 180°/-180° longitude line
	{
		endLon2 = ceil(lonMax + fabs(lonMax - lonMin));
		startLon2 = floor(lonMax);
	}
	for (int i = int(startLat); i > int(endLat); i--)
	{
		string NS_sign = i - 1 >= 0 ? "N" : "S";
		memset(tmp, 0, 512);
		if (abs(i - 1) < 10)
		{
			sprintf(tmp, "0%d", abs(i - 1));
		}
		else
		{
			sprintf(tmp, "%d", abs(i - 1));
		}
		NS_sign = NS_sign + tmp;
		for (int j = int(startLon2); j < int(endLon2); j++)
		{	
			if (j >= 180) j = j - 360;
			string EW_sign = j >= 0 ? "E" : "W";
			memset(tmp, 0, 512);
			if (fabs(j) < 100)
			{
				sprintf(tmp, "0%d", abs(j));
			}
			else
			{
				sprintf(tmp, "%d", abs(j));
			}
			EW_sign = EW_sign + tmp;
			string filename = "Copernicus_DSM_COG_10_" + NS_sign + "_00_" + EW_sign + "_00_DEM.tif";
			name.push_back(filename);
		}
	}
	return 0;
}

int Utils::downloadSRTM(const char* name, const char* DEMpath)
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
	string url = string(SRTMURL) + name;
	string savefile = DEMpath + string("\\") + name;
	std::replace(savefile.begin(), savefile.end(), '/', '\\');
	HRESULT Result = URLDownloadToFileA(NULL, url.c_str(), savefile.c_str(), 0, NULL);
	if (Result != S_OK)
	{
		fprintf(stderr, "downloadSRTM(): download failded!\n");
		return -1;
	}
	return 0;
}

int Utils::downloadCopernicusDEM(const char* name, const char* DEMpath)
{
	bool isConnect;
	DWORD dw;
	isConnect = IsNetworkAlive(&dw);
	if (!isConnect)
	{
		fprintf(stderr, "downloadCopernicusDEM(): network is not connected!\n");
		return -1;
	}
	int ret;
	string folder = name;
	folder = folder.substr(0, folder.length() - 4);
	string url = string(CopernicusDEMURL) + folder + "/" + name;
	string savefile = DEMpath + string("/") + name;
	std::replace(savefile.begin(), savefile.end(), '/', '\\');
	HRESULT Result = URLDownloadToFileA(NULL, url.c_str(), savefile.c_str(), 0, NULL);
	if (Result != S_OK)
	{
		fprintf(stderr, "downloadCopernicusDEM(): download failded!\n");
		return -1;
	}
	return 0;
}

int Utils::CopernicusDEM_geotiffread(const char* filename, Mat& outDEM)
{
	if (!filename)
		return -1;

	GDALAllRegister();   // 注册驱动

	GDALDatasetH hDataset = GDALOpen(filename, GA_ReadOnly);
	if (hDataset == NULL)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): failed to open %s!\n",
			filename);
		return -1;
	}

	int nBand = GDALGetRasterCount(hDataset);
	if (nBand != 1)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): number of Bands != 1\n");
		GDALClose(hDataset);
		return -1;
	}

	GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);
	if (hBand == NULL)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): failed to get band!\n");
		GDALClose(hDataset);
		return -1;
	}

	int xsize = GDALGetRasterBandXSize(hBand);   // cols
	int ysize = GDALGetRasterBandYSize(hBand);   // rows

	if (xsize <= 0 || ysize <= 0)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): band rows and cols error!\n");
		GDALClose(hDataset);
		return -1;
	}

	GDALDataType dataType = GDALGetRasterDataType(hBand);
	/* Copernicus DEM 通常是 GDT_Int16 或 GDT_Float32 */

	float* pbuf = (float*)malloc(sizeof(float) * xsize * ysize);
	if (!pbuf)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): out of memory!\n");
		GDALClose(hDataset);
		return -1;
	}

	if (GDALRasterIO(
		hBand,
		GF_Read,
		0, 0,
		xsize, ysize,
		pbuf,
		xsize, ysize,
		GDT_Float32,   /* 直接转成 float */
		0, 0) != CE_None)
	{
		fprintf(stderr,
			"CopernicusDEM_geotiffread(): RasterIO failed!\n");
		free(pbuf);
		GDALClose(hDataset);
		return -1;
	}

	outDEM.create(ysize, xsize, CV_32F);
	memcpy(outDEM.data, pbuf, sizeof(float) * xsize * ysize);

	free(pbuf);
	GDALClose(hDataset);

	return 0;
}

int Utils::writeOverlayKML(
	double BottomLeft_lon,
	double BottomLeft_lat,
	double BottomRight_lon,
	double BottomRight_lat,
	double TopRight_lon,
	double TopRight_lat, 
	double TopLeft_lon,
	double TopLeft_lat, 
	double Reference_lon,
	double Reference_lat,
	const char* image_file,
	const char* KML_file, 
	const char* Legend_file
)
{
	if (fabs(BottomLeft_lat) > 90.0 ||
		fabs(BottomRight_lat) > 90.0 ||
		fabs(TopLeft_lat) > 90.0 ||
		fabs(TopRight_lat) > 90.0 ||
		fabs(Reference_lat) > 90.0 ||
		fabs(BottomLeft_lon) > 180.0 ||
		fabs(BottomRight_lon) > 180.0 ||
		fabs(TopLeft_lon) > 180.0 ||
		fabs(TopRight_lon) > 180.0 ||
		fabs(Reference_lon) > 180.0 ||
		!image_file ||
		!KML_file
		)
	{
		fprintf(stderr, "writeOverlayKML(): input check failed!\n");
		return -1;
	}
	TiXmlDocument doc;
	TiXmlDeclaration* declaration = new TiXmlDeclaration("1.0", "UTF-8", "yes");
	doc.LinkEndChild(declaration);
	TiXmlElement* Root = new TiXmlElement("kml");
	Root->SetAttribute("xmlns:gx", "http://www.google.com/kml/ext/2.2");
	doc.LinkEndChild(Root);
	TiXmlElement* Document = new TiXmlElement("Document");
	Root->LinkEndChild(Document);

	TiXmlElement* name = new TiXmlElement("name");
	TiXmlText* content = new TiXmlText("SatExplorer Product Map Overlay");
	name->LinkEndChild(content);
	Document->LinkEndChild(name);

	TiXmlElement* Folder = new TiXmlElement("Folder");
	Document->LinkEndChild(Folder);

	TiXmlElement* name2 = new TiXmlElement("name");
	TiXmlText* content2 = new TiXmlText("SatExplorer Product Scene Overlay");
	name2->LinkEndChild(content2);
	Folder->LinkEndChild(name2);

	TiXmlElement* GroundOverlay = new TiXmlElement("GroundOverlay");
	Folder->LinkEndChild(GroundOverlay);

	TiXmlElement* name3 = new TiXmlElement("name");
	TiXmlText* content3 = new TiXmlText("SatExplorer Product Image Overlay");
	name3->LinkEndChild(content3);
	GroundOverlay->LinkEndChild(name3);

	TiXmlElement* Icon = new TiXmlElement("Icon");
	GroundOverlay->LinkEndChild(Icon);

	TiXmlElement* href1 = new TiXmlElement("href");
	TiXmlText* content4 = new TiXmlText(image_file);
	href1->LinkEndChild(content4);
	Icon->LinkEndChild(href1);

	TiXmlElement* LatLonQuad = new TiXmlElement("gx:LatLonQuad");
	GroundOverlay->LinkEndChild(LatLonQuad);

	TiXmlElement* coordinates = new TiXmlElement("coordinates");

	char coordinates_content[1024];
	sprintf(coordinates_content, "%lf,%lf %lf,%lf %lf,%lf %lf,%lf", BottomLeft_lon, BottomLeft_lat,
		BottomRight_lon, BottomRight_lat, TopRight_lon, TopRight_lat, TopLeft_lon, TopLeft_lat);

	TiXmlText* content5 = new TiXmlText(coordinates_content);
	coordinates->LinkEndChild(content5);
	LatLonQuad->LinkEndChild(coordinates);
	
	//如果有图例文件
	if (Legend_file)
	{
		TiXmlElement* Folder2 = new TiXmlElement("Folder");
		Document->LinkEndChild(Folder2);

		TiXmlElement* name22 = new TiXmlElement("name");
		TiXmlText* content22 = new TiXmlText("SatExplorer Product Screen Overlays");
		name22->LinkEndChild(content22);
		Folder2->LinkEndChild(name22);

		TiXmlElement* ScreenOverlay = new TiXmlElement("ScreenOverlay");
		Folder2->LinkEndChild(ScreenOverlay);

		TiXmlElement* name32 = new TiXmlElement("name");
		TiXmlText* content32 = new TiXmlText("SatExplorer Product Image Legend");
		name32->LinkEndChild(content32);
		ScreenOverlay->LinkEndChild(name32);

		TiXmlElement* Icon2 = new TiXmlElement("Icon");
		ScreenOverlay->LinkEndChild(Icon2);

		TiXmlElement* href2 = new TiXmlElement("href");
		TiXmlText* content42 = new TiXmlText(Legend_file);
		href2->LinkEndChild(content42);
		Icon2->LinkEndChild(href2);

		TiXmlElement* overlayXY = new TiXmlElement("overlayXY");
		overlayXY->SetAttribute("x", "0");
		overlayXY->SetAttribute("y", "1");
		overlayXY->SetAttribute("xunits", "fraction");
		overlayXY->SetAttribute("yunits", "fraction");
		ScreenOverlay->LinkEndChild(overlayXY);

		TiXmlElement* screenXY = new TiXmlElement("screenXY");
		screenXY->SetAttribute("x", "0");
		screenXY->SetAttribute("y", "1");
		screenXY->SetAttribute("xunits", "fraction");
		screenXY->SetAttribute("yunits", "fraction");
		ScreenOverlay->LinkEndChild(screenXY);

		TiXmlElement* rotationXY = new TiXmlElement("rotationXY");
		rotationXY->SetAttribute("x", "0");
		rotationXY->SetAttribute("y", "0");
		rotationXY->SetAttribute("xunits", "fraction");
		rotationXY->SetAttribute("yunits", "fraction");
		ScreenOverlay->LinkEndChild(rotationXY);

		TiXmlElement* size = new TiXmlElement("size");
		size->SetAttribute("x", "0");
		size->SetAttribute("y", "0");
		size->SetAttribute("xunits", "fraction");
		size->SetAttribute("yunits", "fraction");
		ScreenOverlay->LinkEndChild(size);
	}

	//参考点标记
	TiXmlElement* Placemark = new TiXmlElement("Placemark");
	Document->LinkEndChild(Placemark);

	TiXmlElement* name23 = new TiXmlElement("name");
	TiXmlText* content23 = new TiXmlText("Reference Point");
	name23->LinkEndChild(content23);
	Placemark->LinkEndChild(name23);

	TiXmlElement* Point = new TiXmlElement("Point");
	Placemark->LinkEndChild(Point);

	TiXmlElement* coordinates_ref = new TiXmlElement("coordinates");
	sprintf(coordinates_content, "%lf,%lf", Reference_lon, Reference_lat);
	TiXmlText* content55 = new TiXmlText(coordinates_content);
	coordinates_ref->LinkEndChild(content55);
	Point->LinkEndChild(coordinates_ref);


	if (!doc.SaveFile(KML_file))
	{
		fprintf(stderr, "writeOverlayKML(): failed to save %s! \n", KML_file);
		return -1;
	}
	return 0;
}

int Utils::S1_subswath_merge(
	const char* IW1_h5file,
	const char* IW2_h5file,
	const char* IW3_h5file, 
	const char* merged_phase_h5file
)
{
	if (!IW1_h5file || !IW2_h5file || !IW3_h5file || !merged_phase_h5file)
	{
		fprintf(stderr, "S1_subswath_merge(): input check failed\n");
		return -1;
	}
	double start1, start2, start3, end1, end2, end3, first_pixel1, first_pixel2, first_pixel3,
		range_spacing, prf;
	int mul_az, mul_rg, mul_az1, mul_rg1, rows1, rows2, rows3, cols1, cols2, cols3;
	FormatConversion conversion;
	string start_time, end_time;
	int ret;
	ret = conversion.read_double_from_h5(IW1_h5file, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW1_h5file, "slant_range_first_pixel", &first_pixel1);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW2_h5file, "slant_range_first_pixel", &first_pixel2);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW3_h5file, "slant_range_first_pixel", &first_pixel3);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	if (first_pixel1 >= first_pixel2 || first_pixel2 >= first_pixel3)
	{
		fprintf(stderr, "S1_subswath_merge(): please rearrange input swath order!\n");
		return -1;
	}
	ret = conversion.read_double_from_h5(IW3_h5file, "range_spacing", &range_spacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW1_h5file, "multilook_az", &mul_az);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW1_h5file, "multilook_rg", &mul_rg);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW2_h5file, "multilook_az", &mul_az1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW2_h5file, "multilook_rg", &mul_rg1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	if (mul_az != mul_az1 || mul_rg != mul_rg1)
	{
		fprintf(stderr, "S1_subswath_merge(): multilook times disagree!\n");
		return -1;
	}
	ret = conversion.read_int_from_h5(IW3_h5file, "multilook_az", &mul_az1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW3_h5file, "multilook_rg", &mul_rg1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	if (mul_az != mul_az1 || mul_rg != mul_rg1)
	{
		fprintf(stderr, "S1_subswath_merge(): multilook times disagree!\n");
		return -1;
	}
	
	ret = conversion.read_int_from_h5(IW1_h5file, "range_len", &cols1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW1_h5file, "azimuth_len", &rows1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW2_h5file, "range_len", &cols2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW2_h5file, "azimuth_len", &rows2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW3_h5file, "range_len", &cols3);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW3_h5file, "azimuth_len", &rows3);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;



	ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start1);
	conversion.utc2gps(end_time.c_str(), &end1);

	ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start2);
	conversion.utc2gps(end_time.c_str(), &end2);

	ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start3);
	conversion.utc2gps(end_time.c_str(), &end3);

	//IW1和IW2拼接

	int x = (double(cols1 * mul_rg) - round((first_pixel2 - first_pixel1) / range_spacing)) / 2.0;
	int tmp = x + (int)round((first_pixel2 - first_pixel1) / range_spacing);
	int col_end_last = tmp / mul_rg;
	int col_start_next = (col_end_last * mul_rg - round((first_pixel2 - first_pixel1) / range_spacing)) / mul_rg;
	int total_cols = col_end_last + (cols2 - col_start_next);
	int row_offset = (start1 - start2) * prf / (double)mul_az;
	int total_rows = (((end1 > end2 ? end1 : end2) - (start1 < start2 ? start1 : start2))*prf + 1) / (double)mul_az + 1;
	Mat phase_tmp(total_rows, total_cols, CV_64F); phase_tmp = 0.0;
	Mat lon_tmp(total_rows, total_cols, CV_32F); lon_tmp = 360.0;
	Mat lat_tmp(total_rows, total_cols, CV_32F); lat_tmp = 360.0;
	Mat phase1, phase2, mapped_lon1, mapped_lon2, mapped_lat1, mapped_lat2;
	ret = conversion.read_array_from_h5(IW1_h5file, "phase", phase1);
	ret = conversion.read_array_from_h5(IW2_h5file, "phase", phase2);
	ret = conversion.read_array_from_h5(IW1_h5file, "mapped_lon", mapped_lon1);
	ret += conversion.read_array_from_h5(IW1_h5file, "mapped_lat", mapped_lat1);
	ret += conversion.read_array_from_h5(IW2_h5file, "mapped_lon", mapped_lon2);
	ret += conversion.read_array_from_h5(IW2_h5file, "mapped_lat", mapped_lat2);
	int temp_ret = ret;
	if (row_offset < 0)
	{
		if (temp_ret == 0)
		{
			mapped_lon1(cv::Range(0, mapped_lon1.rows), cv::Range(0, col_end_last)).copyTo
			(lon_tmp(cv::Range(0, mapped_lon1.rows), cv::Range(0, col_end_last)));
			mapped_lon2(cv::Range(0, mapped_lon2.rows), cv::Range(col_start_next, cols2)).copyTo
			(lon_tmp(cv::Range(-row_offset, mapped_lon2.rows - row_offset), cv::Range(col_end_last, total_cols)));

			mapped_lat1(cv::Range(0, mapped_lat1.rows), cv::Range(0, col_end_last)).copyTo
			(lat_tmp(cv::Range(0, mapped_lat1.rows), cv::Range(0, col_end_last)));
			mapped_lat2(cv::Range(0, mapped_lat2.rows), cv::Range(col_start_next, cols2)).copyTo
			(lat_tmp(cv::Range(-row_offset, mapped_lat2.rows - row_offset), cv::Range(col_end_last, total_cols)));
		}
		phase1(cv::Range(0, phase1.rows), cv::Range(0, col_end_last)).copyTo
		(phase_tmp(cv::Range(0, phase1.rows), cv::Range(0, col_end_last)));
		phase2(cv::Range(0, phase2.rows), cv::Range(col_start_next, cols2)).copyTo
		(phase_tmp(cv::Range(-row_offset, phase2.rows - row_offset), cv::Range(col_end_last, total_cols)));
	}
	else
	{
		if (temp_ret == 0)
		{
			mapped_lon1(cv::Range(0, mapped_lon1.rows), cv::Range(0, col_end_last)).copyTo
			(lon_tmp(cv::Range(row_offset, mapped_lon1.rows + row_offset), cv::Range(0, col_end_last)));
			mapped_lon2(cv::Range(0, mapped_lon2.rows), cv::Range(col_start_next, cols2)).copyTo
			(lon_tmp(cv::Range(0, mapped_lon2.rows), cv::Range(col_end_last, total_cols)));

			mapped_lat1(cv::Range(0, mapped_lat1.rows), cv::Range(0, col_end_last)).copyTo
			(lat_tmp(cv::Range(row_offset, mapped_lat1.rows + row_offset), cv::Range(0, col_end_last)));
			mapped_lat2(cv::Range(0, mapped_lat2.rows), cv::Range(col_start_next, cols2)).copyTo
			(lat_tmp(cv::Range(0, mapped_lat2.rows), cv::Range(col_end_last, total_cols)));
		}
		phase1(cv::Range(0, phase1.rows), cv::Range(0, col_end_last)).copyTo
		(phase_tmp(cv::Range(row_offset, phase1.rows + row_offset), cv::Range(0, col_end_last)));
		phase2(cv::Range(0, phase2.rows), cv::Range(col_start_next, cols2)).copyTo
		(phase_tmp(cv::Range(0, phase2.rows), cv::Range(col_end_last, total_cols)));
	}

	//拼接IW3

	x = (double(total_cols * mul_rg) - round((first_pixel3 - first_pixel1) / range_spacing)) / 2.0;
	tmp = x + (int)round((first_pixel3 - first_pixel1) / range_spacing);
	col_end_last = tmp / mul_rg;
	col_start_next = (col_end_last * mul_rg - round((first_pixel3 - first_pixel1) / range_spacing)) / mul_rg;
	total_cols = col_end_last + (cols3 - col_start_next);
	start1 = start1 < start2 ? start1 : start2;
	end1 = end1 > end2 ? end1 : end2;
	row_offset = (start1 - start3) * prf / (double)mul_az;
	total_rows = (((end1 > end3 ? end1 : end3) - (start1 < start3 ? start1 : start3)) * prf + 1) / (double)mul_az + 1;

	phase1.create(total_rows, total_cols, CV_64F); phase1 = 0.0;
	mapped_lat2.create(total_rows, total_cols, CV_32F); mapped_lat2 = 360.0;
	mapped_lon2.create(total_rows, total_cols, CV_32F); mapped_lon2 = 360.0;
	ret = conversion.read_array_from_h5(IW3_h5file, "phase", phase2);
	ret = conversion.read_array_from_h5(IW3_h5file, "mapped_lon", mapped_lon1);
	ret += conversion.read_array_from_h5(IW3_h5file, "mapped_lat", mapped_lat1);
	temp_ret += ret;
	if (row_offset < 0)
	{
		if (temp_ret == 0)
		{
			lon_tmp(cv::Range(0, lon_tmp.rows), cv::Range(0, col_end_last)).copyTo
			(mapped_lon2(cv::Range(0, lon_tmp.rows), cv::Range(0, col_end_last)));
			mapped_lon1(cv::Range(0, mapped_lon1.rows), cv::Range(col_start_next, cols3)).copyTo
			(mapped_lon2(cv::Range(-row_offset, mapped_lon1.rows - row_offset), cv::Range(col_end_last, total_cols)));

			lat_tmp(cv::Range(0, lat_tmp.rows), cv::Range(0, col_end_last)).copyTo
			(mapped_lat2(cv::Range(0, lat_tmp.rows), cv::Range(0, col_end_last)));
			mapped_lat1(cv::Range(0, mapped_lat1.rows), cv::Range(col_start_next, cols3)).copyTo
			(mapped_lat2(cv::Range(-row_offset, mapped_lat1.rows - row_offset), cv::Range(col_end_last, total_cols)));
		}
		phase_tmp(cv::Range(0, phase_tmp.rows), cv::Range(0, col_end_last)).copyTo
		(phase1(cv::Range(0, phase_tmp.rows), cv::Range(0, col_end_last)));
		phase2(cv::Range(0, phase2.rows), cv::Range(col_start_next, cols3)).copyTo
		(phase1(cv::Range(-row_offset, phase2.rows - row_offset), cv::Range(col_end_last, total_cols)));
	}
	else
	{
		if (temp_ret == 0)
		{
			lon_tmp(cv::Range(0, lon_tmp.rows), cv::Range(0, col_end_last)).copyTo
			(mapped_lon2(cv::Range(row_offset, lon_tmp.rows + row_offset), cv::Range(0, col_end_last)));
			mapped_lon1(cv::Range(0, mapped_lon1.rows), cv::Range(col_start_next, cols3)).copyTo
			(mapped_lon2(cv::Range(0, mapped_lon1.rows), cv::Range(col_end_last, total_cols)));

			lat_tmp(cv::Range(0, lat_tmp.rows), cv::Range(0, col_end_last)).copyTo
			(mapped_lat2(cv::Range(row_offset, lat_tmp.rows + row_offset), cv::Range(0, col_end_last)));
			mapped_lat1(cv::Range(0, mapped_lat1.rows), cv::Range(col_start_next, cols3)).copyTo
			(mapped_lat2(cv::Range(0, mapped_lat1.rows), cv::Range(col_end_last, total_cols)));
		}
		phase_tmp(cv::Range(0, phase_tmp.rows), cv::Range(0, col_end_last)).copyTo
		(phase1(cv::Range(row_offset, phase_tmp.rows + row_offset), cv::Range(0, col_end_last)));
		phase2(cv::Range(0, phase2.rows), cv::Range(col_start_next, cols3)).copyTo
		(phase1(cv::Range(0, phase2.rows), cv::Range(col_end_last, total_cols)));
	}
	ret = conversion.creat_new_h5(merged_phase_h5file);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_array_to_h5(merged_phase_h5file, "phase", phase1);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	if (temp_ret == 0)
	{
		ret = conversion.write_array_to_h5(merged_phase_h5file, "mapped_lon", mapped_lon2);
		ret = conversion.write_array_to_h5(merged_phase_h5file, "mapped_lat", mapped_lat2);
	}

	conversion.write_int_to_h5(merged_phase_h5file, "multilook_az", mul_az);
	conversion.write_int_to_h5(merged_phase_h5file, "multilook_rg", mul_rg);
	conversion.write_int_to_h5(merged_phase_h5file, "azimuth_len", phase1.rows);
	conversion.write_int_to_h5(merged_phase_h5file, "range_len", phase1.cols);

	//写入三个子带的source_1和source_2

	string source_1, source_2;
	ret = conversion.read_str_from_h5(IW1_h5file, "source_1", source_1);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_1_IW1", source_1.c_str());
	ret = conversion.read_str_from_h5(IW1_h5file, "source_2", source_2);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_2_IW1", source_2.c_str());
	ret = conversion.read_str_from_h5(IW2_h5file, "source_1", source_1);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_1_IW2", source_1.c_str());
	ret = conversion.read_str_from_h5(IW2_h5file, "source_2", source_2);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_2_IW2", source_2.c_str());
	ret = conversion.read_str_from_h5(IW3_h5file, "source_1", source_1);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_1_IW3", source_1.c_str());
	ret = conversion.read_str_from_h5(IW3_h5file, "source_1", source_1);
	ret = conversion.write_str_to_h5(merged_phase_h5file, "source_2_IW3", source_2.c_str());

	return 0;
}

int Utils::S1_subswath_merge_slc(const char* IW1_h5file, const char* IW2_h5file, const char* IW3_h5file, const char* merged_phase_h5file)
{
	if (!IW1_h5file || !IW2_h5file || !IW3_h5file || !merged_phase_h5file)
	{
		fprintf(stderr, "S1_subswath_merge_slc(): input check failed\n");
		return -1;
	}
	double start1, start2, start3, end1, end2, end3, first_pixel1, first_pixel2, first_pixel3,
		range_spacing, prf;
	int mul_az, mul_rg, mul_az1, mul_rg1, rows1, rows2, rows3, cols1, cols2, cols3;
	string sensor1, sensor2, sensor3;
	FormatConversion conversion;
	string start_time, end_time;
	int ret;
	ret = conversion.read_str_from_h5(IW1_h5file, "sensor", sensor1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW2_h5file, "sensor", sensor2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW3_h5file, "sensor", sensor3);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	if (sensor1 != "sentinel" || sensor2 != "sentinel" || sensor3 != "sentinel")
	{
		fprintf(stderr, "S1_subswath_merge_slc(): not sentinel1 data!\n");
		return -1;
	}

	ret = conversion.read_double_from_h5(IW1_h5file, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW1_h5file, "slant_range_first_pixel", &first_pixel1);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW2_h5file, "slant_range_first_pixel", &first_pixel2);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(IW3_h5file, "slant_range_first_pixel", &first_pixel3);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	if (first_pixel1 >= first_pixel2 || first_pixel2 >= first_pixel3)
	{
		fprintf(stderr, "S1_subswath_merge_slc(): please rearrange input swath order!\n");
		return -1;
	}
	ret = conversion.read_double_from_h5(IW3_h5file, "range_spacing", &range_spacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;


	ret = conversion.read_int_from_h5(IW1_h5file, "range_len", &cols1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW1_h5file, "azimuth_len", &rows1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW2_h5file, "range_len", &cols2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW2_h5file, "azimuth_len", &rows2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(IW3_h5file, "range_len", &cols3);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(IW3_h5file, "azimuth_len", &rows3);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;



	ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start1);
	conversion.utc2gps(end_time.c_str(), &end1);

	ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start2);
	conversion.utc2gps(end_time.c_str(), &end2);

	ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(start_time.c_str(), &start3);
	conversion.utc2gps(end_time.c_str(), &end3);

	//IW1和IW2拼接

	int x = (double(cols1) - round((first_pixel2 - first_pixel1) / range_spacing)) / 2.0;
	int col_end_last = x + (int)round((first_pixel2 - first_pixel1) / range_spacing);
	int col_start_next = x;
	int total_cols = col_end_last + (cols2 - col_start_next);
	int row_offset = (start1 - start2) * prf;
	int total_rows = (((end1 > end2 ? end1 : end2) - (start1 < start2 ? start1 : start2)) * prf + 1) + 1;
	ComplexMat slc1, slc2, slc3, slc_tmp, slc_tmp2;
	ret = conversion.read_slc_from_h5(IW1_h5file, slc1);
	ret = conversion.read_slc_from_h5(IW2_h5file, slc2);
	if (slc1.type() != slc2.type())
	{
		fprintf(stderr, "S1_subswath_merge_slc(): data type mismatch!\n");
		return -1;
	}
	slc_tmp.re.create(total_rows, total_cols, slc1.type()); slc_tmp.im.create(total_rows, total_cols, slc1.type());
	slc_tmp.re = 0; slc_tmp.im = 0;
	if (row_offset < 0)
	{
		slc1.re(cv::Range(0, slc1.re.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp.re(cv::Range(0, slc1.re.rows), cv::Range(0, col_end_last)));
		slc1.im(cv::Range(0, slc1.im.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp.im(cv::Range(0, slc1.im.rows), cv::Range(0, col_end_last)));

		slc2.re(cv::Range(0, slc2.re.rows), cv::Range(col_start_next, cols2)).copyTo
		(slc_tmp.re(cv::Range(-row_offset, slc2.re.rows - row_offset), cv::Range(col_end_last, total_cols)));
		slc2.im(cv::Range(0, slc2.im.rows), cv::Range(col_start_next, cols2)).copyTo
		(slc_tmp.im(cv::Range(-row_offset, slc2.im.rows - row_offset), cv::Range(col_end_last, total_cols)));
	}
	else
	{
		slc1.re(cv::Range(0, slc1.re.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp.re(cv::Range(row_offset, slc1.re.rows + row_offset), cv::Range(0, col_end_last)));
		slc1.im(cv::Range(0, slc1.im.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp.im(cv::Range(row_offset, slc1.im.rows + row_offset), cv::Range(0, col_end_last)));

		slc2.re(cv::Range(0, slc2.re.rows), cv::Range(col_start_next, cols2)).copyTo
		(slc_tmp.re(cv::Range(0, slc2.re.rows), cv::Range(col_end_last, total_cols)));
		slc2.im(cv::Range(0, slc2.im.rows), cv::Range(col_start_next, cols2)).copyTo
		(slc_tmp.im(cv::Range(0, slc2.im.rows), cv::Range(col_end_last, total_cols)));
	}

	//拼接IW3

	x = (double(total_cols) - round((first_pixel3 - first_pixel1) / range_spacing)) / 2.0;
	col_end_last = x + (int)round((first_pixel3 - first_pixel1) / range_spacing);
	col_start_next = x;
	total_cols = col_end_last + (cols3 - col_start_next);
	double start11 = start1 < start2 ? start1 : start2;
	double end11 = end1 > end2 ? end1 : end2;
	row_offset = (start11 - start3) * prf;
	total_rows = (((end11 > end3 ? end11 : end3) - (start11 < start3 ? start11 : start3)) * prf + 1) + 1;
	slc_tmp2.re.create(total_rows, total_cols, slc1.type()); slc_tmp2.im.create(total_rows, total_cols, slc1.type());
	slc_tmp2.re = 0; slc_tmp2.im = 0;
	ret = conversion.read_slc_from_h5(IW3_h5file, slc2);
	if (slc1.type() != slc2.type())
	{
		fprintf(stderr, "S1_subswath_merge_slc(): data type mismatch!\n");
		return -1;
	}
	if (row_offset < 0)
	{
		slc_tmp.re(cv::Range(0, slc_tmp.re.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp2.re(cv::Range(0, slc_tmp.re.rows), cv::Range(0, col_end_last)));
		slc_tmp.im(cv::Range(0, slc_tmp.im.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp2.im(cv::Range(0, slc_tmp.im.rows), cv::Range(0, col_end_last)));

		slc2.re(cv::Range(0, slc2.re.rows), cv::Range(col_start_next, cols3)).copyTo
		(slc_tmp2.re(cv::Range(-row_offset, slc2.re.rows - row_offset), cv::Range(col_end_last, total_cols)));
		slc2.im(cv::Range(0, slc2.im.rows), cv::Range(col_start_next, cols3)).copyTo
		(slc_tmp2.im(cv::Range(-row_offset, slc2.im.rows - row_offset), cv::Range(col_end_last, total_cols)));
	}
	else
	{
		slc_tmp.re(cv::Range(0, slc_tmp.re.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp2.re(cv::Range(row_offset, slc_tmp.re.rows + row_offset), cv::Range(0, col_end_last)));
		slc_tmp.im(cv::Range(0, slc_tmp.im.rows), cv::Range(0, col_end_last)).copyTo
		(slc_tmp2.im(cv::Range(row_offset, slc_tmp.im.rows + row_offset), cv::Range(0, col_end_last)));

		slc2.re(cv::Range(0, slc2.re.rows), cv::Range(col_start_next, cols3)).copyTo
		(slc_tmp2.re(cv::Range(0, slc2.re.rows), cv::Range(col_end_last, total_cols)));
		slc2.im(cv::Range(0, slc2.im.rows), cv::Range(col_start_next, cols3)).copyTo
		(slc_tmp2.im(cv::Range(0, slc2.im.rows), cv::Range(col_end_last, total_cols)));
	}
	ret = conversion.creat_new_h5(merged_phase_h5file);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(merged_phase_h5file, slc_tmp2);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;


	conversion.write_int_to_h5(merged_phase_h5file, "azimuth_len", slc_tmp2.re.rows);
	conversion.write_int_to_h5(merged_phase_h5file, "range_len", slc_tmp2.re.cols);

	/////写入其它参数
	//拍摄起始时间
	if (start1 < start2 && start1 < start3)
	{
		ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_start_time", start_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_start_time", start_time.c_str());
	}
	else if (start2 < start1 && start2 < start3)
	{
		ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_start_time", start_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_start_time", start_time.c_str());
	}
	else
	{
		ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_start_time", start_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_start_time", start_time.c_str());
	}
	//拍摄结束时间
	if (end1 > end2 && end1 > end3)
	{
		ret = conversion.read_str_from_h5(IW1_h5file, "acquisition_stop_time", end_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_stop_time", end_time.c_str());
	}
	else if (end2 > end1 && end2 > end3)
	{
		ret = conversion.read_str_from_h5(IW2_h5file, "acquisition_stop_time", end_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_stop_time", end_time.c_str());
	}
	else
	{
		ret = conversion.read_str_from_h5(IW3_h5file, "acquisition_stop_time", end_time);
		conversion.write_str_to_h5(merged_phase_h5file, "acquisition_stop_time", end_time.c_str());
	}
	//最近斜距
	conversion.write_double_to_h5(merged_phase_h5file, "slant_range_first_pixel", first_pixel1);

	//采样间隔
	conversion.write_double_to_h5(merged_phase_h5file, "range_spacing", range_spacing);
	conversion.read_double_from_h5(IW3_h5file, "azimuth_spacing", &range_spacing);
	conversion.write_double_to_h5(merged_phase_h5file, "azimuth_spacing", range_spacing);
	//中心下视角
	conversion.read_double_from_h5(IW2_h5file, "incidence_center", &range_spacing);
	conversion.write_double_to_h5(merged_phase_h5file, "incidence_center", range_spacing);
	//prf
	conversion.write_double_to_h5(merged_phase_h5file, "prf", prf);
	//载频
	conversion.read_double_from_h5(IW3_h5file, "carrier_frequency", &range_spacing);
	conversion.write_double_to_h5(merged_phase_h5file, "carrier_frequency", range_spacing);
	//极化、swath、传感器名
	ret = conversion.read_str_from_h5(IW3_h5file, "polarization", end_time);
	conversion.write_str_to_h5(merged_phase_h5file, "polarization", end_time.c_str());
	conversion.write_str_to_h5(merged_phase_h5file, "swath", "IW123");
	//conversion.write_str_to_h5(merged_phase_h5file, "sensor", "sentinel");
	//轨道
	Mat state_vec;
	ret = conversion.read_array_from_h5(IW3_h5file, "state_vec", state_vec);
	conversion.write_array_to_h5(merged_phase_h5file, "state_vec", state_vec);
	//控制点
	Mat gcps1, gcps2, gcps3, lon, lat;
	ret = conversion.read_array_from_h5(IW1_h5file, "gcps", gcps1);
	ret = conversion.read_array_from_h5(IW2_h5file, "gcps", gcps2);
	ret = conversion.read_array_from_h5(IW3_h5file, "gcps", gcps3);
	cv::vconcat(gcps1, gcps2, gcps2);
	cv::vconcat(gcps2, gcps3, gcps3);
	gcps3(cv::Range(0, gcps3.rows), cv::Range(0, 1)).copyTo(lon);
	gcps3(cv::Range(0, gcps3.rows), cv::Range(1, 2)).copyTo(lat);

	double topleft_lon, topright_lon, bottomleft_lon, bottomright_lon,
		topleft_lat, topright_lat, bottomleft_lat, bottomright_lat, lonMin, lonMax, latMin, latMax;

	cv::minMaxLoc(lon, &lonMin, &lonMax);
	cv::minMaxLoc(lat, &latMin, &latMax);
	double extra = 5.0 / 6000;
	lonMin = lonMin - extra * 100;
	lonMax = lonMax + extra * 100;
	latMin = latMin - extra * 100;
	latMax = latMax + extra * 100;
	topleft_lon = lonMin;
	topright_lon = lonMax;
	bottomleft_lon = lonMin;
	bottomright_lon = lonMax;
	topleft_lat = latMax;
	topright_lat = latMax;
	bottomleft_lat = latMin;
	bottomright_lat = latMin;

	conversion.write_double_to_h5(merged_phase_h5file, "topLeftLat", topleft_lat);
	conversion.write_double_to_h5(merged_phase_h5file, "topLeftLon", topleft_lon);
	conversion.write_double_to_h5(merged_phase_h5file, "topRightLat", topright_lat);
	conversion.write_double_to_h5(merged_phase_h5file, "topRightLon", topright_lon);
	conversion.write_double_to_h5(merged_phase_h5file, "bottomLeftLat", bottomleft_lat);
	conversion.write_double_to_h5(merged_phase_h5file, "bottomLeftLon", bottomleft_lon);
	conversion.write_double_to_h5(merged_phase_h5file, "bottomRightLat", bottomright_lat);
	conversion.write_double_to_h5(merged_phase_h5file, "bottomRightLon", bottomright_lon);


	return 0;
}

int Utils::S1_frame_merge(vector<string>& h5files, const char* merged_phase_h5)
{
	if (h5files.size() < 2 || !merged_phase_h5)
	{
		fprintf(stderr, "S1_frame_merge(): input check failed!\n");
		return -1;
	}
	Mat phase;
	int num_files = h5files.size();
	//根据每个拍摄frame的拍摄起始时间对文件排序（从小到大）
	int ret; string start_time; double start;
	Mat stime(1, num_files, CV_64F), order;
	FormatConversion conversion;
	for (int i = 0; i < num_files; i++)
	{
		ret = conversion.read_str_from_h5(h5files[i].c_str(), "acquisition_start_time", start_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time.c_str(), &start);
		if (return_check(ret, "utc2gps()", error_head)) return -1;
		stime.at<double>(0, i) = start;
	}
	cv::sortIdx(stime, order, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	//检查多视倍数和最近斜距是否相同
	double start1, start2, start3, end1, end2, end3, first_pixel1, first_pixel2, first_pixel3,
		range_spacing, prf;
	int mul_az, mul_rg, mul_az1, mul_rg1;
	ret = conversion.read_int_from_h5(h5files[0].c_str(), "multilook_az", &mul_az);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(h5files[0].c_str(), "multilook_rg", &mul_rg);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5files[0].c_str(), "slant_range_first_pixel", &first_pixel1);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5files[0].c_str(), "range_spacing", &range_spacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	for (int i = 1; i < num_files; i++)
	{
		ret = conversion.read_int_from_h5(h5files[i].c_str(), "multilook_az", &mul_az1);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		ret = conversion.read_int_from_h5(h5files[i].c_str(), "multilook_rg", &mul_rg1);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		if (mul_az != mul_az1 || mul_rg != mul_rg1)
		{
			fprintf(stderr, "S1_frame_merge(): multilook times disagree!\n");
			return -1;
		}
		ret = conversion.read_double_from_h5(h5files[i].c_str(), "slant_range_first_pixel", &first_pixel2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		if (fabs(first_pixel2 - first_pixel1) > 0.1)
		{
			fprintf(stderr, "S1_frame_merge(): slant_range_first_pixel disagree!\n");
			return -1;
		}
	}

	Mat phase1;
	ret = conversion.creat_new_h5(merged_phase_h5);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5files[order.at<int>(0)].c_str(), "phase", phase);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(h5files[order.at<int>(0)].c_str(), "acquisition_stop_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &end1);
	if (return_check(ret, "utc2gps()", error_head)) return -1;
	ret = conversion.read_str_from_h5(h5files[order.at<int>(0)].c_str(), "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	if (return_check(ret, "utc2gps()", error_head)) return -1;
	conversion.write_str_to_h5(merged_phase_h5, "acquisition_start_time", start_time.c_str());
	ret = conversion.read_str_from_h5(h5files[order.at<int>(num_files - 1)].c_str(), "acquisition_stop_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.write_str_to_h5(merged_phase_h5, "acquisition_stop_time", start_time.c_str());
	ret = conversion.read_double_from_h5(h5files[order.at<int>(0)].c_str(), "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	for (int i = 1; i < num_files; i++)
	{
		ret = conversion.read_array_from_h5(h5files[order.at<int>(i)].c_str(), "phase", phase1);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		if (phase1.cols != phase.cols)
		{
			fprintf(stderr, "S1_frame_merge(): frame cols mismatch!\n");
			return -1;
		}
		ret = conversion.read_str_from_h5(h5files[order.at<int>(i)].c_str(), "acquisition_start_time", start_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time.c_str(), &start2);
		if (return_check(ret, "utc2gps()", error_head)) return -1;
		if (end1 < start2)
		{
			fprintf(stderr, "S1_frame_merge(): no overlap area between frames!\n");
			return -1;
		}
		int last_upper_count = ((start2 - start) * prf + 1 + ((end1 - start2) * prf + 1) / 2.0) / (double)mul_az;
		int tmp = last_upper_count * mul_az - ((start2 - start) * prf + 1);
		tmp = (end1 - start2) * prf + 1 - tmp;
		int next_lower_start = round((double)tmp / (double)mul_az);

		phase(cv::Range(0, last_upper_count), cv::Range(0, phase.cols)).copyTo(phase);
		phase1(cv::Range(next_lower_start, phase1.rows), cv::Range(0, phase1.cols)).copyTo(phase1);
		cv::vconcat(phase, phase1, phase);

		ret = conversion.read_str_from_h5(h5files[order.at<int>(i)].c_str(), "acquisition_stop_time", start_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time.c_str(), &end1);
		if (return_check(ret, "utc2gps()", error_head)) return -1;
	}
	conversion.write_array_to_h5(merged_phase_h5, "phase", phase);
	conversion.write_int_to_h5(merged_phase_h5, "multilook_az", mul_az);
	conversion.write_int_to_h5(merged_phase_h5, "multilook_rg", mul_rg);
	conversion.write_int_to_h5(merged_phase_h5, "azimuth_len", phase.rows);
	conversion.write_int_to_h5(merged_phase_h5, "range_len", phase.cols);
	conversion.write_double_to_h5(merged_phase_h5, "prf", prf);
	conversion.write_double_to_h5(merged_phase_h5, "range_spacing", range_spacing);
	conversion.write_double_to_h5(merged_phase_h5, "slant_range_first_pixel", first_pixel1);
	return 0;
}

int Utils::S1_frame_merge(const char* frame1_h5, const char* frame2_h5, const char* outframe_h5)
{
	if (!frame1_h5 || !frame2_h5 || !outframe_h5)
	{
		fprintf(stderr, "S1_frame_merge(): input check failed!\n");
		return -1;
	}
	//检查是否属于同一轨道相邻frame
	int ret;
	FormatConversion conversion;
	double start1, start2, end1, end2, prf, slant_range_first_pixel1, slant_range_first_pixel2;
	string start_time1, end_time1, start_time2, end_time2;
	ret = conversion.read_str_from_h5(frame1_h5, "acquisition_start_time", start_time1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time1.c_str(), &start1);
	ret = conversion.read_str_from_h5(frame1_h5, "acquisition_stop_time", end_time1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(end_time1.c_str(), &end1);
	ret = conversion.read_str_from_h5(frame2_h5, "acquisition_start_time", start_time2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time2.c_str(), &start2);
	ret = conversion.read_str_from_h5(frame2_h5, "acquisition_stop_time", end_time2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(end_time2.c_str(), &end2);
	if (start1 < start2)
	{
		if (start2 > end1)
		{
			fprintf(stderr, "S1_frame_merge(): not adjacent frames!\n");
			return -1;
		}
	}
	else
	{
		if (start1 > end2)
		{
			fprintf(stderr, "S1_frame_merge(): not adjacent frames!\n");
			return -1;
		}
	}
	ret = conversion.read_double_from_h5(frame1_h5, "slant_range_first_pixel", &slant_range_first_pixel1);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(frame2_h5, "slant_range_first_pixel", &slant_range_first_pixel2);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	if (fabs(slant_range_first_pixel1 - slant_range_first_pixel2) > 0.1)
	{
		fprintf(stderr, "S1_frame_merge(): not the same track!\n");
		return -1;
	}

	ret = conversion.creat_new_h5(outframe_h5);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//融合azimuthFmRateList
	Mat azimuthFmRateList1, azimuthFmRateList2;
	ret = conversion.read_array_from_h5(frame1_h5, "azimuthFmRateList", azimuthFmRateList1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "azimuthFmRateList", azimuthFmRateList2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		int t_remain = azimuthFmRateList1.rows - 1;
		for (int i = 0; i < azimuthFmRateList1.rows; i++)
		{
			if (azimuthFmRateList2.at<double>(0, 0) <= azimuthFmRateList1.at<double>(i, 0))
			{
				t_remain = i; break;
			}
		}
		azimuthFmRateList1(cv::Range(0, t_remain), cv::Range(0, azimuthFmRateList1.cols)).copyTo(azimuthFmRateList1);
		cv::vconcat(azimuthFmRateList1, azimuthFmRateList2, azimuthFmRateList1);
		conversion.write_array_to_h5(outframe_h5, "azimuthFmRateList", azimuthFmRateList1);
	}
	else
	{
		int t_remain = azimuthFmRateList2.rows - 1;
		for (int i = 0; i < azimuthFmRateList2.rows; i++)
		{
			if ((azimuthFmRateList1.at<double>(0, 0) <= azimuthFmRateList2.at<double>(i, 0)))
			{
				t_remain = i; break;
			}
		}
		azimuthFmRateList2(cv::Range(0, t_remain), cv::Range(0, azimuthFmRateList2.cols)).copyTo(azimuthFmRateList2);
		cv::vconcat(azimuthFmRateList2, azimuthFmRateList1, azimuthFmRateList1);
		conversion.write_array_to_h5(outframe_h5, "azimuthFmRateList", azimuthFmRateList1);
	}
	//融合dcEstimateList
	Mat dcEstimateList1, dcEstimateList2;
	ret = conversion.read_array_from_h5(frame1_h5, "dcEstimateList", dcEstimateList1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "dcEstimateList", dcEstimateList2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		int t_remain = dcEstimateList1.rows - 1;
		for (int i = 0; i < dcEstimateList1.rows; i++)
		{
			if ((dcEstimateList2.at<double>(0, 0) <= dcEstimateList1.at<double>(i, 0)))
			{
				t_remain = i; break;
			}
		}
		dcEstimateList1(cv::Range(0, t_remain), cv::Range(0, dcEstimateList1.cols)).copyTo(dcEstimateList1);
		cv::vconcat(dcEstimateList1, dcEstimateList2, dcEstimateList1);
		conversion.write_array_to_h5(outframe_h5, "dcEstimateList", dcEstimateList1);
	}
	else
	{
		int t_remain = dcEstimateList2.rows - 1;
		for (int i = 0; i < dcEstimateList2.rows; i++)
		{
			if ((dcEstimateList1.at<double>(0, 0) <= dcEstimateList2.at<double>(i, 0)))
			{
				t_remain = i; break;
			}
		}
		dcEstimateList2(cv::Range(0, t_remain), cv::Range(0, dcEstimateList2.cols)).copyTo(dcEstimateList2);
		cv::vconcat(dcEstimateList2, dcEstimateList1, dcEstimateList1);
		conversion.write_array_to_h5(outframe_h5, "dcEstimateList", dcEstimateList1);
	}
	//融合state_vec
	Mat state_vec1, state_vec2;
	ret = conversion.read_array_from_h5(frame1_h5, "state_vec", state_vec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "state_vec", state_vec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		int t_remain = state_vec1.rows - 1;
		for (int i = 0; i < state_vec1.rows; i++)
		{
			if (state_vec2.at<double>(0, 0) <= state_vec1.at<double>(i, 0))
			{
				t_remain = i; break;
			}
		}
		if (t_remain != 0)
		{
			state_vec1(cv::Range(0, t_remain), cv::Range(0, state_vec1.cols)).copyTo(state_vec1);
			cv::vconcat(state_vec1, state_vec2, state_vec1);
		}
		else
		{
			state_vec2.copyTo(state_vec1);
		}
		conversion.write_array_to_h5(outframe_h5, "state_vec", state_vec1);
	}
	else
	{
		int t_remain = state_vec2.rows - 1;
		for (int i = 0; i < state_vec2.rows; i++)
		{
			if (state_vec1.at<double>(0, 0) <= state_vec2.at<double>(i, 0))
			{
				t_remain = i; break;
			}
		}
		if (t_remain != 0)
		{
			state_vec2(cv::Range(0, t_remain), cv::Range(0, state_vec2.cols)).copyTo(state_vec2);
			cv::vconcat(state_vec2, state_vec1, state_vec1);
		}
		else
		{

		}
		
		conversion.write_array_to_h5(outframe_h5, "state_vec", state_vec1);
	}
	//融合fine_state_vec
	Mat fine_state_vec1, fine_state_vec2;
	ret = conversion.read_array_from_h5(frame1_h5, "fine_state_vec", fine_state_vec1);
	ret += conversion.read_array_from_h5(frame2_h5, "fine_state_vec", fine_state_vec2);
	if (ret == 0)
	{
		if (start1 < start2)
		{
			int t_remain = fine_state_vec1.rows - 1;
			for (int i = 0; i < fine_state_vec1.rows; i++)
			{
				if (fine_state_vec2.at<double>(0, 0) <= fine_state_vec1.at<double>(i, 0))
				{
					t_remain = i; break;
				}
			}
			if (t_remain != 0)
			{
				fine_state_vec1(cv::Range(0, t_remain), cv::Range(0, fine_state_vec1.cols)).copyTo(fine_state_vec1);
				cv::vconcat(fine_state_vec1, fine_state_vec2, fine_state_vec1);
			}
			else
			{
				fine_state_vec2.copyTo(fine_state_vec1);
			}
			conversion.write_array_to_h5(outframe_h5, "fine_state_vec", fine_state_vec1);
		}
		else
		{
			int t_remain = fine_state_vec2.rows - 1;
			for (int i = 0; i < fine_state_vec2.rows; i++)
			{
				if (fine_state_vec1.at<double>(0, 0) <= fine_state_vec2.at<double>(i, 0))
				{
					t_remain = i; break;
				}
			}
			if (t_remain != 0)
			{
				fine_state_vec2(cv::Range(0, t_remain), cv::Range(0, fine_state_vec2.cols)).copyTo(fine_state_vec2);
				cv::vconcat(fine_state_vec2, fine_state_vec1, fine_state_vec1);
			}
			
			conversion.write_array_to_h5(outframe_h5, "fine_state_vec", fine_state_vec1);
		}
	}
	
	//融合gcps
	Mat gcps1, gcps2;
	int rows1, rows2;
	ret = conversion.read_int_from_h5(frame1_h5, "azimuth_len", &rows1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(frame2_h5, "azimuth_len", &rows2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame1_h5, "gcps", gcps1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "gcps", gcps2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		int t;
		for (int i = 0; i < gcps2.rows; i++)
		{
			if (fabs(gcps2.at<double>(i, 3)) > 1.0)
			{
				t = i; break;
			}
		}
		gcps2(cv::Range(t, gcps2.rows), cv::Range(0, gcps2.cols)).copyTo(gcps2);
		Mat temp = gcps2(cv::Range(0, gcps2.rows), cv::Range(3, 4)) + gcps1.at<double>(gcps1.rows - 1, 3);
		temp.copyTo(gcps2(cv::Range(0, gcps2.rows), cv::Range(3, 4)));
		cv::vconcat(gcps1, gcps2, gcps1);
		conversion.write_array_to_h5(outframe_h5, "gcps", gcps1);
	}
	else
	{
		int t;
		for (int i = 0; i < gcps1.rows; i++)
		{
			if (fabs(gcps1.at<double>(i, 3)) > 1.0)
			{
				t = i; break;
			}
		}
		gcps1(cv::Range(t, gcps1.rows), cv::Range(0, gcps1.cols)).copyTo(gcps1);
		Mat temp = gcps1(cv::Range(0, gcps1.rows), cv::Range(3, 4)) + gcps2.at<double>(gcps2.rows - 1, 3);
		temp.copyTo(gcps1(cv::Range(0, gcps1.rows), cv::Range(3, 4)));
		cv::vconcat(gcps2, gcps1, gcps1);
		conversion.write_array_to_h5(outframe_h5, "gcps", gcps1);
	}
	//融合burstAzimuthTime
	Mat burstAzimuthTime1, burstAzimuthTime2;
	ret = conversion.read_array_from_h5(frame1_h5, "burstAzimuthTime", burstAzimuthTime1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "burstAzimuthTime", burstAzimuthTime2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		cv::vconcat(burstAzimuthTime1, burstAzimuthTime2, burstAzimuthTime1);
	}
	else
	{
		cv::vconcat(burstAzimuthTime2, burstAzimuthTime1, burstAzimuthTime1);
	}
	conversion.write_array_to_h5(outframe_h5, "burstAzimuthTime", burstAzimuthTime1);
	//融合firstValidLine
	Mat firstValidLine1, firstValidLine2;
	ret = conversion.read_array_from_h5(frame1_h5, "firstValidLine", firstValidLine1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "firstValidLine", firstValidLine2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		cv::vconcat(firstValidLine1, firstValidLine2, firstValidLine1);
	}
	else
	{
		cv::vconcat(firstValidLine2, firstValidLine1, firstValidLine1);
	}
	conversion.write_array_to_h5(outframe_h5, "firstValidLine", firstValidLine1);
	//融合firstValidSample
	Mat firstValidSample1, firstValidSample2;
	ret = conversion.read_array_from_h5(frame1_h5, "firstValidSample", firstValidSample1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "firstValidSample", firstValidSample2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		cv::vconcat(firstValidSample1, firstValidSample2, firstValidSample1);
	}
	else
	{
		cv::vconcat(firstValidSample2, firstValidSample1, firstValidSample1);
	}
	conversion.write_array_to_h5(outframe_h5, "firstValidSample", firstValidSample1);
	//融合lastValidLine
	Mat lastValidLine1, lastValidLine2;
	ret = conversion.read_array_from_h5(frame1_h5, "lastValidLine", lastValidLine1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "lastValidLine", lastValidLine2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		cv::vconcat(lastValidLine1, lastValidLine2, lastValidLine1);
	}
	else
	{
		cv::vconcat(lastValidLine2, lastValidLine1, lastValidLine1);
	}
	conversion.write_array_to_h5(outframe_h5, "lastValidLine", lastValidLine1);
	//融合lastValidSample
	Mat lastValidSample1, lastValidSample2;
	ret = conversion.read_array_from_h5(frame1_h5, "lastValidSample", lastValidSample1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(frame2_h5, "lastValidSample", lastValidSample2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (start1 < start2)
	{
		cv::vconcat(lastValidSample1, lastValidSample2, lastValidSample1);
	}
	else
	{
		cv::vconcat(lastValidSample2, lastValidSample1, lastValidSample1);
	}
	conversion.write_array_to_h5(outframe_h5, "lastValidSample", lastValidSample1);
	//融合拍摄时间
	if (start1 < start2)
	{
		conversion.write_str_to_h5(outframe_h5, "acquisition_start_time", start_time1.c_str());
		conversion.write_str_to_h5(outframe_h5, "acquisition_stop_time", end_time2.c_str());
	}
	else
	{
		conversion.write_str_to_h5(outframe_h5, "acquisition_start_time", start_time2.c_str());
		conversion.write_str_to_h5(outframe_h5, "acquisition_stop_time", end_time1.c_str());
	}
	//融合azimuthSteeringRate
	double azimuthSteeringRate;
	ret = conversion.read_double_from_h5(frame1_h5, "azimuthSteeringRate", &azimuthSteeringRate);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	conversion.write_double_to_h5(outframe_h5, "azimuthSteeringRate", azimuthSteeringRate);
	//融合azimuth_len，range_len
	int azimuth_len1, range_len1, azimuth_len2, range_len2;
	ret = conversion.read_int_from_h5(frame1_h5, "azimuth_len", &azimuth_len1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(frame1_h5, "range_len", &range_len1);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(frame2_h5, "azimuth_len", &azimuth_len2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(frame2_h5, "range_len", &range_len2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	azimuth_len1 += azimuth_len2;
	conversion.write_int_to_h5(outframe_h5, "azimuth_len", azimuth_len1);
	int range_len = range_len1 >= range_len2 ? range_len1 : range_len2;
	conversion.write_int_to_h5(outframe_h5, "range_len", range_len1);
	//融合azimuth_spacing，range_spacing
	double azimuth_spacing, range_spacing;
	ret = conversion.read_double_from_h5(frame1_h5, "azimuth_spacing", &azimuth_spacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	conversion.write_double_to_h5(outframe_h5, "azimuth_spacing", azimuth_spacing);
	ret = conversion.read_double_from_h5(frame1_h5, "range_spacing", &range_spacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	conversion.write_double_to_h5(outframe_h5, "range_spacing", range_spacing);
	//burstCount, carrier_frequency, heading, incidence_center, linesPerBurst, prf, samplesPerBurst, 
	double carrier_frequency, heading, incidence_center, slant_range_first_pixel;
	int burstCount1, burstCount2, linesPerBurst, samplesPerBurst;
	conversion.read_double_from_h5(frame1_h5, "carrier_frequency", &carrier_frequency);
	conversion.read_double_from_h5(frame1_h5, "slant_range_first_pixel", &slant_range_first_pixel);
	conversion.read_double_from_h5(frame1_h5, "incidence_center", &incidence_center);
	conversion.read_double_from_h5(frame1_h5, "heading", &heading);
	conversion.read_double_from_h5(frame1_h5, "prf", &prf);
	conversion.read_int_from_h5(frame1_h5, "burstCount", &burstCount1);
	conversion.read_int_from_h5(frame2_h5, "burstCount", &burstCount2);
	conversion.read_int_from_h5(frame1_h5, "linesPerBurst", &linesPerBurst);
	//conversion.read_int_from_h5(frame1_h5, "samplesPerBurst", &samplesPerBurst);
	conversion.write_double_to_h5(outframe_h5, "slant_range_first_pixel", slant_range_first_pixel);
	conversion.write_double_to_h5(outframe_h5, "carrier_frequency", carrier_frequency);
	conversion.write_double_to_h5(outframe_h5, "incidence_center", incidence_center);
	conversion.write_double_to_h5(outframe_h5, "heading", heading);
	conversion.write_double_to_h5(outframe_h5, "prf", prf);
	conversion.write_int_to_h5(outframe_h5, "burstCount", burstCount1 + burstCount2);
	conversion.write_int_to_h5(outframe_h5, "linesPerBurst", linesPerBurst);
	conversion.write_int_to_h5(outframe_h5, "samplesPerBurst", range_len);

	//sensor, swath, polarization, orbit_dir, imaging_mode;
	string sensor, swath, polarization, orbit_dir, imaging_mode;
	conversion.read_str_from_h5(frame1_h5, "sensor", sensor);
	conversion.read_str_from_h5(frame1_h5, "swath", swath);
	conversion.read_str_from_h5(frame1_h5, "polarization", polarization);
	conversion.read_str_from_h5(frame1_h5, "orbit_dir", orbit_dir);
	conversion.read_str_from_h5(frame1_h5, "imaging_mode", imaging_mode);

	conversion.write_str_to_h5(outframe_h5, "sensor", sensor.c_str());
	conversion.write_str_to_h5(outframe_h5, "swath", swath.c_str());
	conversion.write_str_to_h5(outframe_h5, "polarization", polarization.c_str());
	conversion.write_str_to_h5(outframe_h5, "orbit_dir", orbit_dir.c_str());
	conversion.write_str_to_h5(outframe_h5, "imaging_mode", imaging_mode.c_str());

	//s_re, s_im
	Mat s_re, s_re2;
	conversion.read_array_from_h5(frame1_h5, "s_re", s_re);
	conversion.read_array_from_h5(frame2_h5, "s_re", s_re2);
	if (range_len1 != range_len2)
	{
		if (range_len1 > range_len2)
		{
			cv::copyMakeBorder(s_re2, s_re2, 0, 0, 0, range_len1 - range_len2, BORDER_CONSTANT, cv::Scalar(0));
		}
		else
		{
			cv::copyMakeBorder(s_re, s_re, 0, 0, 0, range_len2 - range_len1, BORDER_CONSTANT, cv::Scalar(0));
		}
	}
	if (start1 < start2)
	{
		cv::vconcat(s_re, s_re2, s_re);
	}
	else
	{
		cv::vconcat(s_re2, s_re, s_re);
	}
	conversion.write_array_to_h5(outframe_h5, "s_re", s_re);

	conversion.read_array_from_h5(frame1_h5, "s_im", s_re);
	conversion.read_array_from_h5(frame2_h5, "s_im", s_re2);
	if (range_len1 != range_len2)
	{
		if (range_len1 > range_len2)
		{
			cv::copyMakeBorder(s_re2, s_re2, 0, 0, 0, range_len1 - range_len2, BORDER_CONSTANT, cv::Scalar(0));
		}
		else
		{
			cv::copyMakeBorder(s_re, s_re, 0, 0, 0, range_len2 - range_len1, BORDER_CONSTANT, cv::Scalar(0));
		}
	}
	if (start1 < start2)
	{
		cv::vconcat(s_re, s_re2, s_re);
	}
	else
	{
		cv::vconcat(s_re2, s_re, s_re);
	}
	conversion.write_array_to_h5(outframe_h5, "s_im", s_re);


	//拟合系数

	{
		Mat lon_coefficient, lat_coefficient, inc_coefficient, row_coefficient, col_coefficient;
		int numberOfSamples = range_len;
		double mean_lon, mean_lat, mean_inc, max_lon, max_lat, max_inc, min_lon, min_lat, min_inc;
		Mat lon, lat, inc, row, col, gcps;
		gcps1.copyTo(gcps);
		gcps(cv::Range(0, gcps.rows), cv::Range(0, 1)).copyTo(lon);
		gcps(cv::Range(0, gcps.rows), cv::Range(1, 2)).copyTo(lat);
		gcps(cv::Range(0, gcps.rows), cv::Range(3, 4)).copyTo(row);
		gcps(cv::Range(0, gcps.rows), cv::Range(4, 5)).copyTo(col);
		gcps(cv::Range(0, gcps.rows), cv::Range(5, 6)).copyTo(inc);
		mean_lon = cv::mean(lon)[0];
		mean_lat = cv::mean(lat)[0];
		mean_inc = cv::mean(inc)[0];
		cv::minMaxLoc(lon, &min_lon, &max_lon);
		cv::minMaxLoc(lat, &min_lat, &max_lat);
		cv::minMaxLoc(inc, &min_inc, &max_inc);
		lon = (lon - mean_lon) / (max_lon - min_lon + 1e-10);
		lat = (lat - mean_lat) / (max_lat - min_lat + 1e-10);
		inc = (inc - mean_inc) / (max_inc - min_inc + 1e-10);
		row = (row + 1 - double(numberOfSamples) * 0.5) / (double(numberOfSamples) + 1e-10);//sentinel行列起点为0，+1统一为1.
		col = (col + 1 - double(numberOfSamples) * 0.5) / (double(numberOfSamples) + 1e-10);

		//拟合经度

		Mat A, B, b, temp, coefficient, error, eye, b_t, a, a_t;
		double rms;
		lon.copyTo(b);
		A = Mat::ones(lon.rows, 25, CV_64F);
		row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
		temp = row.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

		col.copyTo(temp);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

		cv::transpose(A, temp);
		B = temp * b;
		A.copyTo(a);
		cv::transpose(a, a_t);
		A = temp * A;
		rms = -1.0;
		if (cv::invert(A, error, cv::DECOMP_LU) > 0)
		{
			cv::transpose(b, b_t);
			error = b_t * b - (b_t * a) * error * (a_t * b);
			//error = b_t * (eye - a * error * a_t) * b;
			rms = sqrt(error.at<double>(0, 0) / double(b.rows));
		}
		if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
		{
			temp.create(1, 32, CV_64F);
			temp.at<double>(0, 0) = mean_lon;
			temp.at<double>(0, 1) = max_lon - min_lon + 1e-10;
			temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 4) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 5) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 31) = rms;
			cv::transpose(coefficient, coefficient);
			coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
			temp.copyTo(lon_coefficient);
		}

		//拟合纬度

		lat.copyTo(b);

		A = Mat::ones(lon.rows, 25, CV_64F);
		row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
		temp = row.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

		col.copyTo(temp);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

		col.copyTo(temp);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
		temp = temp.mul(row);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

		cv::transpose(A, temp);
		B = temp * b;
		A.copyTo(a);
		cv::transpose(a, a_t);
		A = temp * A;
		rms = -1.0;
		if (cv::invert(A, error, cv::DECOMP_LU) > 0)
		{
			cv::transpose(b, b_t);
			error = b_t * b - (b_t * a) * error * (a_t * b);
			//error = b_t * (eye - a * error * a_t) * b;
			rms = sqrt(error.at<double>(0, 0) / double(b.rows));
		}
		if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
		{
			temp.create(1, 32, CV_64F);
			temp.at<double>(0, 0) = mean_lat;
			temp.at<double>(0, 1) = max_lat - min_lat + 1e-10;
			temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 4) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 5) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 31) = rms;
			cv::transpose(coefficient, coefficient);
			coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
			temp.copyTo(lat_coefficient);
		}

		//拟合下视角

		inc.copyTo(b);
		A = Mat::ones(inc.rows, 6, CV_64F);
		col.copyTo(A(cv::Range(0, inc.rows), cv::Range(1, 2)));
		temp = col.mul(col);
		temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(2, 3)));
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(3, 4)));
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(4, 5)));
		temp = temp.mul(col);
		temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(5, 6)));
		cv::transpose(A, temp);
		B = temp * b;
		A.copyTo(a);
		cv::transpose(a, a_t);
		A = temp * A;
		rms = -1.0;
		if (cv::invert(A, error, cv::DECOMP_LU) > 0)
		{
			cv::transpose(b, b_t);
			error = b_t * b - (b_t * a) * error * (a_t * b);
			//error = b_t * (eye - a * error * a_t) * b;
			rms = sqrt(error.at<double>(0, 0) / double(b.rows));
		}
		if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
		{
			temp.create(1, 11, CV_64F);
			temp.at<double>(0, 0) = mean_inc;
			temp.at<double>(0, 1) = max_inc - min_inc + 1e-10;
			temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 10) = rms;
			cv::transpose(coefficient, coefficient);
			coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(4, 10)));
			temp.copyTo(inc_coefficient);
		}

		//拟合行坐标

		row.copyTo(b);
		A = Mat::ones(lon.rows, 25, CV_64F);
		lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
		temp = lon.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

		lat.copyTo(temp);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

		cv::transpose(A, temp);
		B = temp * b;
		A.copyTo(a);
		cv::transpose(a, a_t);
		A = temp * A;
		rms = -1.0;
		if (cv::invert(A, error, cv::DECOMP_LU) > 0)
		{
			cv::transpose(b, b_t);
			error = b_t * b - (b_t * a) * error * (a_t * b);
			//error = b_t * (eye - a * error * a_t) * b;
			rms = sqrt(error.at<double>(0, 0) / double(b.rows));
		}
		if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
		{
			temp.create(1, 32, CV_64F);
			temp.at<double>(0, 0) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 1) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 2) = mean_lon;
			temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
			temp.at<double>(0, 4) = mean_lat;
			temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
			temp.at<double>(0, 31) = rms;
			cv::transpose(coefficient, coefficient);
			coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
			temp.copyTo(row_coefficient);
		}

		//拟合列坐标

		col.copyTo(b);
		A = Mat::ones(lon.rows, 25, CV_64F);
		lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
		temp = lon.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

		lat.copyTo(temp);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

		lat.copyTo(temp);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp = temp.mul(lat);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
		temp = temp.mul(lon);
		temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

		cv::transpose(A, temp);
		B = temp * b;
		A.copyTo(a);
		cv::transpose(a, a_t);
		A = temp * A;
		rms = -1.0;
		if (cv::invert(A, error, cv::DECOMP_LU) > 0)
		{
			cv::transpose(b, b_t);
			error = b_t * b - (b_t * a) * error * (a_t * b);
			//error = b_t * (eye - a * error * a_t) * b;
			rms = sqrt(error.at<double>(0, 0) / double(b.rows));
		}
		if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
		{
			temp.create(1, 32, CV_64F);
			temp.at<double>(0, 0) = double(numberOfSamples) * 0.5;
			temp.at<double>(0, 1) = double(numberOfSamples) + 1e-10;
			temp.at<double>(0, 2) = mean_lon;
			temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
			temp.at<double>(0, 4) = mean_lat;
			temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
			temp.at<double>(0, 31) = rms;
			cv::transpose(coefficient, coefficient);
			coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
			temp.copyTo(col_coefficient);
		}
		conversion.write_array_to_h5(outframe_h5, "lon_coefficient", lon_coefficient);
		conversion.write_array_to_h5(outframe_h5, "lat_coefficient", lat_coefficient);
		conversion.write_array_to_h5(outframe_h5, "inc_coefficient", inc_coefficient);
		conversion.write_array_to_h5(outframe_h5, "row_coefficient", row_coefficient);
		conversion.write_array_to_h5(outframe_h5, "col_coefficient", col_coefficient);
	}
	
	return 0;
}

int Utils::SAR2UTM(
	Mat& mapped_lon,
	Mat& mapped_lat,
	Mat& phase, 
	Mat& mapped_phase,
	int interpolation_method,
	double* lon_east,
	double* lon_west,
	double* lat_north,
	double* lat_south
)
{
	if (mapped_lat.size() != mapped_lon.size() ||
		mapped_lat.size() != phase.size() ||
		phase.rows < 2 ||
		phase.cols < 2 ||
		phase.type() != CV_64F ||
		mapped_lat.type() != mapped_lon.type() ||
		(mapped_lon.type() != CV_32F && mapped_lon.type() != CV_64F)
		)
	{
		fprintf(stderr, "SAR2UTM(): input check failed!\n");
		return -1;
	}
	//确定经纬度覆盖范围
	double max_lon = -380.0, min_lon = 380.0, max_lat = -380.0, min_lat = 180.0;
	if (mapped_lat.type() == CV_32F)
	{
		for (int i = 0; i < mapped_lat.rows; i++)
		{
			for (int j = 0; j < mapped_lat.cols; j++)
			{
				if (mapped_lat.at<float>(i, j) < 350.0)
				{
					min_lat = min_lat > mapped_lat.at<float>(i, j) ? mapped_lat.at<float>(i, j) : min_lat;
					max_lat = max_lat < mapped_lat.at<float>(i, j) ? mapped_lat.at<float>(i, j) : max_lat;
					min_lon = min_lon > mapped_lon.at<float>(i, j) ? mapped_lon.at<float>(i, j) : min_lon;
					max_lon = max_lon < mapped_lon.at<float>(i, j) ? mapped_lon.at<float>(i, j) : max_lon;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < mapped_lat.rows; i++)
		{
			for (int j = 0; j < mapped_lat.cols; j++)
			{
				if (mapped_lat.at<double>(i, j) < 350.0)
				{
					min_lat = min_lat > mapped_lat.at<double>(i, j) ? mapped_lat.at<double>(i, j) : min_lat;
					max_lat = max_lat < mapped_lat.at<double>(i, j) ? mapped_lat.at<double>(i, j) : max_lat;
					min_lon = min_lon > mapped_lon.at<double>(i, j) ? mapped_lon.at<double>(i, j) : min_lon;
					max_lon = max_lon < mapped_lon.at<double>(i, j) ? mapped_lon.at<double>(i, j) : max_lon;
				}
			}
		}
	}
	//cv::minMaxLoc(mapped_lon, &min_lon, &max_lon);
	//cv::minMaxLoc(mapped_lat, &min_lat, &max_lat);
	double west = max_lon - min_lon > 180.0 ? max_lon : min_lon;
	if (lon_west) *lon_west = west;
	
	//确定经纬度采样间隔
	double lon_interval, lat_interval;
	Mat temp1, temp2;
	int rows_start = mapped_lon.rows / 4;
	int rows_end = rows_start + mapped_lon.rows / 4;
	mapped_lon(cv::Range(rows_start, rows_end), cv::Range(0, mapped_lon.cols)).copyTo(temp1);
	mapped_lon(cv::Range(rows_start + 1, rows_end + 1), cv::Range(0, mapped_lon.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = temp1 / 180.0 * PI;
	wrap(temp1, temp1);
	temp1 = temp1 / PI * 180.0;
	temp1 = cv::abs(temp1);
	lon_interval = cv::mean(temp1)[0];
	double sigma1, sigma2;
	std(temp1, &sigma1);

	mapped_lon(cv::Range(rows_start, rows_end), cv::Range(0, mapped_lon.cols - 1)).copyTo(temp1);
	mapped_lon(cv::Range(rows_start, rows_end), cv::Range(1, mapped_lon.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = temp1 / 180.0 * PI;
	wrap(temp1, temp1);
	temp1 = temp1 / PI * 180.0;
	temp1 = cv::abs(temp1);
	std(temp1, &sigma2);
	lon_interval = (cv::mean(temp1)[0] + fabs(lon_interval)) / 2.0;
	lon_interval = lon_interval + 1.0 * (sigma1 + sigma2) / 2.0;

	mapped_lat(cv::Range(rows_start, rows_end), cv::Range(0, mapped_lat.cols)).copyTo(temp1);
	mapped_lat(cv::Range(rows_start + 1, rows_end + 1), cv::Range(0, mapped_lat.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = cv::abs(temp1);
	lat_interval = cv::mean(temp1)[0];
	std(temp1, &sigma1);

	mapped_lat(cv::Range(rows_start, rows_end), cv::Range(0, mapped_lat.cols - 1)).copyTo(temp1);
	mapped_lat(cv::Range(rows_start, rows_end), cv::Range(1, mapped_lat.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = cv::abs(temp1);
	lat_interval = (cv::mean(temp1)[0] + lat_interval) / 2.0;
	std(temp1, &sigma2);
	lat_interval = lat_interval + 1.0 * (sigma1 + sigma2) / 2.0;

	int rows = phase.rows; int cols = phase.cols;
	

	//计算UTM坐标系相位尺寸
	int UTM_rows = (max_lat - min_lat) / lat_interval;
	UTM_rows += 2;
	double south = max_lat - (double)(UTM_rows - 1) * lat_interval;
	double north = max_lat;
	if (lat_north) *lat_north = north;
	if (lat_south) *lat_south = south;
	double max_lon_temp = max_lon - min_lon;
	max_lon_temp = max_lon_temp > 180.0 ? 360.0 - max_lon_temp : max_lon_temp;
	int UTM_cols = max_lon_temp / lon_interval;
	UTM_cols += 2;
	double east = west + (double)(UTM_cols - 1) * lon_interval;
	east = east > 180.0 ? east - 360.0 : east;
	if (lon_east) *lon_east = east;
	Mat b_filled(UTM_rows, UTM_cols, CV_8U); b_filled = 0;
	mapped_phase.create(UTM_rows, UTM_cols, CV_64F); mapped_phase = 0.0;
	//开始地理编码
	if (mapped_lat.type() == CV_32F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<float>(i, j);
				if (lon > 350.0) continue;
				lat = mapped_lat.at<float>(i, j);
				row = (int)round((lat - min_lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_phase.at<double>(row, col) = phase.at<double>(i, j);
				b_filled.at<uchar>(row, col) = 1;
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<double>(i, j);
				if (lon > 350.0) continue;
				lat = mapped_lat.at<double>(i, j);
				row = (int)round((lat - min_lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_phase.at<double>(row, col) = phase.at<double>(i, j);
				b_filled.at<uchar>(row, col) = 1;
			}
		}
	}
	

	//插值
	if (interpolation_method == 0)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < UTM_rows; i++)
		{
			for (int j = 0; j < UTM_cols; j++)
			{
				if (b_filled.at<uchar>(i, j) != 0) continue;
				int up, down, left, right, up_count, down_count, left_count, right_count;
				double value1, value2, ratio1, ratio2;
				//寻找上面有值的点
				up = i;
				while (true)
				{
					up--;
					if (up < 0) break;
					if (b_filled.at<uchar>(up, j) != 0) break;
				}
				//寻找下面有值的点
				down = i;
				while (true)
				{
					down++;
					if (down > UTM_rows - 1) break;
					if (b_filled.at<uchar>(down, j) != 0) break;
				}
				//寻找左边有值的点
				left = j;
				while (true)
				{
					left--;
					if (left < 0) break;
					if (b_filled.at<uchar>(i, left) != 0) break;
				}
				//寻找右边有值的点
				right = j;
				while (true)
				{
					right++;
					if (right > UTM_cols - 1) break;
					if (b_filled.at<uchar>(i, right) != 0) break;
				}

				//上下左右都有值
				if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
				{
					int x_i = i, x_j = j;
					int down_distance = down - i;
					int up_distance = i - up;
					int left_distance = j - left;
					int right_distance = right - j;
					if (down_distance < up_distance && down_distance < left_distance && down_distance < right_distance)
					{
						x_i = down; x_j = j;
					}
					else if (up_distance < down_distance && up_distance < left_distance && up_distance < right_distance)
					{
						x_i = up; x_j = j;
					}
					else if (left_distance < down_distance && left_distance < up_distance && left_distance < right_distance)
					{
						x_i = i; x_j = left;
					}
					else
					{
						x_i = i; x_j = right;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
				//上下有值
				if (up >= 0 && down <= UTM_rows - 1)
				{
					int x_i = i, x_j = j;
					int down_distance = down - i;
					int up_distance = i - up;
					if (down_distance < up_distance)
					{
						x_i = down; x_j = j;
					}
					else
					{
						x_i = down; x_j = j;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
				//左右有值
				if (left >= 0 && right <= UTM_cols - 1)
				{
					int x_i = i, x_j = j;
					int left_distance = j - left;
					int right_distance = right - j;
					if (left_distance < right_distance)
					{
						x_i = i; x_j = left;
					}
					else
					{
						x_i = i; x_j = right;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < UTM_rows; i++)
		{
			for (int j = 0; j < UTM_cols; j++)
			{
				if (b_filled.at<uchar>(i, j) != 0) continue;
				int up, down, left, right, up_count, down_count, left_count, right_count;
				double value1, value2, ratio1, ratio2;
				//寻找上面有值的点
				up = i;
				while (true)
				{
					up--;
					if (up < 0) break;
					if (b_filled.at<uchar>(up, j) != 0) break;
				}
				//寻找下面有值的点
				down = i;
				while (true)
				{
					down++;
					if (down > UTM_rows - 1) break;
					if (b_filled.at<uchar>(down, j) != 0) break;
				}
				//寻找左边有值的点
				left = j;
				while (true)
				{
					left--;
					if (left < 0) break;
					if (b_filled.at<uchar>(i, left) != 0) break;
				}
				//寻找右边有值的点
				right = j;
				while (true)
				{
					right++;
					if (right > UTM_cols - 1) break;
					if (b_filled.at<uchar>(i, right) != 0) break;
				}

				//上下左右都有值
				if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
				{
					
					ratio1 = double(j - left) / double(right - left);
					value1 = double(mapped_phase.at<double>(i, left)) +
						double(mapped_phase.at<double>(i, right) - mapped_phase.at<double>(i, right)) * ratio1;
					ratio2 = double(i - up) / double(down - up);
					value2 = double(mapped_phase.at<double>(up, j)) +
						double(mapped_phase.at<double>(down, j) - mapped_phase.at<double>(up, j)) * ratio2;
					mapped_phase.at<double>(i, j) = (value1 + value2) / 2.0;
					continue;
				}
				//上下有值
				if (up >= 0 && down <= UTM_rows - 1)
				{
					ratio2 = double(i - up) / double(down - up);
					value2 = double(mapped_phase.at<double>(up, j)) +
						double(mapped_phase.at<double>(down, j) - mapped_phase.at<double>(up, j)) * ratio2;
					mapped_phase.at<double>(i, j) = value2;
					continue;
				}
				//左右有值
				if (left >= 0 && right <= UTM_cols - 1)
				{
					ratio1 = double(j - left) / double(right - left);
					value1 = double(mapped_phase.at<double>(i, left)) +
						double(mapped_phase.at<double>(i, right) - mapped_phase.at<double>(i, right)) * ratio1;
					mapped_phase.at<double>(i, j) = value1;
					continue;
				}
			}
		}
	}
	

	cv::flip(mapped_phase, mapped_phase, 0);
	return 0;
}

int Utils::SAR2UTM(Mat& mapped_lon, Mat& mapped_lat, Mat& phase, Mat& mapped_phase, double grid_size, int interpolation_method, double* lon_east, double* lon_west, double* lat_north, double* lat_south)
{
	if (mapped_lat.size() != mapped_lon.size() ||
		mapped_lat.size() != phase.size() ||
		phase.rows < 2 ||
		phase.cols < 2 ||
		phase.type() != CV_64F ||
		mapped_lat.type() != mapped_lon.type() ||
		(mapped_lon.type() != CV_32F && mapped_lon.type() != CV_64F)
		)
	{
		fprintf(stderr, "SAR2UTM(): input check failed!\n");
		return -1;
	}
	//确定经纬度覆盖范围
	double max_lon = -380.0, min_lon = 380.0, max_lat = -380.0, min_lat = 180.0;
	if (mapped_lat.type() == CV_32F)
	{
		for (int i = 0; i < mapped_lat.rows; i++)
		{
			for (int j = 0; j < mapped_lat.cols; j++)
			{
				if (mapped_lat.at<float>(i, j) < 350.0)
				{
					min_lat = min_lat > mapped_lat.at<float>(i, j) ? mapped_lat.at<float>(i, j) : min_lat;
					max_lat = max_lat < mapped_lat.at<float>(i, j) ? mapped_lat.at<float>(i, j) : max_lat;
					min_lon = min_lon > mapped_lon.at<float>(i, j) ? mapped_lon.at<float>(i, j) : min_lon;
					max_lon = max_lon < mapped_lon.at<float>(i, j) ? mapped_lon.at<float>(i, j) : max_lon;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < mapped_lat.rows; i++)
		{
			for (int j = 0; j < mapped_lat.cols; j++)
			{
				if (mapped_lat.at<double>(i, j) < 350.0)
				{
					min_lat = min_lat > mapped_lat.at<double>(i, j) ? mapped_lat.at<double>(i, j) : min_lat;
					max_lat = max_lat < mapped_lat.at<double>(i, j) ? mapped_lat.at<double>(i, j) : max_lat;
					min_lon = min_lon > mapped_lon.at<double>(i, j) ? mapped_lon.at<double>(i, j) : min_lon;
					max_lon = max_lon < mapped_lon.at<double>(i, j) ? mapped_lon.at<double>(i, j) : max_lon;
				}
			}
		}
	}
	//cv::minMaxLoc(mapped_lon, &min_lon, &max_lon);
	//cv::minMaxLoc(mapped_lat, &min_lat, &max_lat);
	double west = max_lon - min_lon > 180.0 ? max_lon : min_lon;
	if (lon_west) *lon_west = west;

	//确定经纬度采样间隔
	double lon_interval, lat_interval;
	lon_interval = 5.0 / 6000.0 / (90.0 / grid_size);
	lat_interval = lon_interval;

	int rows = phase.rows; int cols = phase.cols;


	//计算UTM坐标系相位尺寸
	int UTM_rows = (max_lat - min_lat) / lat_interval;
	UTM_rows += 2;
	double south = max_lat - (double)(UTM_rows - 1) * lat_interval;
	double north = max_lat;
	if (lat_north) *lat_north = north;
	if (lat_south) *lat_south = south;
	double max_lon_temp = max_lon - min_lon;
	max_lon_temp = max_lon_temp > 180.0 ? 360.0 - max_lon_temp : max_lon_temp;
	int UTM_cols = max_lon_temp / lon_interval;
	UTM_cols += 2;
	double east = west + (double)(UTM_cols - 1) * lon_interval;
	east = east > 180.0 ? east - 360.0 : east;
	if (lon_east) *lon_east = east;
	Mat b_filled(UTM_rows, UTM_cols, CV_8U); b_filled = 0;
	mapped_phase.create(UTM_rows, UTM_cols, CV_64F); mapped_phase = 0.0;
	//开始地理编码
	if (mapped_lat.type() == CV_32F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<float>(i, j);
				if (lon > 350.0) continue;
				lat = mapped_lat.at<float>(i, j);
				row = (int)round((max_lat - lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_phase.at<double>(row, col) = phase.at<double>(i, j);
				b_filled.at<uchar>(row, col) = b_filled.at<uchar>(row, col) + 1;
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<double>(i, j);
				if (lon > 350.0) continue;
				lat = mapped_lat.at<double>(i, j);
				row = (int)round((max_lat - lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_phase.at<double>(row, col) = phase.at<double>(i, j);
				b_filled.at<uchar>(row, col) = b_filled.at<uchar>(row, col) + 1;
			}
		}
	}

	//插值
	if (interpolation_method == 0)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < UTM_rows; i++)
		{
			for (int j = 0; j < UTM_cols; j++)
			{
				if (b_filled.at<uchar>(i, j) != 0) continue;
				int up, down, left, right, up_count, down_count, left_count, right_count;
				double value1, value2, ratio1, ratio2;
				//寻找上面有值的点
				up = i;
				while (true)
				{
					up--;
					if (up < 0) break;
					if (b_filled.at<uchar>(up, j) != 0) break;
				}
				//寻找下面有值的点
				down = i;
				while (true)
				{
					down++;
					if (down > UTM_rows - 1) break;
					if (b_filled.at<uchar>(down, j) != 0) break;
				}
				//寻找左边有值的点
				left = j;
				while (true)
				{
					left--;
					if (left < 0) break;
					if (b_filled.at<uchar>(i, left) != 0) break;
				}
				//寻找右边有值的点
				right = j;
				while (true)
				{
					right++;
					if (right > UTM_cols - 1) break;
					if (b_filled.at<uchar>(i, right) != 0) break;
				}

				//上下左右都有值
				if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
				{
					int x_i = i, x_j = j;
					int down_distance = down - i;
					int up_distance = i - up;
					int left_distance = j - left;
					int right_distance = right - j;
					if (down_distance < up_distance && down_distance < left_distance && down_distance < right_distance)
					{
						x_i = down; x_j = j;
					}
					else if (up_distance < down_distance && up_distance < left_distance && up_distance < right_distance)
					{
						x_i = up; x_j = j;
					}
					else if (left_distance < down_distance && left_distance < up_distance && left_distance < right_distance)
					{
						x_i = i; x_j = left;
					}
					else
					{
						x_i = i; x_j = right;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
				//上下有值
				if (up >= 0 && down <= UTM_rows - 1)
				{
					int x_i = i, x_j = j;
					int down_distance = down - i;
					int up_distance = i - up;
					if (down_distance < up_distance)
					{
						x_i = down; x_j = j;
					}
					else
					{
						x_i = down; x_j = j;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
				//左右有值
				if (left >= 0 && right <= UTM_cols - 1)
				{
					int x_i = i, x_j = j;
					int left_distance = j - left;
					int right_distance = right - j;
					if (left_distance < right_distance)
					{
						x_i = i; x_j = left;
					}
					else
					{
						x_i = i; x_j = right;
					}
					mapped_phase.at<double>(i, j) = mapped_phase.at<double>(x_i, x_j);
					continue;
				}
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < UTM_rows; i++)
		{
			for (int j = 0; j < UTM_cols; j++)
			{
				if (b_filled.at<uchar>(i, j) != 0) continue;
				int up, down, left, right, up_count, down_count, left_count, right_count;
				double value1, value2, ratio1, ratio2;
				//寻找上面有值的点
				up = i;
				while (true)
				{
					up--;
					if (up < 0) break;
					if (b_filled.at<uchar>(up, j) != 0) break;
				}
				//寻找下面有值的点
				down = i;
				while (true)
				{
					down++;
					if (down > UTM_rows - 1) break;
					if (b_filled.at<uchar>(down, j) != 0) break;
				}
				//寻找左边有值的点
				left = j;
				while (true)
				{
					left--;
					if (left < 0) break;
					if (b_filled.at<uchar>(i, left) != 0) break;
				}
				//寻找右边有值的点
				right = j;
				while (true)
				{
					right++;
					if (right > UTM_cols - 1) break;
					if (b_filled.at<uchar>(i, right) != 0) break;
				}

				//上下左右都有值
				if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
				{

					ratio1 = double(j - left) / double(right - left);
					value1 = double(mapped_phase.at<double>(i, left)) +
						double(mapped_phase.at<double>(i, right) - mapped_phase.at<double>(i, right)) * ratio1;
					ratio2 = double(i - up) / double(down - up);
					value2 = double(mapped_phase.at<double>(up, j)) +
						double(mapped_phase.at<double>(down, j) - mapped_phase.at<double>(up, j)) * ratio2;
					mapped_phase.at<double>(i, j) = (value1 + value2) / 2.0;
					continue;
				}
				//上下有值
				if (up >= 0 && down <= UTM_rows - 1)
				{
					ratio2 = double(i - up) / double(down - up);
					value2 = double(mapped_phase.at<double>(up, j)) +
						double(mapped_phase.at<double>(down, j) - mapped_phase.at<double>(up, j)) * ratio2;
					mapped_phase.at<double>(i, j) = value2;
					continue;
				}
				//左右有值
				if (left >= 0 && right <= UTM_cols - 1)
				{
					ratio1 = double(j - left) / double(right - left);
					value1 = double(mapped_phase.at<double>(i, left)) +
						double(mapped_phase.at<double>(i, right) - mapped_phase.at<double>(i, right)) * ratio1;
					mapped_phase.at<double>(i, j) = value1;
					continue;
				}
			}
		}
	}


	//cv::flip(mapped_phase, mapped_phase, 0);
	return 0;
}

int Utils::SAR2UTM(
	Mat& mapped_lon,
	Mat& mapped_lat,
	ComplexMat& slc, 
	ComplexMat& mapped_slc, 
	int interpolation_method,
	double* lon_east,
	double* lon_west,
	double* lat_north,
	double* lat_south
)
{
	if (mapped_lat.size() != mapped_lon.size() ||
		mapped_lat.size() != slc.re.size() ||
		mapped_lat.size() != slc.im.size() ||
		slc.GetRows() < 2 ||
		slc.GetCols() < 2 ||
		(slc.type() != CV_16S && slc.type() != CV_32F) ||
		mapped_lat.type() != CV_32F ||
		mapped_lon.type() != CV_32F
		)
	{
		fprintf(stderr, "SAR2UTM(): input check failed!\n");
		return -1;
	}
	//确定经纬度覆盖范围
	double max_lon, min_lon, max_lat, min_lat;
	cv::minMaxLoc(mapped_lon, &min_lon, &max_lon);
	cv::minMaxLoc(mapped_lat, &min_lat, &max_lat);
	double west = max_lon - min_lon > 180.0 ? max_lon : min_lon;
	if (lon_west) *lon_west = west;
	//确定经纬度采样间隔
	double lon_interval, lat_interval;
	Mat temp1, temp2;
	mapped_lon(cv::Range(0, mapped_lon.rows - 1), cv::Range(0, mapped_lon.cols)).copyTo(temp1);
	mapped_lon(cv::Range(1, mapped_lon.rows), cv::Range(0, mapped_lon.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = temp1 / 180.0 * PI;
	wrap(temp1, temp1);
	temp1 = temp1 / PI * 180.0;
	temp1 = cv::abs(temp1);
	lon_interval = cv::mean(temp1)[0];
	double sigma1, sigma2;
	std(temp1, &sigma1);

	mapped_lon(cv::Range(0, mapped_lon.rows), cv::Range(0, mapped_lon.cols - 1)).copyTo(temp1);
	mapped_lon(cv::Range(0, mapped_lon.rows), cv::Range(1, mapped_lon.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = temp1 / 180.0 * PI;
	wrap(temp1, temp1);
	temp1 = temp1 / PI * 180.0;
	temp1 = cv::abs(temp1);
	std(temp1, &sigma2);
	lon_interval = (cv::mean(temp1)[0] + fabs(lon_interval)) / 2.0;
	lon_interval = lon_interval + 1.0 * (sigma1 + sigma2) / 2.0;

	mapped_lat(cv::Range(0, mapped_lat.rows - 1), cv::Range(0, mapped_lat.cols)).copyTo(temp1);
	mapped_lat(cv::Range(1, mapped_lat.rows), cv::Range(0, mapped_lat.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = cv::abs(temp1);
	lat_interval = cv::mean(temp1)[0];
	std(temp1, &sigma1);

	mapped_lat(cv::Range(0, mapped_lat.rows), cv::Range(0, mapped_lat.cols - 1)).copyTo(temp1);
	mapped_lat(cv::Range(0, mapped_lat.rows), cv::Range(1, mapped_lat.cols)).copyTo(temp2);
	temp1 = temp2 - temp1;
	temp1 = cv::abs(temp1);
	lat_interval = (cv::mean(temp1)[0] + lat_interval) / 2.0;
	std(temp1, &sigma2);
	lat_interval = lat_interval + 1.0 * (sigma1 + sigma2) / 2.0;

	int rows = slc.GetRows(); int cols = slc.GetCols();

	
	
	//计算UTM坐标系相位尺寸
	int UTM_rows = (max_lat - min_lat) / lat_interval;
	UTM_rows += 2;
	double south = max_lat - (double)(UTM_rows - 1) * lat_interval;
	double north = max_lat;
	if (lat_north) *lat_north = north;
	if (lat_south) *lat_south = south;
	double max_lon_temp = max_lon - min_lon;
	max_lon_temp = max_lon_temp > 180.0 ? 360.0 - max_lon_temp : max_lon_temp;
	int UTM_cols = max_lon_temp / lon_interval;
	UTM_cols += 2;
	double east = west + (double)(UTM_cols - 1) * lon_interval;
	east = east > 180.0 ? east - 360.0 : east;
	if (lon_east) *lon_east = east;
	Mat b_filled(UTM_rows, UTM_cols, CV_8U); b_filled = 0;
	mapped_slc.re.create(UTM_rows, UTM_cols, slc.type()); mapped_slc.im.create(UTM_rows, UTM_cols, slc.type());
	mapped_slc.re = 0;
	mapped_slc.im = 0;
	//开始地理编码
	if (slc.type() == CV_16S)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<float>(i, j);
				lat = mapped_lat.at<float>(i, j);
				row = (int)round((lat - min_lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_slc.re.at<short>(row, col) = slc.re.at<short>(i, j);
				mapped_slc.im.at<short>(row, col) = slc.im.at<short>(i, j);
				b_filled.at<uchar>(row, col) = 1;
			}
		}
	}
	else
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double lon, lat;
				int row, col;
				lon = mapped_lon.at<float>(i, j);
				lat = mapped_lat.at<float>(i, j);
				row = (int)round((lat - min_lat) / lat_interval);
				lon = fabs(lon - west);
				lon = lon > 180.0 ? 360.0 - lon : lon;
				col = (int)round(lon / lon_interval);
				mapped_slc.re.at<float>(row, col) = slc.re.at<float>(i, j);
				mapped_slc.im.at<float>(row, col) = slc.im.at<float>(i, j);
				b_filled.at<uchar>(row, col) = 1;
			}
		}
	}
	

	//插值
	if (interpolation_method == 0)
	{
		if (slc.type() == CV_16S)
		{
			for (int i = 0; i < UTM_rows; i++)
			{
				for (int j = 0; j < UTM_cols; j++)
				{
					if (b_filled.at<uchar>(i, j) != 0) continue;
					int up, down, left, right, up_count, down_count, left_count, right_count;
					double value1, value2, ratio1, ratio2;
					//寻找上面有值的点
					up = i;
					while (true)
					{
						up--;
						if (up < 0) break;
						if (b_filled.at<uchar>(up, j) != 0) break;
					}
					//寻找下面有值的点
					down = i;
					while (true)
					{
						down++;
						if (down > UTM_rows - 1) break;
						if (b_filled.at<uchar>(down, j) != 0) break;
					}
					//寻找左边有值的点
					left = j;
					while (true)
					{
						left--;
						if (left < 0) break;
						if (b_filled.at<uchar>(i, left) != 0) break;
					}
					//寻找右边有值的点
					right = j;
					while (true)
					{
						right++;
						if (right > UTM_cols - 1) break;
						if (b_filled.at<uchar>(i, right) != 0) break;
					}

					//上下左右都有值
					if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
					{
						int x_i = i, x_j = j;
						int down_distance = down - i;
						int up_distance = i - up;
						int left_distance = j - left;
						int right_distance = right - j;
						if (down_distance < up_distance && down_distance < left_distance && down_distance < right_distance)
						{
							x_i = down; x_j = j;
						}
						else if (up_distance < down_distance && up_distance < left_distance && up_distance < right_distance)
						{
							x_i = up; x_j = j;
						}
						else if (left_distance < down_distance && left_distance < up_distance && left_distance < right_distance)
						{
							x_i = i; x_j = left;
						}
						else
						{
							x_i = i; x_j = right;
						}
						mapped_slc.re.at<short>(i, j) = mapped_slc.re.at<short>(x_i, x_j);
						mapped_slc.im.at<short>(i, j) = mapped_slc.im.at<short>(x_i, x_j);
						continue;
					}
					//上下有值
					if (up >= 0 && down <= UTM_rows - 1)
					{
						int x_i = i, x_j = j;
						int down_distance = down - i;
						int up_distance = i - up;
						if (down_distance < up_distance)
						{
							x_i = down; x_j = j;
						}
						else
						{
							x_i = down; x_j = j;
						}
						mapped_slc.re.at<short>(i, j) = mapped_slc.re.at<short>(x_i, x_j);
						mapped_slc.im.at<short>(i, j) = mapped_slc.im.at<short>(x_i, x_j);
						continue;
					}
					//左右有值
					if (left >= 0 && right <= UTM_cols - 1)
					{
						int x_i = i, x_j = j;
						int left_distance = j - left;
						int right_distance = right - j;
						if (left_distance < right_distance)
						{
							x_i = i; x_j = left;
						}
						else
						{
							x_i = i; x_j = right;
						}
						mapped_slc.re.at<short>(i, j) = mapped_slc.re.at<short>(x_i, x_j);
						mapped_slc.im.at<short>(i, j) = mapped_slc.im.at<short>(x_i, x_j);
						continue;
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < UTM_rows; i++)
			{
				for (int j = 0; j < UTM_cols; j++)
				{
					if (b_filled.at<uchar>(i, j) != 0) continue;
					int up, down, left, right, up_count, down_count, left_count, right_count;
					double value1, value2, ratio1, ratio2;
					//寻找上面有值的点
					up = i;
					while (true)
					{
						up--;
						if (up < 0) break;
						if (b_filled.at<uchar>(up, j) != 0) break;
					}
					//寻找下面有值的点
					down = i;
					while (true)
					{
						down++;
						if (down > UTM_rows - 1) break;
						if (b_filled.at<uchar>(down, j) != 0) break;
					}
					//寻找左边有值的点
					left = j;
					while (true)
					{
						left--;
						if (left < 0) break;
						if (b_filled.at<uchar>(i, left) != 0) break;
					}
					//寻找右边有值的点
					right = j;
					while (true)
					{
						right++;
						if (right > UTM_cols - 1) break;
						if (b_filled.at<uchar>(i, right) != 0) break;
					}

					//上下左右都有值
					if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
					{
						int x_i = i, x_j = j;
						int down_distance = down - i;
						int up_distance = i - up;
						int left_distance = j - left;
						int right_distance = right - j;
						if (down_distance < up_distance && down_distance < left_distance && down_distance < right_distance)
						{
							x_i = down; x_j = j;
						}
						else if (up_distance < down_distance && up_distance < left_distance && up_distance < right_distance)
						{
							x_i = up; x_j = j;
						}
						else if (left_distance < down_distance && left_distance < up_distance && left_distance < right_distance)
						{
							x_i = i; x_j = left;
						}
						else
						{
							x_i = i; x_j = right;
						}
						mapped_slc.re.at<float>(i, j) = mapped_slc.re.at<float>(x_i, x_j);
						mapped_slc.im.at<float>(i, j) = mapped_slc.im.at<float>(x_i, x_j);
						continue;
					}
					//上下有值
					if (up >= 0 && down <= UTM_rows - 1)
					{
						int x_i = i, x_j = j;
						int down_distance = down - i;
						int up_distance = i - up;
						if (down_distance < up_distance)
						{
							x_i = down; x_j = j;
						}
						else
						{
							x_i = down; x_j = j;
						}
						mapped_slc.re.at<float>(i, j) = mapped_slc.re.at<float>(x_i, x_j);
						mapped_slc.im.at<float>(i, j) = mapped_slc.im.at<float>(x_i, x_j);
						continue;
					}
					//左右有值
					if (left >= 0 && right <= UTM_cols - 1)
					{
						int x_i = i, x_j = j;
						int left_distance = j - left;
						int right_distance = right - j;
						if (left_distance < right_distance)
						{
							x_i = i; x_j = left;
						}
						else
						{
							x_i = i; x_j = right;
						}
						mapped_slc.re.at<float>(i, j) = mapped_slc.re.at<float>(x_i, x_j);
						mapped_slc.im.at<float>(i, j) = mapped_slc.im.at<float>(x_i, x_j);
						continue;
					}
				}
			}
		}
		
	}
	else
	{
		if (slc.type() == CV_16S)
		{
			for (int i = 0; i < UTM_rows; i++)
			{
				for (int j = 0; j < UTM_cols; j++)
				{
					if (b_filled.at<uchar>(i, j) != 0) continue;
					int up, down, left, right, up_count, down_count, left_count, right_count;
					double value1, value2, ratio1, ratio2;
					//寻找上面有值的点
					up = i;
					while (true)
					{
						up--;
						if (up < 0) break;
						if (b_filled.at<uchar>(up, j) != 0) break;
					}
					//寻找下面有值的点
					down = i;
					while (true)
					{
						down++;
						if (down > UTM_rows - 1) break;
						if (b_filled.at<uchar>(down, j) != 0) break;
					}
					//寻找左边有值的点
					left = j;
					while (true)
					{
						left--;
						if (left < 0) break;
						if (b_filled.at<uchar>(i, left) != 0) break;
					}
					//寻找右边有值的点
					right = j;
					while (true)
					{
						right++;
						if (right > UTM_cols - 1) break;
						if (b_filled.at<uchar>(i, right) != 0) break;
					}

					//上下左右都有值
					if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
					{

						ratio1 = double(j - left) / double(right - left);
						value1 = double(mapped_slc.re.at<short>(i, left)) +
							double(mapped_slc.re.at<short>(i, right) - mapped_slc.re.at<short>(i, right)) * ratio1;
						ratio2 = double(i - up) / double(down - up);
						value2 = double(mapped_slc.re.at<short>(up, j)) +
							double(mapped_slc.re.at<short>(down, j) - mapped_slc.re.at<short>(up, j)) * ratio2;
						mapped_slc.re.at<short>(i, j) = (value1 + value2) / 2.0;

						value1 = double(mapped_slc.im.at<short>(i, left)) +
							double(mapped_slc.im.at<short>(i, right) - mapped_slc.im.at<short>(i, right)) * ratio1;
						value2 = double(mapped_slc.im.at<short>(up, j)) +
							double(mapped_slc.im.at<short>(down, j) - mapped_slc.im.at<short>(up, j)) * ratio2;
						mapped_slc.im.at<short>(i, j) = (value1 + value2) / 2.0;
						continue;
					}
					//上下有值
					if (up >= 0 && down <= UTM_rows - 1)
					{
						ratio2 = double(i - up) / double(down - up);
						value2 = double(mapped_slc.re.at<short>(up, j)) +
							double(mapped_slc.re.at<short>(down, j) - mapped_slc.re.at<short>(up, j)) * ratio2;
						mapped_slc.re.at<short>(i, j) = value2;

						value2 = double(mapped_slc.im.at<short>(up, j)) +
							double(mapped_slc.im.at<short>(down, j) - mapped_slc.im.at<short>(up, j)) * ratio2;
						mapped_slc.im.at<short>(i, j) = value2;
						continue;
					}
					//左右有值
					if (left >= 0 && right <= UTM_cols - 1)
					{
						ratio1 = double(j - left) / double(right - left);
						value1 = double(mapped_slc.re.at<short>(i, left)) +
							double(mapped_slc.re.at<short>(i, right) - mapped_slc.re.at<short>(i, right)) * ratio1;
						mapped_slc.re.at<short>(i, j) = value1;

						value1 = double(mapped_slc.im.at<short>(i, left)) +
							double(mapped_slc.im.at<short>(i, right) - mapped_slc.im.at<short>(i, right)) * ratio1;
						mapped_slc.im.at<short>(i, j) = value1;
						continue;
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < UTM_rows; i++)
			{
				for (int j = 0; j < UTM_cols; j++)
				{
					if (b_filled.at<uchar>(i, j) != 0) continue;
					int up, down, left, right, up_count, down_count, left_count, right_count;
					double value1, value2, ratio1, ratio2;
					//寻找上面有值的点
					up = i;
					while (true)
					{
						up--;
						if (up < 0) break;
						if (b_filled.at<uchar>(up, j) != 0) break;
					}
					//寻找下面有值的点
					down = i;
					while (true)
					{
						down++;
						if (down > UTM_rows - 1) break;
						if (b_filled.at<uchar>(down, j) != 0) break;
					}
					//寻找左边有值的点
					left = j;
					while (true)
					{
						left--;
						if (left < 0) break;
						if (b_filled.at<uchar>(i, left) != 0) break;
					}
					//寻找右边有值的点
					right = j;
					while (true)
					{
						right++;
						if (right > UTM_cols - 1) break;
						if (b_filled.at<uchar>(i, right) != 0) break;
					}

					//上下左右都有值
					if (left >= 0 && right <= UTM_cols - 1 && up >= 0 && down <= UTM_rows - 1)
					{

						ratio1 = double(j - left) / double(right - left);
						value1 = double(mapped_slc.re.at<float>(i, left)) +
							double(mapped_slc.re.at<float>(i, right) - mapped_slc.re.at<float>(i, right)) * ratio1;
						ratio2 = double(i - up) / double(down - up);
						value2 = double(mapped_slc.re.at<float>(up, j)) +
							double(mapped_slc.re.at<float>(down, j) - mapped_slc.re.at<float>(up, j)) * ratio2;
						mapped_slc.re.at<float>(i, j) = (value1 + value2) / 2.0;

						value1 = double(mapped_slc.im.at<float>(i, left)) +
							double(mapped_slc.im.at<float>(i, right) - mapped_slc.im.at<float>(i, right)) * ratio1;
						value2 = double(mapped_slc.im.at<float>(up, j)) +
							double(mapped_slc.im.at<float>(down, j) - mapped_slc.im.at<float>(up, j)) * ratio2;
						mapped_slc.im.at<float>(i, j) = (value1 + value2) / 2.0;
						continue;
					}
					//上下有值
					if (up >= 0 && down <= UTM_rows - 1)
					{
						ratio2 = double(i - up) / double(down - up);
						value2 = double(mapped_slc.re.at<float>(up, j)) +
							double(mapped_slc.re.at<float>(down, j) - mapped_slc.re.at<float>(up, j)) * ratio2;
						mapped_slc.re.at<float>(i, j) = value2;

						value2 = double(mapped_slc.im.at<float>(up, j)) +
							double(mapped_slc.im.at<float>(down, j) - mapped_slc.im.at<float>(up, j)) * ratio2;
						mapped_slc.im.at<float>(i, j) = value2;
						continue;
					}
					//左右有值
					if (left >= 0 && right <= UTM_cols - 1)
					{
						ratio1 = double(j - left) / double(right - left);
						value1 = double(mapped_slc.re.at<float>(i, left)) +
							double(mapped_slc.re.at<float>(i, right) - mapped_slc.re.at<float>(i, right)) * ratio1;
						mapped_slc.re.at<float>(i, j) = value1;

						value1 = double(mapped_slc.im.at<float>(i, left)) +
							double(mapped_slc.im.at<float>(i, right) - mapped_slc.im.at<float>(i, right)) * ratio1;
						mapped_slc.im.at<float>(i, j) = value1;
						continue;
					}
				}
			}
		}
	}


	cv::flip(mapped_slc.re, mapped_slc.re, 0);
	cv::flip(mapped_slc.im, mapped_slc.im, 0);
	return 0;
}

int Utils::geocode(
	Mat& DEM84,
	Mat& input,
	double mapped_resolution_x,
	double mapped_resolution_y,
	Mat& mapped_result,
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
	double lon_spacing,
	double lat_spacing,
	double* lon_east,
	double* lon_west,
	double* lat_north, 
	double* lat_south
)
{
	if (DEM84.empty() ||
		DEM84.type() != CV_16S ||
		(input.type() != CV_32F && input.type() != CV_64F)||
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
		fprintf(stderr, "geocode(): input check failed!\n");
		return -1;
	}
	//计算DEM插值倍数
	int interp_times_x, interp_times_y;
	double a = 6378137, b = 6356752;
	double C_short = (a + b) * PI;//经线一圈长度
	double C_long = a * 2 * PI;
	double lon_per_meter = 360.0 / C_short;//经线上每米多少度
	double lat_per_meter = 360.0 / (C_long * cos(lat_upperleft / 180.0 * PI));//纬线上每米多少度
	interp_times_x = lon_spacing / lon_per_meter / mapped_resolution_x;
	interp_times_y = lat_spacing / lat_per_meter / mapped_resolution_y;
	//84坐标系DEM插值
	Mat DEM, stateVector_interp;
	interp_times_x = interp_times_x < 1 ? 1 : interp_times_x;
	interp_times_y = interp_times_y < 1 ? 1 : interp_times_y;
	cv::resize(DEM84, DEM, cv::Size(DEM84.cols * interp_times_x, DEM84.rows * interp_times_y));
	lon_spacing = lon_spacing / (double)interp_times_x;
	lat_spacing = lat_spacing / (double)interp_times_y;
	//初始化轨道类
	orbitStateVectors stateVectors(stateVector, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;

	int DEM_rows = DEM.rows; int DEM_cols = DEM.cols;
	double dopplerFrequency = 0.0;
	mapped_result.create(DEM.rows, DEM.cols, input.type()); mapped_result = 0.0;
	//地理编码
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
			int azimuthIndex = floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
			int rangeIndex = floor((distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing);
			azimuthIndex = azimuthIndex - offset_row;
			rangeIndex = rangeIndex - offset_col;
			if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 2 || rangeIndex < 0 || rangeIndex > sceneWidth - 2)
			{

			}
			else
			{
				//双线性插值计算
				double ratio_x = (distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing - floor((distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing);
				double ratio_y = (zeroDopplerTime - acquisitionStartTime) / time_interval - floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
				
				if (input.type() == CV_32F)
				{
					double upper = (double)input.at<float>(azimuthIndex, rangeIndex) + double(input.at<float>(azimuthIndex, rangeIndex + 1)
						- input.at<float>(azimuthIndex, rangeIndex)) * ratio_x;
					double lower = (double)input.at<float>(azimuthIndex + 1, rangeIndex) + double(input.at<float>(azimuthIndex + 1, rangeIndex + 1)
						- input.at<float>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_result.at<float>(i, j) = upper + (lower - upper) * ratio_y;

					/*mapped_result.at<float>(i, j) = input.at<float>(azimuthIndex + 1, rangeIndex + 1) + input.at<float>(azimuthIndex + 1, rangeIndex) +
						input.at<float>(azimuthIndex, rangeIndex + 1) + input.at<float>(azimuthIndex, rangeIndex);
					mapped_result.at<float>(i, j) = mapped_result.at<float>(i, j) / 4.0;*/
				}
				else
				{
					double upper = (double)input.at<double>(azimuthIndex, rangeIndex) + double(input.at<double>(azimuthIndex, rangeIndex + 1)
						- input.at<double>(azimuthIndex, rangeIndex)) * ratio_x;
					double lower = (double)input.at<double>(azimuthIndex + 1, rangeIndex) + double(input.at<double>(azimuthIndex + 1, rangeIndex + 1)
						- input.at<double>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_result.at<double>(i, j) = upper + (lower - upper) * ratio_y;

					/*mapped_result.at<double>(i, j) = input.at<double>(azimuthIndex + 1, rangeIndex + 1) + input.at<double>(azimuthIndex + 1, rangeIndex) +
						input.at<double>(azimuthIndex, rangeIndex + 1) + input.at<double>(azimuthIndex, rangeIndex);
					mapped_result.at<double>(i, j) = mapped_result.at<double>(i, j) / 4.0;*/
				}
			}
		}
	}
	if (lat_north && lat_south && lon_east && lon_west)
	{
		*lat_north = lat_upperleft;
		*lat_south = lat_upperleft - (double)(DEM.rows - 1) * lat_spacing;
		*lon_east = lon_upperleft;
		*lon_west = lon_upperleft + (double)(DEM.cols - 1) * lon_spacing;
		*lon_west = *lon_west > 180.0 ? (*lon_west - 360.0) : *lon_west;
	}
	return 0;
}

int Utils::geocode(
	Mat& DEM84,
	ComplexMat& slc,
	double mapped_resolution_x,
	double mapped_resolution_y,
	ComplexMat& mapped_slc, 
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
	double lon_spacing,
	double lat_spacing,
	double* lon_east,
	double* lon_west,
	double* lat_north,
	double* lat_south
)
{
	if (DEM84.empty() ||
		DEM84.type() != CV_16S ||
		(slc.type() != CV_32F && slc.type() != CV_16S) ||
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
		fprintf(stderr, "geocode(): input check failed!\n");
		return -1;
	}
	//计算DEM插值倍数
	int interp_times_x, interp_times_y;
	double a = 6378137, b = 6356752;
	double C_short = (a + b) * PI;//经线一圈长度
	double C_long = a * 2 * PI;
	double lon_per_meter = 360.0 / C_short;//经线上每米多少度
	double lat_per_meter = 360.0 / (C_long * cos(lat_upperleft / 180.0 * PI));//纬线上每米多少度
	interp_times_x = lon_spacing / lon_per_meter / mapped_resolution_x;
	interp_times_y = lat_spacing / lat_per_meter / mapped_resolution_y;
	//84坐标系DEM插值
	Mat DEM, stateVector_interp;
	interp_times_x = interp_times_x < 1 ? 1 : interp_times_x;
	interp_times_y = interp_times_y < 1 ? 1 : interp_times_y;
	cv::resize(DEM84, DEM, cv::Size(DEM84.cols * interp_times_x, DEM84.rows * interp_times_y));
	lon_spacing = lon_spacing / (double)interp_times_x;
	lat_spacing = lat_spacing / (double)interp_times_y;
	//初始化轨道类
	orbitStateVectors stateVectors(stateVector, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;

	int DEM_rows = DEM.rows; int DEM_cols = DEM.cols;
	double dopplerFrequency = 0.0;
	mapped_slc.re.create(DEM.rows, DEM.cols, slc.type()); mapped_slc.re = 0.0;
	mapped_slc.im.create(DEM.rows, DEM.cols, slc.type()); mapped_slc.im = 0.0;
	//地理编码
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
			int azimuthIndex = floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
			int rangeIndex = floor((distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing);
			azimuthIndex = azimuthIndex - offset_row;
			rangeIndex = rangeIndex - offset_col;
			if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 2 || rangeIndex < 0 || rangeIndex > sceneWidth - 2)
			{

			}
			else
			{
				//双线性插值计算
				double ratio_x = (distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing - floor((distance - nearRangeTime * VEL_C * 0.5) / rangeSpacing);
				double ratio_y = (zeroDopplerTime - acquisitionStartTime) / time_interval - floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
				if (slc.type() == CV_32F)
				{
					double upper = (double)slc.re.at<float>(azimuthIndex, rangeIndex) + double(slc.re.at<float>(azimuthIndex, rangeIndex + 1)
						- slc.re.at<float>(azimuthIndex, rangeIndex)) * ratio_x;
					double lower = (double)slc.re.at<float>(azimuthIndex + 1, rangeIndex) + double(slc.re.at<float>(azimuthIndex + 1, rangeIndex + 1)
						- slc.re.at<float>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_slc.re.at<float>(i, j) = upper + (lower - upper) * ratio_y;

					upper = (double)slc.im.at<float>(azimuthIndex, rangeIndex) + double(slc.im.at<float>(azimuthIndex, rangeIndex + 1)
						- slc.im.at<float>(azimuthIndex, rangeIndex)) * ratio_x;
					lower = (double)slc.im.at<float>(azimuthIndex + 1, rangeIndex) + double(slc.im.at<float>(azimuthIndex + 1, rangeIndex + 1)
						- slc.im.at<float>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_slc.im.at<float>(i, j) = upper + (lower - upper) * ratio_y;

					/*mapped_slc.im.at<float>(i, j) = slc.im.at<float>(azimuthIndex + 1, rangeIndex + 1) + slc.im.at<float>(azimuthIndex + 1, rangeIndex) +
						slc.im.at<float>(azimuthIndex, rangeIndex + 1) + slc.im.at<float>(azimuthIndex, rangeIndex);
					mapped_slc.im.at<float>(i, j) = mapped_slc.im.at<float>(i, j) / 4.0;

					mapped_slc.re.at<float>(i, j) = slc.re.at<float>(azimuthIndex + 1, rangeIndex + 1) + slc.re.at<float>(azimuthIndex + 1, rangeIndex) +
						slc.re.at<float>(azimuthIndex, rangeIndex + 1) + slc.re.at<float>(azimuthIndex, rangeIndex);
					mapped_slc.re.at<float>(i, j) = mapped_slc.re.at<float>(i, j) / 4.0;*/
				}
				else
				{
					double upper = (double)slc.re.at<short>(azimuthIndex, rangeIndex) + double(slc.re.at<short>(azimuthIndex, rangeIndex + 1)
						- slc.re.at<short>(azimuthIndex, rangeIndex)) * ratio_x;
					double lower = (double)slc.re.at<short>(azimuthIndex + 1, rangeIndex) + double(slc.re.at<short>(azimuthIndex + 1, rangeIndex + 1)
						- slc.re.at<short>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_slc.re.at<short>(i, j) = upper + (lower - upper) * ratio_y;

					upper = (double)slc.im.at<short>(azimuthIndex, rangeIndex) + double(slc.im.at<short>(azimuthIndex, rangeIndex + 1)
						- slc.im.at<short>(azimuthIndex, rangeIndex)) * ratio_x;
					lower = (double)slc.im.at<short>(azimuthIndex + 1, rangeIndex) + double(slc.im.at<short>(azimuthIndex + 1, rangeIndex + 1)
						- slc.im.at<short>(azimuthIndex + 1, rangeIndex)) * ratio_x;
					mapped_slc.im.at<short>(i, j) = upper + (lower - upper) * ratio_y;

					/*mapped_slc.im.at<short>(i, j) = slc.im.at<short>(azimuthIndex + 1, rangeIndex + 1) + slc.im.at<short>(azimuthIndex + 1, rangeIndex) +
						slc.im.at<short>(azimuthIndex, rangeIndex + 1) + slc.im.at<short>(azimuthIndex, rangeIndex);
					mapped_slc.im.at<short>(i, j) = double(mapped_slc.im.at<short>(i, j)) / 4.0;

					mapped_slc.re.at<short>(i, j) = slc.re.at<short>(azimuthIndex + 1, rangeIndex + 1) + slc.re.at<short>(azimuthIndex + 1, rangeIndex) +
						slc.re.at<short>(azimuthIndex, rangeIndex + 1) + slc.re.at<short>(azimuthIndex, rangeIndex);
					mapped_slc.re.at<short>(i, j) = double(mapped_slc.re.at<short>(i, j)) / 4.0;*/
				}
			}
		}
	}
	if (lat_north && lat_south && lon_east && lon_west)
	{
		*lat_north = lat_upperleft;
		*lat_south = lat_upperleft - (double)(DEM.rows - 1) * lat_spacing;
		*lon_east = lon_upperleft;
		*lon_west = lon_upperleft + (double)(DEM.cols - 1) * lon_spacing;
		*lon_west = *lon_west > 180.0 ? (*lon_west - 360.0) : *lon_west;
	}
	return 0;
}

int Utils::geo_transformation(
	const char* grille_file,
	Mat DTM, double xllcorner,
	double yllcorner,
	Mat& prior_DTM,
	Mat& mapped_DTM,
	int SAR_extent_x,
	int SAR_extent_y
)
{
	if (!grille_file || DTM.empty() || xllcorner < 0 || yllcorner < 0)
	{
		fprintf(stderr, "geo_transformation(): input check failed!\n");
		return -1;
	}
	int ret;
	if (DTM.type() != CV_64F)
	{
		DTM.convertTo(DTM, CV_64F);
	}
	FormatConversion conversion;
	Mat row_matrix, col_matrix, utm_x, utm_y;
	vector<Mat> lon_matrix;
	vector<Mat> lat_matrix;
	vector<double> height_vector;
	mapped_DTM.create(SAR_extent_y, SAR_extent_x, CV_64F);
	mapped_DTM = -1.0;
	ret = read_grille(grille_file, row_matrix, col_matrix, lon_matrix, lat_matrix, height_vector);
	if (return_check(ret, "read_grille()", error_head)) return -1;
	//cvmat2bin("D:\\working_dir\\projects\\software\\InSAR\\bin\\row_matrix.bin", row_matrix);
	//cvmat2bin("D:\\working_dir\\projects\\software\\InSAR\\bin\\col_matrix.bin", col_matrix);
	//经纬度转UTM
	for (int i = 0; i < height_vector.size(); i++)
	{
		ret = lonlat2utm(lon_matrix[i], lat_matrix[i], utm_x, utm_y);
		if (return_check(ret, "lonlat2utm()", error_head)) return -1;
		//conversion.creat_new_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5");
		//conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5", "X", utm_x);
		//conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5", "Y", utm_y);
		utm_x.copyTo(lon_matrix[i]);
		utm_y.copyTo(lat_matrix[i]);
	}
	int DTM_rows = DTM.rows;
	int DTM_cols = DTM.cols;
	//DTM逐点转换
	Mat DTM_mapped_X, DTM_mapped_Y;
	DTM_mapped_X = Mat::zeros(DTM_rows, DTM_cols, CV_64F);
	DTM_mapped_Y = Mat::zeros(DTM_rows, DTM_cols, CV_64F);
	DTM_mapped_X = -1;
	DTM_mapped_Y = -1;
	volatile int count = 0;
	Mat mask = Mat::zeros(DTM_rows, DTM_cols, CV_8UC1);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < SAR_extent_y; i++)
	{
		for (int j = 0; j < SAR_extent_x; j++)
		{
			//首先确定DTM值是否有效
			double h = prior_DTM.at<double>(i, j);
			//double h = 20.6;
			if (h < -9000)
			{
				continue;
			}
			mask.at<uchar>(i, j) = 1;
			//定位相邻的网格层
			int low_ix, high_ix;
			for (int k = 0; k < height_vector.size() - 1; k++)
			{
				if (h >= height_vector[k] && h <= height_vector[k + 1])
				{
					low_ix = k;
					high_ix = k + 1;
					//对下层网格寻找定位点
					bool located1 = false;
					for (int ii = 0; ii < row_matrix.rows - 1; ii++)
					{
						for (int jj = 0; jj < row_matrix.cols - 1; jj++)
						{
							//located1 = false; located2 = false;
							//判断该点是否在四个点中间
							double Mx, My, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy;
							Ax = row_matrix.at<double>(ii, jj); Ay = col_matrix.at<double>(ii, jj);
							Bx = row_matrix.at<double>(ii + 1, jj); Dy = col_matrix.at<double>(ii, jj + 1);

							if (i >= Ax && 
								i <= Bx &&
								j >= Ay &&
								j <= Dy &&
								!located1
								)
							{
								located1 = true;
								//线性插值得到在下层网格上的UTM坐标
								//UTM_x插值
								double UTM_x_upleft = lon_matrix[low_ix].at<double>(ii, jj);
								double UTM_x_upright = lon_matrix[low_ix].at<double>(ii, jj + 1);
								double UTM_x_lowleft = lon_matrix[low_ix].at<double>(ii + 1, jj);
								double UTM_x_lowright = lon_matrix[low_ix].at<double>(ii + 1, jj + 1);
								double upper = UTM_x_upleft + (UTM_x_upright - UTM_x_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) * 
									(j - col_matrix.at<double>(ii, jj));
								double lower = UTM_x_lowleft + (UTM_x_lowright - UTM_x_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
									(j - col_matrix.at<double>(ii + 1, jj));
								double UTM_x_final_lower = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
									(i - row_matrix.at<double>(ii + 1, jj));

								//UTM_y插值
								double UTM_y_upleft = lat_matrix[low_ix].at<double>(ii, jj);
								double UTM_y_upright = lat_matrix[low_ix].at<double>(ii, jj + 1);
								double UTM_y_lowleft = lat_matrix[low_ix].at<double>(ii + 1, jj);
								double UTM_y_lowright = lat_matrix[low_ix].at<double>(ii + 1, jj + 1);
								upper = UTM_y_upleft + (UTM_y_upright - UTM_y_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
									(j - col_matrix.at<double>(ii, jj));
								lower = UTM_y_lowleft + (UTM_y_lowright - UTM_y_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
									(j - col_matrix.at<double>(ii + 1, jj));
								double UTM_y_final_lower = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
									(i - row_matrix.at<double>(ii + 1, jj));

								//线性插值得到在上层网格上的UTM坐标
								//UTM_x插值
								UTM_x_upleft = lon_matrix[high_ix].at<double>(ii, jj);
								UTM_x_upright = lon_matrix[high_ix].at<double>(ii, jj + 1);
								UTM_x_lowleft = lon_matrix[high_ix].at<double>(ii + 1, jj);
								UTM_x_lowright = lon_matrix[high_ix].at<double>(ii + 1, jj + 1);
								upper = UTM_x_upleft + (UTM_x_upright - UTM_x_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
									(j - col_matrix.at<double>(ii, jj));
								lower = UTM_x_lowleft + (UTM_x_lowright - UTM_x_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
									(j - col_matrix.at<double>(ii + 1, jj));
								double UTM_x_final_higher = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
									(i - row_matrix.at<double>(ii + 1, jj));

								//UTM_y插值
								UTM_y_upleft = lat_matrix[high_ix].at<double>(ii, jj);
								UTM_y_upright = lat_matrix[high_ix].at<double>(ii, jj + 1);
								UTM_y_lowleft = lat_matrix[high_ix].at<double>(ii + 1, jj);
								UTM_y_lowright = lat_matrix[high_ix].at<double>(ii + 1, jj + 1);
								upper = UTM_y_upleft + (UTM_y_upright - UTM_y_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
									(j - col_matrix.at<double>(ii, jj));
								lower = UTM_y_lowleft + (UTM_y_lowright - UTM_y_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
									(j - col_matrix.at<double>(ii + 1, jj));
								double UTM_y_final_higher = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
									(i - row_matrix.at<double>(ii + 1, jj));



								//上下两层之间插值得到UTM_x和UTM_y
								double UTM_x_final = UTM_x_final_lower + (UTM_x_final_higher - UTM_x_final_lower) / (height_vector[high_ix] - height_vector[low_ix]) *
									(h - height_vector[low_ix]);
								double UTM_y_final = UTM_y_final_lower + (UTM_y_final_higher - UTM_y_final_lower) / (height_vector[high_ix] - height_vector[low_ix]) *
									(h - height_vector[low_ix]);



								//通过插值得到的UTM_x和UTM_y再次插值得到DTM
								if (UTM_x_final >= xllcorner &&
									UTM_x_final < (xllcorner + DTM_cols - 1) &&
									UTM_y_final >= yllcorner &&
									UTM_y_final < (yllcorner + DTM_rows - 1)
									)
								{
									int row = floor((yllcorner + DTM_rows - 1) - UTM_y_final);
									int col = floor(UTM_x_final - xllcorner);
									if (DTM.at<double>(row, col) > -10 &&
										DTM.at<double>(row + 1, col) > -10 &&
										DTM.at<double>(row + 1, col + 1) > -10 &&
										DTM.at<double>(row, col + 1) > -10
										)
									{
										upper = DTM.at<double>(row, col) + (DTM.at<double>(row, col + 1) - DTM.at<double>(row, col)) / 1 * ((yllcorner + DTM_rows - 1) - UTM_y_final - row);
										lower = DTM.at<double>(row + 1, col) + (DTM.at<double>(row + 1, col + 1) - DTM.at<double>(row + 1, col)) / 1 *
											((yllcorner + DTM_rows - 1) - UTM_y_final - row);
										mapped_DTM.at<double>(i, j) = lower + (upper - lower) / 1 * (UTM_x_final - xllcorner - col);
									}
									
								}
								break;
							}
						}
						if (located1)break;
					}
				}
			}
		}
		count++;
		if (count % 10 == 0)
		{
			printf("\r估计进度1：%lf%%", double(count) / double(SAR_extent_y) * 100.0);
			fflush(stdout);
		}
	}
	
//	conversion.creat_new_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\DTM_mapped.h5");
//	conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\DTM_mapped.h5", "X", DTM_mapped_X);
//	conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\DTM_mapped.h5", "Y", DTM_mapped_Y);
//
//	conversion.read_array_from_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\DTM_mapped.h5", "X", DTM_mapped_X);
//	conversion.read_array_from_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\DTM_mapped.h5", "Y", DTM_mapped_Y);
//
//	
//	//最邻近插值法得到SAR坐标系DTM
////#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < DTM_rows; i++)
//	{
//		for (int j = 0; j < DTM_cols; j++)
//		{
//			if (DTM_mapped_X.at<double>(i, j) > 0 && DTM_mapped_Y.at<double>(i, j) > 0 && DTM.at<double>(i, j) > -1)
//			{
//				int col = round(DTM_mapped_Y.at<double>(i, j));
//				int row = round(DTM_mapped_X.at<double>(i, j));
//				mapped_DTM.at<double>(row, col) = DTM.at<double>(i, j);
//			}
//		}
//	}
	
//	//插值得到SAR坐标系的DTM
//	
//	double minX1, maxX1, minY1, maxY1;
//	cv::minMaxIdx(DTM_mapped_X, &minX1, &maxX1, NULL, NULL, mask);
//	cv::minMaxIdx(DTM_mapped_Y, &minY1, &maxY1, NULL, NULL, mask);
//	int start_y = floor(minY1) < 0 ? 0 : floor(minY1);
//	int end_y = ceil(maxY1);
//	int start_x = floor(minX1) < 0 ? 0 : floor(minX1);
//	int end_x = ceil(maxX1);
//	count = 0;
//#pragma omp parallel for schedule(guided)
//	for (int i = start_y; i < end_y; i++)
//	{
//		for (int j = start_x; j < end_x; j++)
//		{
//			bool located = false;
//			for (int ii = 0; ii < DTM_rows - 1; ii++)
//			{
//				for (int jj = 0; jj < DTM_cols - 1; jj++)
//				{
//					Mat X_tmp = Mat::zeros(1, 4, CV_64F), Y_tmp = Mat::zeros(1, 4, CV_64F);
//					X_tmp.at<double>(0, 0) = DTM_mapped_X.at<double>(ii, jj); Y_tmp.at<double>(0, 0) = DTM_mapped_Y.at<double>(ii, jj);
//					X_tmp.at<double>(0, 1) = DTM_mapped_X.at<double>(ii, jj + 1); Y_tmp.at<double>(0, 1) = DTM_mapped_Y.at<double>(ii, jj + 1);
//					X_tmp.at<double>(0, 2) = DTM_mapped_X.at<double>(ii + 1, jj + 1); Y_tmp.at<double>(0, 2) = DTM_mapped_Y.at<double>(ii + 1, jj + 1);
//					X_tmp.at<double>(0, 3) = DTM_mapped_X.at<double>(ii + 1, jj); Y_tmp.at<double>(0, 3) = DTM_mapped_Y.at<double>(ii + 1, jj);
//					double minX, minY, maxX, maxY;
//					cv::minMaxIdx(X_tmp, &minX, &maxX);
//					cv::minMaxIdx(Y_tmp, &minY, &maxY);
//					if (j >= minX &&
//						j <= maxX &&
//						i >= minY &&
//						i <= maxY
//						)
//					{
//						located = true;
//						//(反距离权重法插值)
//						double d1 = sqrt((j - DTM_mapped_X.at<double>(ii, jj)) * (j - DTM_mapped_X.at<double>(ii, jj)) +
//							(i - DTM_mapped_Y.at<double>(ii, jj)) * (i - DTM_mapped_Y.at<double>(ii, jj)));
//						double d2 = sqrt((j - DTM_mapped_X.at<double>(ii + 1, jj)) * (j - DTM_mapped_X.at<double>(ii + 1, jj)) +
//							(i - DTM_mapped_Y.at<double>(ii + 1, jj)) * (i - DTM_mapped_Y.at<double>(ii + 1, jj)));
//						double d3 = sqrt((j - DTM_mapped_X.at<double>(ii, jj + 1)) * (j - DTM_mapped_X.at<double>(ii, jj + 1)) +
//							(i - DTM_mapped_Y.at<double>(ii, jj + 1)) * (i - DTM_mapped_Y.at<double>(ii, jj + 1)));
//						double d4 = sqrt((j - DTM_mapped_X.at<double>(ii + 1, jj + 1)) * (j - DTM_mapped_X.at<double>(ii + 1, jj + 1)) +
//							(i - DTM_mapped_Y.at<double>(ii + 1, jj + 1)) * (i - DTM_mapped_Y.at<double>(ii + 1, jj + 1)));
//
//						mapped_DTM.at<double>(j, i) = (DTM.at<double>(ii, jj) / d1 + DTM.at<double>(ii + 1, jj) / d2 + DTM.at<double>(ii, jj + 1) / d3 + DTM.at<double>(ii + 1, jj + 1) / d4) /
//							(1 / d1 + 1 / d2 + 1 / d3 + 1 / d4);
//						break;
//					}
//				}
//				if (located)break;
//			}
//			count++;
//			if (count % 1 == 0)
//			{
//				printf("\r估计进度2：%lf%%", double(count) / double(end_y - start_y + 1) / double(end_x - start_x + 1) * 100.0);
//				fflush(stdout);
//			}
//		}
//		/*count++;
//		if (count % 10 == 0)
//		{
//			printf("\r估计进度2：%lf%%", double(count) / double(end_y - start_y + 1) * 100.0);
//			fflush(stdout);
//		}*/
//	}
	

	return 0;
}

int Utils::geo_transformation(
	const char* grille_file,
	Mat DTM,
	double lon_upleft,
	double lat_upleft,
	double lon_interval,
	double lat_interval,
	Mat& prior_DTM,
	Mat& mapped_DTM,
	int SAR_extent_x,
	int SAR_extent_y
)
{
//	if (!grille_file || DTM.empty() || fabs(lon_upleft) > 180.0 || fabs(lat_upleft) > 90.0)
//	{
//		fprintf(stderr, "geo_transformation(): input check failed!\n");
//		return -1;
//	}
//	int ret;
//	if (DTM.type() != CV_64F)
//	{
//		DTM.convertTo(DTM, CV_64F);
//	}
//	FormatConversion conversion;
//	Mat row_matrix, col_matrix, utm_x, utm_y;
//	vector<Mat> lon_matrix;
//	vector<Mat> lat_matrix;
//	vector<double> height_vector;
//	mapped_DTM.create(SAR_extent_y, SAR_extent_x, CV_64F);
//	mapped_DTM = -1.0;
//	ret = read_grille(grille_file, row_matrix, col_matrix, lon_matrix, lat_matrix, height_vector);
//	if (return_check(ret, "read_grille()", error_head)) return -1;
//	//cvmat2bin("D:\\working_dir\\projects\\software\\InSAR\\bin\\row_matrix.bin", row_matrix);
//	//cvmat2bin("D:\\working_dir\\projects\\software\\InSAR\\bin\\col_matrix.bin", col_matrix);
//	//经纬度转UTM
//	for (int i = 0; i < height_vector.size(); i++)
//	{
//		ret = lonlat2utm(lon_matrix[i], lat_matrix[i], utm_x, utm_y);
//		if (return_check(ret, "lonlat2utm()", error_head)) return -1;
//		//conversion.creat_new_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5");
//		//conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5", "X", utm_x);
//		//conversion.write_array_to_h5("D:\\working_dir\\projects\\software\\InSAR\\bin\\utm_x.h5", "Y", utm_y);
//		utm_x.copyTo(lon_matrix[i]);
//		utm_y.copyTo(lat_matrix[i]);
//	}
//
//	
//
//	int DTM_rows = DTM.rows;
//	int DTM_cols = DTM.cols;
//	//DTM逐点转换
//	Mat DTM_mapped_X, DTM_mapped_Y;
//	DTM_mapped_X = Mat::zeros(DTM_rows, DTM_cols, CV_64F);
//	DTM_mapped_Y = Mat::zeros(DTM_rows, DTM_cols, CV_64F);
//	DTM_mapped_X = -1;
//	DTM_mapped_Y = -1;
//	volatile int count = 0;
//	//Mat mask = Mat::zeros(DTM_rows, DTM_cols, CV_8UC1);
//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < SAR_extent_y; i++)
//	{
//		const char* path[] = { "D:\\softwarepackages\\release-1928-x64-gdal-3-3-1-mapserver-7-6-4\\bin\\proj7\\share" ,nullptr };
//		OSRSetPROJSearchPaths(path);
//		OGRSpatialReference monUtm;
//		monUtm.SetWellKnownGeogCS("WGS84");
//		monUtm.SetUTM(22, 1);
//		OGRSpatialReference monGeo;
//		monGeo.SetWellKnownGeogCS("WGS84");
//		OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&monUtm, &monGeo);
//
//		for (int j = 0; j < SAR_extent_x; j++)
//		{
//			//首先确定DTM值是否有效
//			double h = prior_DTM.at<double>(i, j);
//			//double h = 20.6;
//			if (h < -9000)
//			{
//				continue;
//			}
//			//mask.at<uchar>(i, j) = 1;
//			//定位相邻的网格层
//			int low_ix, high_ix;
//			for (int k = 0; k < height_vector.size() - 1; k++)
//			{
//				if (h >= height_vector[k] && h <= height_vector[k + 1])
//				{
//					low_ix = k;
//					high_ix = k + 1;
//					//对下层网格寻找定位点
//					bool located1 = false;
//					for (int ii = 0; ii < row_matrix.rows - 1; ii++)
//					{
//						for (int jj = 0; jj < row_matrix.cols - 1; jj++)
//						{
//							//located1 = false; located2 = false;
//							//判断该点是否在四个点中间
//							double Mx, My, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy;
//							Ax = row_matrix.at<double>(ii, jj); Ay = col_matrix.at<double>(ii, jj);
//							Bx = row_matrix.at<double>(ii + 1, jj); Dy = col_matrix.at<double>(ii, jj + 1);
//
//							if (i >= Ax &&
//								i <= Bx &&
//								j >= Ay &&
//								j <= Dy &&
//								!located1
//								)
//							{
//								located1 = true;
//								//线性插值得到在下层网格上的UTM坐标
//								//UTM_x插值
//								double UTM_x_upleft = lon_matrix[low_ix].at<double>(ii, jj);
//								double UTM_x_upright = lon_matrix[low_ix].at<double>(ii, jj + 1);
//								double UTM_x_lowleft = lon_matrix[low_ix].at<double>(ii + 1, jj);
//								double UTM_x_lowright = lon_matrix[low_ix].at<double>(ii + 1, jj + 1);
//								double upper = UTM_x_upleft + (UTM_x_upright - UTM_x_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
//									(j - col_matrix.at<double>(ii, jj));
//								double lower = UTM_x_lowleft + (UTM_x_lowright - UTM_x_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
//									(j - col_matrix.at<double>(ii + 1, jj));
//								double UTM_x_final_lower = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
//									(i - row_matrix.at<double>(ii + 1, jj));
//
//								//UTM_y插值
//								double UTM_y_upleft = lat_matrix[low_ix].at<double>(ii, jj);
//								double UTM_y_upright = lat_matrix[low_ix].at<double>(ii, jj + 1);
//								double UTM_y_lowleft = lat_matrix[low_ix].at<double>(ii + 1, jj);
//								double UTM_y_lowright = lat_matrix[low_ix].at<double>(ii + 1, jj + 1);
//								upper = UTM_y_upleft + (UTM_y_upright - UTM_y_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
//									(j - col_matrix.at<double>(ii, jj));
//								lower = UTM_y_lowleft + (UTM_y_lowright - UTM_y_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
//									(j - col_matrix.at<double>(ii + 1, jj));
//								double UTM_y_final_lower = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
//									(i - row_matrix.at<double>(ii + 1, jj));
//
//								//线性插值得到在上层网格上的UTM坐标
//								//UTM_x插值
//								UTM_x_upleft = lon_matrix[high_ix].at<double>(ii, jj);
//								UTM_x_upright = lon_matrix[high_ix].at<double>(ii, jj + 1);
//								UTM_x_lowleft = lon_matrix[high_ix].at<double>(ii + 1, jj);
//								UTM_x_lowright = lon_matrix[high_ix].at<double>(ii + 1, jj + 1);
//								upper = UTM_x_upleft + (UTM_x_upright - UTM_x_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
//									(j - col_matrix.at<double>(ii, jj));
//								lower = UTM_x_lowleft + (UTM_x_lowright - UTM_x_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
//									(j - col_matrix.at<double>(ii + 1, jj));
//								double UTM_x_final_higher = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
//									(i - row_matrix.at<double>(ii + 1, jj));
//
//								//UTM_y插值
//								UTM_y_upleft = lat_matrix[high_ix].at<double>(ii, jj);
//								UTM_y_upright = lat_matrix[high_ix].at<double>(ii, jj + 1);
//								UTM_y_lowleft = lat_matrix[high_ix].at<double>(ii + 1, jj);
//								UTM_y_lowright = lat_matrix[high_ix].at<double>(ii + 1, jj + 1);
//								upper = UTM_y_upleft + (UTM_y_upright - UTM_y_upleft) / (col_matrix.at<double>(ii, jj + 1) - col_matrix.at<double>(ii, jj)) *
//									(j - col_matrix.at<double>(ii, jj));
//								lower = UTM_y_lowleft + (UTM_y_lowright - UTM_y_lowleft) / (col_matrix.at<double>(ii + 1, jj + 1) - col_matrix.at<double>(ii + 1, jj)) *
//									(j - col_matrix.at<double>(ii + 1, jj));
//								double UTM_y_final_higher = lower + (upper - lower) / (row_matrix.at<double>(ii, jj) - row_matrix.at<double>(ii + 1, jj)) *
//									(i - row_matrix.at<double>(ii + 1, jj));
//
//
//
//								//上下两层之间插值得到UTM_x和UTM_y
//								double UTM_x_final = UTM_x_final_lower + (UTM_x_final_higher - UTM_x_final_lower) / (height_vector[high_ix] - height_vector[low_ix]) *
//									(h - height_vector[low_ix]);
//								double UTM_y_final = UTM_y_final_lower + (UTM_y_final_higher - UTM_y_final_lower) / (height_vector[high_ix] - height_vector[low_ix]) *
//									(h - height_vector[low_ix]);
//
//								double lat_x, lon_y;
//								lat_x = UTM_x_final;
//								lon_y = UTM_y_final;
//								int reprojected = coordTrans->Transform(1, &lat_x, &lon_y);
//
//
//								//通过插值得到的lat_x和lon_y再次插值得到DTM
//								if (lat_x <= lat_upleft &&
//									lat_x >= (lat_upleft - (DTM_rows - 1)*lat_interval) &&
//									lon_y <= (lon_upleft + (DTM_cols - 1)*lon_interval) &&
//									lon_y >= lon_upleft
//									)
//								{
//									int row = floor((lat_upleft - lat_x) / lat_interval);
//									double delta_row = ((lat_upleft - lat_x) / lat_interval - row);
//									int col = floor((lon_y - lon_upleft) / lon_interval);
//									double delta_col = ((lon_y - lon_upleft) / lon_interval - col);
//									if (row < 0 || row >= DTM_rows - 1 || col < 0 || col >= DTM_cols - 1) continue;
//									if (DTM.at<double>(row, col) > -1000 &&
//										DTM.at<double>(row + 1, col) > -1000 &&
//										DTM.at<double>(row + 1, col + 1) > -1000 &&
//										DTM.at<double>(row, col + 1) > -1000
//										)
//									{
//										upper = DTM.at<double>(row, col) + (DTM.at<double>(row, col + 1) - DTM.at<double>(row, col)) * delta_row;
//										lower = DTM.at<double>(row + 1, col) + (DTM.at<double>(row + 1, col + 1) - DTM.at<double>(row + 1, col)) * delta_row;
//										mapped_DTM.at<double>(i, j) = lower + (upper - lower) * delta_col;
//									}
//
//								}
//								break;
//							}
//						}
//						if (located1)break;
//					}
//				}
//			}
//		}
//		count++;
//		if (count % 10 == 0)
//		{
//			printf("\r估计进度1：%lf%%", double(count) / double(SAR_extent_y) * 100.0);
//			fflush(stdout);
//		}
//	}
	return 0;
}

int Utils::read_grille(
	const char* grille_file,
	Mat& row_matrix, 
	Mat& col_matrix,
	vector<Mat>& lon_matrix,
	vector<Mat>& lat_matrix, 
	vector<double>& height_vector
)
{
	if (!grille_file)
	{
		fprintf(stderr, "read_grille(): input check failed!\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(grille_file, "rt");
	if (!fp)
	{
		fprintf(stderr, "read_grille(): cannot open grille file!\n");
		return -1;
	}
	char* data = (char*)malloc(1024);
	char* ptr;
	double val;
	fgets(data, 1024, fp);
	fgets(data, 1024, fp);
	fgets(data, 1024, fp);
	memset(data, 0, 1024);
	fgets(data, 1024, fp);
	int rows = strtol(data, &ptr, 10);
	fgets(data, 1024, fp);
	int cols = strtol(data, &ptr, 10);
	fgets(data, 1024, fp);
	int height_num = strtol(data, &ptr, 10);
	size_t total = rows * cols * height_num;
	Mat zeroMat = Mat::zeros(rows, cols, CV_64F);
	zeroMat.copyTo(row_matrix);
	zeroMat.copyTo(col_matrix);
	lon_matrix.resize(height_num);
	lat_matrix.resize(height_num);
	height_vector.resize(height_num);
	for (int i = 0; i < height_num; i++)
	{
		zeroMat.copyTo(lon_matrix[i]);
		zeroMat.copyTo(lat_matrix[i]);
	}
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			for (int k = 0; k < height_num; k++)
			{
				fgets(data, 1024, fp);
				row_matrix.at<double>(i, j) = strtod(data, &ptr);
				col_matrix.at<double>(i, j) = strtod(ptr, &ptr);
				height_vector[k] = strtod(ptr, &ptr);
				lon_matrix[k].at<double>(i, j) = strtod(ptr, &ptr);
				lat_matrix[k].at<double>(i, j) = strtod(ptr, &ptr);
			}
		}
	}
	fclose(fp);
	free(data);
	return 0;
}

int Utils::lonlat2utm(Mat lon, Mat lat, Mat& UTM_X, Mat& UTM_Y)
{
	//if (lon.empty() || lat.empty() || lon.size() != lat.size() || lon.type() != lat.type() || lon.type() != CV_64F)
	//{
	//	fprintf(stderr, "lonlat2utm(): input check failed!\n");
	//	return -1;
	//}
	//const char* path[] = { "D:\\softwarepackages\\release-1928-x64-gdal-3-3-1-mapserver-7-6-4\\bin\\proj7\\share" ,nullptr };
	//OSRSetPROJSearchPaths(path);
	//OGRSpatialReference monUtm;
	//monUtm.SetWellKnownGeogCS("WGS84");
	//monUtm.SetUTM(22, true);

	//OGRSpatialReference monGeo;
	//monGeo.SetWellKnownGeogCS("WGS84");

	//OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&monGeo, &monUtm);

	//int nr = lon.rows;
	//int nc = lon.cols;
	//lon.copyTo(UTM_X);
	//lon.copyTo(UTM_Y);
	//for (int i = 0; i < nr; i++)
	//{
	//	for (int j = 0; j < nc; j++)
	//	{
	//		double x, y;
	//		x = lat.at<double>(i, j);
	//		y = lon.at<double>(i, j);
	//		int reprojected = coordTrans->Transform(1, &x, &y);
	//		UTM_X.at<double>(i, j) = x;
	//		UTM_Y.at<double>(i, j) = y;
	//	}
	//}
	///*int reprojected = coordTrans->Transform(1, &x, &y);*/
	//// If OK, print the coords.

	////delete coordTrans;
	////coordTrans = OGRCreateCoordinateTransformation(&monGeo, &monUtm);
	////reprojected = coordTrans->Transform(1, &x, &y);

	//// If OK, Print the coords.
	//delete coordTrans;
	return 0;
}

int Utils::read_LVIS(
	vector<string>& LVIS2_filelist,
	Mat& DTM,
	Mat& DSM,
	Mat& RH100,
	Mat& RH95,
	Mat& lat,
	Mat& lon
)
{
	if (LVIS2_filelist.empty())
	{
		fprintf(stderr, "read_LVIS(): input check failed!\n");
		return -1;
	}
	int n_files = LVIS2_filelist.size();
	int LVIS_cols = 42;
	vector<double> dtm, dsm, rh100, rh95, lo, la;
	FILE* fp = NULL;
	char str[1024];
	char* ret = NULL;
	char* ptr = NULL;
	double val;
	for (int i = 0; i < n_files; i++)
	{
		fp = fopen(LVIS2_filelist[i].c_str(), "rt");
		if (!fp)
		{
			fprintf(stderr, "read_LVIS(): cannot open %s!\n", LVIS2_filelist[i].c_str());
			return -1;
		}
		for (int j = 0; j < 18; j++)
		{
			ret = fgets(str, 1024, fp);
		}
		while (ret)
		{
			ret=fgets(str, 1024, fp);
			vector<double> tmp;
			ptr = str;
			for (int k = 0; k < LVIS_cols; k++)
			{
				val = strtod(ptr, &ptr);
				tmp.push_back(val);
			}
			lo.push_back(tmp[3]);
			la.push_back(tmp[4]);
			dtm.push_back(tmp[5]);
			dsm.push_back(tmp[11]);
			rh100.push_back(tmp[33]);
			rh95.push_back(tmp[28]);
			tmp.clear();
		}
		fclose(fp);
		fp = NULL;
	}
	DTM.create(dtm.size(), 1, CV_64F);
	DTM = 0;
	DTM.copyTo(DSM);
	DTM.copyTo(lon);
	DTM.copyTo(lat);
	DTM.copyTo(RH100);
	DTM.copyTo(RH95);
	for (int i = 0; i < dtm.size(); i++)
	{
		DTM.at<double>(i, 0) = dtm[i];
		DSM.at<double>(i, 0) = dsm[i];
		lon.at<double>(i, 0) = lo[i];
		lat.at<double>(i, 0) = la[i];
		RH100.at<double>(i, 0) = rh100[i];
		RH95.at<double>(i, 0) = rh95[i];
	}
	return 0;
}

int Utils::geo2sar_DLR(
	double east_min,
	double north_min,
	double east_max,
	double north_max,
	int projection_zone,
	double pixel_spacing,
	Mat& sr2geo_az,
	Mat& sr2geo_rg,
	Mat& sr2geo_h_ref,
	Mat& sr2geo3d_rg_o1,
	Mat& sr2geo3d_rg_o2,
	Mat& sr2geo3d_az_o1,
	Mat& sr2geo3d_az_o2,
	Mat& DTM,
	Mat& DTM_lon,
	Mat& DTM_lat,
	Mat& mapped_slc_rg,
	Mat& mapped_slc_az,
	Mat& mapped_h0
)
{
	//if (pixel_spacing < 0 ||
	//	sr2geo_az.empty() || sr2geo_az.type() != CV_64F ||
	//	sr2geo_az.size() != sr2geo_rg.size() || sr2geo_az.size() != sr2geo_h_ref.size() || sr2geo_az.size() != sr2geo3d_rg_o1.size() || 
	//	sr2geo_az.size() != sr2geo3d_rg_o2.size() || sr2geo_az.size() != sr2geo3d_az_o1.size() || sr2geo_az.size() != sr2geo3d_az_o2.size() ||
	//	DTM.empty() || DTM.type() != CV_64F ||
	//	DTM.size() != DTM_lon.size() || DTM.size() != DTM_lat.size()
	//	)
	//{
	//	fprintf(stderr, "geo2sar_DLR(): input check failed!\n");
	//	return -1;
	//}
	////确定场景的经纬度范围
	//Mat lon, lat;
	//lon.create(4, 1, CV_64F); lon = 0.0; lon.copyTo(lat);
	//double utm_x, utm_y;
	//const char* path[] = { "D:\\softwarepackages\\release-1928-x64-gdal-3-3-1-mapserver-7-6-4\\bin\\proj7\\share" ,nullptr };
	//OSRSetPROJSearchPaths(path);
	//OGRSpatialReference monUtm;
	//monUtm.SetWellKnownGeogCS("WGS84");
	//monUtm.SetUTM(abs(projection_zone), projection_zone>0);
	//OGRSpatialReference monGeo;
	//monGeo.SetWellKnownGeogCS("WGS84");
	//OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&monUtm, &monGeo);
	//double x, y;
	//x = east_min; y = north_min;
	//int reprojected = coordTrans->Transform(1, &x, &y);
	//lon.at<double>(0, 0) = y; lat.at<double>(0, 0) = x;

	//x = east_min; y = north_max;
	//reprojected = coordTrans->Transform(1, &x, &y);
	//lon.at<double>(1, 0) = y; lat.at<double>(1, 0) = x;

	//x = east_max; y = north_min;
	//reprojected = coordTrans->Transform(1, &x, &y);
	//lon.at<double>(2, 0) = y; lat.at<double>(2, 0) = x;

	//x = east_max; y = north_max;
	//reprojected = coordTrans->Transform(1, &x, &y);
	//lon.at<double>(3, 0) = y; lat.at<double>(3, 0) = x;

	//double lon_min, lon_max, lat_min, lat_max;
	//cv::minMaxLoc(lon, &lon_min, &lon_max);
	//cv::minMaxLoc(lat, &lat_min, &lat_max);
	//delete coordTrans;

	//
	//mapped_slc_rg.create(DTM.rows, DTM.cols, CV_64F); mapped_slc_rg = -1;
	//mapped_slc_az.create(DTM.rows, DTM.cols, CV_64F); mapped_slc_az = -1;
	//mapped_h0.create(DTM.rows, DTM.cols, CV_64F); mapped_h0 = -1;
	//Mat UTM_x, UTM_y;
	//DTM_lat.copyTo(UTM_x);
	//DTM_lon.copyTo(UTM_y);
	//OGRCoordinateTransformation* coordTrans2 = OGRCreateCoordinateTransformation(&monGeo , &monUtm);
	////将DTM的经纬度转换为UTM
	//for (int i = 0; i < DTM.rows; i++)
	//{
	//	for (int j = 0; j < DTM.cols; j++)
	//	{
	//		if (DTM.at<double>(i, j) < -9000.0) continue;
	//		double xx, yy;
	//		xx = DTM_lat.at<double>(i, j);
	//		yy = DTM_lon.at<double>(i, j);
	//		int reprojected = coordTrans2->Transform(1, &xx, &yy);
	//		UTM_x.at<double>(i, j) = xx;
	//		UTM_y.at<double>(i, j) = yy;
	//	}
	//}
	//delete coordTrans2;
	////开始投影
	//for (int i = 0; i < DTM.rows; i++)
	//{
	//	for (int j = 0; j < DTM.cols; j++)
	//	{
	//		if (DTM.at<double>(i, j) < -9000.0) continue;
	//		double utm_x, utm_y, h;
	//		utm_x = UTM_x.at<double>(i, j);
	//		utm_y = UTM_y.at<double>(i, j);
	//		h = DTM.at<double>(i, j);
	//		if (DTM_lon.at<double>(i, j) <= lon_min || DTM_lon.at<double>(i, j) >= lon_max ||
	//			DTM_lat.at<double>(i, j) <= lat_min || DTM_lat.at<double>(i, j) >= lat_max
	//			)
	//		{
	//			continue;
	//		}
	//		//插值得到2D地理编码下的距离方位坐标，参考DEM和3D地理编码系数
	//		double easting, northing, rg0, az0, h0, rg_o1, rg_o2, az_o1, az_o2, upper, lower, left, right;
	//		int row, col;
	//		easting = (utm_x - east_min) / pixel_spacing;
	//		northing = (utm_y - north_min) / pixel_spacing;

	//		row = sr2geo_h_ref.rows - 1 - int(floor(northing));
	//		col = int(floor(easting));
	//		if (row < 1 || row > sr2geo_h_ref.rows - 1 || col < 0 || col >= sr2geo_h_ref.cols - 1)
	//		{
	//			continue;
	//		}
	//		if (sr2geo_h_ref.at<double>(row, col) < -9000 ||
	//			sr2geo_h_ref.at<double>(row, col + 1) < -9000 ||
	//			sr2geo_h_ref.at<double>(row - 1, col) < -9000 ||
	//			sr2geo_h_ref.at<double>(row - 1, col + 1) < -9000
	//			)
	//		{
	//			continue;
	//		}
	//			


	//		lower = sr2geo_rg.at<double>(row, col) + (easting - double(col)) * (sr2geo_rg.at<double>(row, col + 1) - sr2geo_rg.at<double>(row, col));
	//		upper = sr2geo_rg.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo_rg.at<double>(row - 1, col + 1) - sr2geo_rg.at<double>(row - 1, col));
	//		rg0 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo_az.at<double>(row, col) + (easting - double(col)) * (sr2geo_az.at<double>(row, col + 1) - sr2geo_az.at<double>(row, col));
	//		upper = sr2geo_az.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo_az.at<double>(row - 1, col + 1) - sr2geo_az.at<double>(row - 1, col));
	//		az0 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo_h_ref.at<double>(row, col) + (easting - double(col)) * (sr2geo_h_ref.at<double>(row, col + 1) - sr2geo_h_ref.at<double>(row, col));
	//		upper = sr2geo_h_ref.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo_h_ref.at<double>(row - 1, col + 1) - sr2geo_h_ref.at<double>(row - 1, col));
	//		h0 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo3d_rg_o1.at<double>(row, col) + (easting - double(col)) * (sr2geo3d_rg_o1.at<double>(row, col + 1) - sr2geo3d_rg_o1.at<double>(row, col));
	//		upper = sr2geo3d_rg_o1.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo3d_rg_o1.at<double>(row - 1, col + 1) - sr2geo3d_rg_o1.at<double>(row - 1, col));
	//		rg_o1 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo3d_rg_o2.at<double>(row, col) + (easting - double(col)) * (sr2geo3d_rg_o2.at<double>(row, col + 1) - sr2geo3d_rg_o2.at<double>(row, col));
	//		upper = sr2geo3d_rg_o2.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo3d_rg_o2.at<double>(row - 1, col + 1) - sr2geo3d_rg_o2.at<double>(row - 1, col));
	//		rg_o2 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo3d_az_o1.at<double>(row, col) + (easting - double(col)) * (sr2geo3d_az_o1.at<double>(row, col + 1) - sr2geo3d_az_o1.at<double>(row, col));
	//		upper = sr2geo3d_az_o1.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo3d_az_o1.at<double>(row - 1, col + 1) - sr2geo3d_az_o1.at<double>(row - 1, col));
	//		az_o1 = lower + (upper - lower) * (northing - floor(northing));

	//		lower = sr2geo3d_az_o2.at<double>(row, col) + (easting - double(col)) * (sr2geo3d_az_o2.at<double>(row, col + 1) - sr2geo3d_az_o2.at<double>(row, col));
	//		upper = sr2geo3d_az_o2.at<double>(row - 1, col) + (easting - double(col)) * (sr2geo3d_az_o2.at<double>(row - 1, col + 1) - sr2geo3d_az_o2.at<double>(row - 1, col));
	//		az_o2 = lower + (upper - lower) * (northing - floor(northing));


	//		//left = sr2geo_rg.at<double>(row, col) + (northing - floor(northing)) * (sr2geo_rg.at<double>(row - 1, col) - sr2geo_rg.at<double>(row, col));
	//		//right = sr2geo_rg.at<double>(row, col + 1) + (northing - floor(northing)) * (sr2geo_rg.at<double>(row - 1, col + 1) - sr2geo_rg.at<double>(row, col + 1));
	//		//rg0 = left + (right - left) * (easting - double(col));

	//		//left = sr2geo_az.at<double>(row, col) + (northing - floor(northing)) * (sr2geo_az.at<double>(row - 1, col) - sr2geo_az.at<double>(row, col));
	//		//right = sr2geo_az.at<double>(row, col + 1) + (northing - floor(northing)) * (sr2geo_az.at<double>(row - 1, col + 1) - sr2geo_az.at<double>(row, col + 1));
	//		//az0 = left + (right - left) * (easting - double(col));

	//		//left = sr2geo_h_ref.at<double>(row, col) + (northing - floor(northing)) * (sr2geo_h_ref.at<double>(row - 1, col) - sr2geo_h_ref.at<double>(row, col));
	//		//right = sr2geo_h_ref.at<double>(row, col + 1) + (northing - floor(northing)) * (sr2geo_h_ref.at<double>(row - 1, col + 1) - sr2geo_h_ref.at<double>(row, col + 1));
	//		//h0 = left + (right - left) * (easting - double(col));

	//		double delta_h, rg_new, az_new;
	//		delta_h = h - h0;
	//		rg_new = rg0 + rg_o1 * delta_h + rg_o2 * delta_h * delta_h;
	//		az_new = az0 + az_o1 * delta_h + az_o2 * delta_h * delta_h;
	//
	//		mapped_slc_rg.at<double>(i, j) = rg0;
	//		mapped_slc_az.at<double>(i, j) = az0;
	//		mapped_h0.at<double>(i, j) = h0;
	//	}
	//}
	return 0;
}

// 根据经纬度获取大地水准面高差
double Utils::getGeoidHeight(const std::string& geoidFilePath, double lon, double lat) {
	// 注册 GDAL 驱动
	GDALAllRegister();

	/* 打开 Geoid 文件 */
	GDALDatasetH hDataset = GDALOpen(geoidFilePath.c_str(), GA_ReadOnly);
	if (hDataset == NULL)
	{
		fprintf(stderr,
			"无法打开 Geoid 文件: %s\n",
			geoidFilePath);
		return 0.0;
	}

	/* 获取第一个波段 */
	GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);
	if (hBand == NULL)
	{
		fprintf(stderr, "无法获取波段数据\n");
		GDALClose(hDataset);
		return 0.0;
	}

	/* 获取地理变换 */
	double adfGeoTransform[6];
	if (GDALGetGeoTransform(hDataset, adfGeoTransform) != CE_None)
	{
		fprintf(stderr, "无法获取地理变换信息\n");
		GDALClose(hDataset);
		return 0.0;
	}

	/* 经纬度 → 像素坐标 */
	double x = (lon - adfGeoTransform[0]) / adfGeoTransform[1];
	double y = (lat - adfGeoTransform[3]) / adfGeoTransform[5];

	/* 读取 Geoid Height */
	float geoidHeight = 0.0f;
	if (GDALRasterIO(
		hBand,
		GF_Read,
		(int)x, (int)y,
		1, 1,
		&geoidHeight,
		1, 1,
		GDT_Float32,
		0, 0) != CE_None)
	{
		fprintf(stderr, "无法读取 Geoid 数据\n");
		GDALClose(hDataset);
		return 0.0;
	}

	GDALClose(hDataset);
	return (double)geoidHeight;
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






