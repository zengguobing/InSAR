// Registration.cpp : 定义 DLL 应用程序的导出函数。
//

#include<string.h>
#include "stdafx.h"
#include<math.h>
#include"..\include\Registration.h"
#include"..\include\FormatConversion.h"
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

Registration::Registration()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "REGISTRATION_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "REGISTRATION_DLL_ERROR: error happens when using parallel computing in function: ");
}

Registration::~Registration()
{
}
int Registration::fft2(Mat& Src, Mat& Dst)
{
	if (Src.rows < 1 ||
		Src.cols < 1 ||
		Src.channels() != 1||
		Src.type() != CV_64F)
	{
		fprintf(stderr, "fft2(): input check failed!\n\n");
		return -1;
	}
	Mat planes[] = { Mat_<double>(Src), Mat::zeros(Src.size(), CV_64F) };
	Mat complexImg;
	merge(planes, 2, complexImg);
	dft(complexImg, Dst);
	return 0;
}

int Registration::fftshift2(Mat& matrix)
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

int Registration::real_coherent(ComplexMat& Master, ComplexMat& Slave, int* offset_row, int* offset_col)
{
	if (Master.GetRows() < 1 ||
		Master.GetCols() < 1 ||
		Master.GetRows() != Slave.GetRows() ||
		Master.GetCols() != Slave.GetCols())
	{
		fprintf(stderr, "real_coherent(): input check failed!\n\n");
		return -1;
	}
	int ret;
	Mat img1;
	Mat img2;
	img1 = Master.GetMod();
	img2 = Slave.GetMod();
	Mat im1fft;
	Mat im2fft;
	ret = fft2(img1, im1fft);
	if (return_check(ret, "fft2(*, *)", error_head)) return -1;
	ret = fft2(img2, im2fft);
	if (return_check(ret, "fft2(*, *)", error_head)) return -1;
	Mat spectrum;
	mulSpectrums(im1fft, im2fft, spectrum, 0, true);

	Mat result;
	idft(spectrum, result, DFT_REAL_OUTPUT);//需要显示图像时可以用DFT_SCALE

	ret = fftshift2(result);
	if (return_check(ret, "fftshift2(*)", error_head)) return -1;
	normalize(result, result, 0, 1, NORM_MINMAX);

	int r = result.rows / 2;
	int c = result.cols / 2;
	Point peak_loc;
	minMaxLoc(result, NULL, NULL, NULL, &peak_loc);

	*offset_row = r - peak_loc.y;
	*offset_col = c - peak_loc.x;
	return 0;
}

int Registration::registration_pixel(ComplexMat& Master, ComplexMat& Slave, int* move_r, int* move_c)
{
	if (Master.GetRows() < 1 ||
		Master.GetCols() < 1 ||
		Master.GetRows() != Slave.GetRows() ||
		Master.GetCols() != Slave.GetCols())
	{
		fprintf(stderr, "registration_pixel(): input check failed!\n\n");
		return -1;
	}
	int offset_rows, offset_cols, ret;
	offset_cols = 0;
	offset_rows = 0;
	ret = real_coherent(Master, Slave, &offset_rows, &offset_cols);//相关函数求取偏移量
	if (move_r)*move_r = offset_rows;
	if (move_c)*move_c = offset_cols;
	if (return_check(ret, "real_coherent(*, *, *, *)", error_head)) return -1;
	////搬移与裁剪

	//行偏移（竖直移动）
	ComplexMat image_master_mid;
	ComplexMat image_slave_mid;
	int nr = Slave.GetRows();
	int nc = Slave.GetCols();
	//////////////////////////////检查粗配准偏移量是否超过图像大小////////////////////////////////////
	if ((offset_rows > 0 ? offset_rows : -offset_rows) >= nr ||
		(offset_cols > 0 ? offset_cols : -offset_cols) >= nc)
	{
		fprintf(stderr, "registration_pixel(): registration offset exceed size of images!\n\n");
		return -1;
	}
	if (offset_rows >= 0)
	{
		//实部
		Slave.re(Range(offset_rows, nr), Range(0, nc)).copyTo(image_slave_mid.re);//辅图像向上搬移

		Master.re(Range(0, nr - offset_rows), Range(0, nc)).copyTo(image_master_mid.re);//裁剪主图像


		//虚部
		Slave.im(Range(offset_rows, nr), Range(0, nc)).copyTo(image_slave_mid.im);//辅图像向上搬移

		Master.im(Range(0, nr - offset_rows), Range(0, nc)).copyTo(image_master_mid.im);//裁剪

	}

	else
	{
		//实部
		Slave.re(Range(0, nr + offset_rows), Range(0, nc)).copyTo(image_slave_mid.re);//辅图像向下搬移

		Master.re(Range(-offset_rows, nr), Range(0, nc)).copyTo(image_master_mid.re);//裁剪


		//虚部
		Slave.im(Range(0, nr + offset_rows), Range(0, nc)).copyTo(image_slave_mid.im);//辅图像向下搬移

		Master.im(Range(-offset_rows, nr), Range(0, nc)).copyTo(image_master_mid.im);//裁剪

	}

	ComplexMat image_master_regis;
	ComplexMat image_slave_regis;
	int nr1 = image_master_mid.re.rows;
	int nc1 = image_master_mid.re.cols;
	//列偏移（水平移动）
	if (offset_cols >= 0)//辅图像向左搬移
	{
		//实部
		image_slave_mid.re(Range(0, nr1), Range(offset_cols, nc1)).copyTo(image_slave_regis.re);

		image_master_mid.re(Range(0, nr1), Range(0, nc1 - offset_cols)).copyTo(image_master_regis.re);

		//虚部
		image_slave_mid.im(Range(0, nr1), Range(offset_cols, nc1)).copyTo(image_slave_regis.im);

		image_master_mid.im(Range(0, nr1), Range(0, nc1 - offset_cols)).copyTo(image_master_regis.im);

	}

	else//辅图像向右搬移
	{
		//实部
		image_slave_mid.re(Range(0, nr1), Range(0, nc1 + offset_cols)).copyTo(image_slave_regis.re);

		image_master_mid.re(Range(0, nr1), Range(-offset_cols, nc1)).copyTo(image_master_regis.re);

		//虚部
		image_slave_mid.im(Range(0, nr1), Range(0, nc1 + offset_cols)).copyTo(image_slave_regis.im);

		image_master_mid.im(Range(0, nr1), Range(-offset_cols, nc1)).copyTo(image_master_regis.im);
	}

	Master.re = image_master_regis.re;

	Master.im = image_master_regis.im;

	Slave.re = image_slave_regis.re;

	Slave.im = image_slave_regis.im;
	return 0;
}

int Registration::interp_paddingzero(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, int interp_times)
{
	if (InputMatrix.GetRows() < 2 ||
		InputMatrix.GetCols() < 2 ||
		interp_times < 2)
	{
		fprintf(stderr,"interp_paddingzero(): input check failed!\n\n");
		return -1;
	}

	int nr = InputMatrix.GetRows();
	int nc = InputMatrix.GetCols();

	OutputMatrix.re = Mat::zeros(interp_times * nr, interp_times * nc, CV_64F);

	OutputMatrix.im = Mat::zeros(interp_times * nr, interp_times * nc, CV_64F);
	Mat re, im;
	InputMatrix.re.copyTo(re);
	InputMatrix.im.copyTo(im);
	Mat planes[] = { Mat_<double>(re), Mat_<double>(im) };

	Mat complexImg;

	merge(planes, 2, complexImg);

	dft(complexImg, complexImg, DFT_COMPLEX_OUTPUT);

	split(complexImg, planes);

	planes[0](Range(0, nr / 2), Range(0, nc / 2)).copyTo(OutputMatrix.re(Range(0, nr / 2), Range(0, nc / 2)));
	planes[1](Range(0, nr / 2), Range(0, nc / 2)).copyTo(OutputMatrix.im(Range(0, nr / 2), Range(0, nc / 2)));

	planes[0](Range(nr / 2, nr), Range(0, nc / 2)).copyTo(OutputMatrix.re(Range(nr * interp_times - nr / 2, nr * interp_times), Range(0, nc / 2)));
	planes[1](Range(nr / 2, nr), Range(0, nc / 2)).copyTo(OutputMatrix.im(Range(nr * interp_times - nr / 2, nr * interp_times), Range(0, nc / 2)));

	planes[0](Range(0, nr / 2), Range(nc / 2, nc)).copyTo(OutputMatrix.re(Range(0, nr / 2), Range(nc * interp_times - nc / 2, nc * interp_times)));
	planes[1](Range(0, nr / 2), Range(nc / 2, nc)).copyTo(OutputMatrix.im(Range(0, nr / 2), Range(nc * interp_times - nc / 2, nc * interp_times)));

	planes[0](Range(nr / 2, nr), Range(nc / 2, nc)).copyTo(OutputMatrix.re(Range(nr * interp_times - nr / 2, nr * interp_times), Range(nc * interp_times - nc / 2, nc * interp_times)));
	planes[1](Range(nr / 2, nr), Range(nc / 2, nc)).copyTo(OutputMatrix.im(Range(nr * interp_times - nr / 2, nr * interp_times), Range(nc * interp_times - nc / 2, nc * interp_times)));

	Mat planes1[] = { Mat_<double>(OutputMatrix.re), Mat_<double>(OutputMatrix.im) };

	merge(planes1, 2, complexImg);

	idft(complexImg, complexImg);

	split(complexImg, planes1);

	OutputMatrix.re = planes1[0];
	OutputMatrix.im = planes1[1];
	return 0;
}

int Registration::interp_cubic(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, double offset_row, double offset_col)
{
	int nr = InputMatrix.GetRows();
	int nc = InputMatrix.GetCols();//输入矩阵尺寸
	if (nr < 2 || nc < 2 || InputMatrix.re.type() != CV_64F)
	{
		fprintf(stderr, "interp_cubic(): input check failed!\n\n");
		return -1;
	}
	ComplexMat new_image_slave;
	new_image_slave.re = Mat::zeros(nr + 3, nc + 3, CV_64F);
	new_image_slave.im = Mat::zeros(nr + 3, nc + 3, CV_64F);

	//扩充矩阵(扩展三行三列)

	//实部

	InputMatrix.re(Range(0, 1), Range(0, 1)).copyTo(new_image_slave.re(Range(0, 1), Range(0, 1)));

	InputMatrix.re(Range(0, 1), Range(0, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(1, nc + 1)));

	InputMatrix.re(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(nc + 2, nc + 3)));

	InputMatrix.re(Range(0, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, nc)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(1, nc + 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(nc + 2, nc + 3)));

	new_image_slave.re(Range(nr + 1, nr + 2), Range(0, nc + 3)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(0, nc + 3)));

	InputMatrix.re(Range(0, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(nc + 2, nc + 3)));

	new_image_slave.re(Range(0, nr + 3), Range(nc + 1, nc + 2)).copyTo(new_image_slave.re(Range(0, nr + 3), Range(nc + 2, nc + 3)));



	InputMatrix.re(Range(0, nr), Range(0, nc)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(1, nc + 1)));


	//虚部

	InputMatrix.im(Range(0, 1), Range(0, 1)).copyTo(new_image_slave.im(Range(0, 1), Range(0, 1)));

	InputMatrix.im(Range(0, 1), Range(0, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(1, nc + 1)));

	InputMatrix.im(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(nc + 2, nc + 3)));

	InputMatrix.im(Range(0, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, nc)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(1, nc + 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(nc + 2, nc + 3)));

	new_image_slave.im(Range(nr + 1, nr + 2), Range(0, nc + 3)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(0, nc + 3)));

	InputMatrix.im(Range(0, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(nc + 2, nc + 3)));

	new_image_slave.im(Range(0, nr + 3), Range(nc + 1, nc + 2)).copyTo(new_image_slave.im(Range(0, nr + 3), Range(nc + 2, nc + 3)));

	//内点
	InputMatrix.im(Range(0, nr), Range(0, nc)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(1, nc + 1)));


	//行权
	double row_weight[4];
	row_weight[0] = WeightCalculation(1.0 + offset_row);
	row_weight[1] = WeightCalculation(offset_row);
	row_weight[2] = WeightCalculation(1.0 - offset_row);
	row_weight[3] = WeightCalculation(2.0 - offset_row);

	Mat Row_weight(4, 1, CV_64F, row_weight);

	//列权
	double col_weight[4];
	col_weight[0] = WeightCalculation(1.0 + offset_col);
	col_weight[1] = WeightCalculation(offset_col);
	col_weight[2] = WeightCalculation(1.0 - offset_col);
	col_weight[3] = WeightCalculation(2.0 - offset_col);

	Mat Col_weight(1, 4, CV_64F, col_weight);

	
	//权矩阵
	Mat Weight = Row_weight * Col_weight;

	//实部虚部分别插值
	Mat image_slave_regis_re = Mat::zeros(nr, nc, CV_64F);
	Mat image_slave_regis_im = Mat::zeros(nr, nc, CV_64F);

	double temp;
#pragma omp parallel for schedule(guided) \
	private(temp)
	for (int i = 0; i <= nr - 1; i++)
	{
		for (int j = 0; j <= nc - 1; j++)
		{
			temp = new_image_slave.re(Range(i, i + 4), Range(j, j + 4)).dot(Weight);

			image_slave_regis_re.at<double>(i, j) = temp;


			temp = new_image_slave.im(Range(i, i + 4), Range(j, j + 4)).dot(Weight);

			image_slave_regis_im.at<double>(i, j) = temp;
		}

	}
	image_slave_regis_re.copyTo(OutputMatrix.re);
	image_slave_regis_im.copyTo(OutputMatrix.im);
	return 0;
}

int Registration::interp_cubic(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, Mat& Coefficient)
{
	int nr = InputMatrix.GetRows();
	int nc = InputMatrix.GetCols();//输入矩阵尺寸
	if (nr < 2 || nc < 2 || InputMatrix.re.type() != CV_64F)
	{
		fprintf(stderr, "interp_cubic(): input check failed!\n\n");
		return -1;
	}

	ComplexMat new_image_slave;
	new_image_slave.re = Mat::zeros(nr + 3, nc + 3, CV_64F);
	new_image_slave.im = Mat::zeros(nr + 3, nc + 3, CV_64F);

	//扩充矩阵(扩展三行三列)

	//实部

	InputMatrix.re(Range(0, 1), Range(0, 1)).copyTo(new_image_slave.re(Range(0, 1), Range(0, 1)));

	InputMatrix.re(Range(0, 1), Range(0, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(1, nc + 1)));

	InputMatrix.re(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(0, 1), Range(nc + 2, nc + 3)));

	InputMatrix.re(Range(0, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(0, 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(0, nc)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(1, nc + 1)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 1, nr + 2), Range(nc + 2, nc + 3)));

	new_image_slave.re(Range(nr + 1, nr + 2), Range(0, nc + 3)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(0, nc + 3)));

	InputMatrix.re(Range(0, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(nc + 1, nc + 2)));

	InputMatrix.re(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.re(Range(nr + 2, nr + 3), Range(nc + 2, nc + 3)));

	new_image_slave.re(Range(0, nr + 3), Range(nc + 1, nc + 2)).copyTo(new_image_slave.re(Range(0, nr + 3), Range(nc + 2, nc + 3)));



	InputMatrix.re(Range(0, nr), Range(0, nc)).copyTo(new_image_slave.re(Range(1, nr + 1), Range(1, nc + 1)));


	//虚部

	InputMatrix.im(Range(0, 1), Range(0, 1)).copyTo(new_image_slave.im(Range(0, 1), Range(0, 1)));

	InputMatrix.im(Range(0, 1), Range(0, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(1, nc + 1)));

	InputMatrix.im(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(0, 1), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(0, 1), Range(nc + 2, nc + 3)));

	InputMatrix.im(Range(0, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, 1)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(0, 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(0, nc)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(1, nc + 1)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 1, nr + 2), Range(nc + 2, nc + 3)));

	new_image_slave.im(Range(nr + 1, nr + 2), Range(0, nc + 3)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(0, nc + 3)));

	InputMatrix.im(Range(0, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(nc + 1, nc + 2)));

	InputMatrix.im(Range(nr - 1, nr), Range(nc - 1, nc)).copyTo(new_image_slave.im(Range(nr + 2, nr + 3), Range(nc + 2, nc + 3)));

	new_image_slave.im(Range(0, nr + 3), Range(nc + 1, nc + 2)).copyTo(new_image_slave.im(Range(0, nr + 3), Range(nc + 2, nc + 3)));

	//内点
	InputMatrix.im(Range(0, nr), Range(0, nc)).copyTo(new_image_slave.im(Range(1, nr + 1), Range(1, nc + 1)));




	
	Mat image_slave_regis_re = Mat::zeros(nr, nc, CV_64F);
	Mat image_slave_regis_im = Mat::zeros(nr, nc, CV_64F);

	
	int ret;
	volatile bool parallel_flag = true;
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (int i = 0; i <= nr - 1; i++)
	{
		if (!parallel_flag) continue;
		
		for (int j = 0; j <= nc - 1; j++)
		{
			if (!parallel_flag) continue;
			double temp, offset_row, offset_col;
			Mat Row_weight = Mat::zeros(4, 1, CV_64F);
			Mat Col_weight = Mat::zeros(1, 4, CV_64F);
			offset_row = 0;
			offset_col = 0;
			ret = every_subpixel_move(i + 1, j + 1, Coefficient, &offset_row, &offset_col);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			Row_weight.at<double>(0, 0) = WeightCalculation(offset_row + 1.0);
			Row_weight.at<double>(1, 0) = WeightCalculation(offset_row);
			Row_weight.at<double>(2, 0) = WeightCalculation(1.0 - offset_row);
			Row_weight.at<double>(3, 0) = WeightCalculation(2.0 - offset_row);

			Col_weight.at<double>(0, 0) = WeightCalculation(1.0 + offset_col);
			Col_weight.at<double>(0, 1) = WeightCalculation(offset_col);
			Col_weight.at<double>(0, 2) = WeightCalculation(1.0 - offset_col);
			Col_weight.at<double>(0, 3) = WeightCalculation(2.0 - offset_col);
			temp = new_image_slave.re(Range(i, i + 4), Range(j, j + 4)).dot(Row_weight * Col_weight);

			image_slave_regis_re.at<double>(i, j) = temp;


			temp = new_image_slave.im(Range(i, i + 4), Range(j, j + 4)).dot(Row_weight * Col_weight);

			image_slave_regis_im.at<double>(i, j) = temp;
		}

	}
	if(parallel_check(parallel_flag, "interp_cubic()", parallel_error_head)) return -1;
	image_slave_regis_re.copyTo(OutputMatrix.re);
	image_slave_regis_im.copyTo(OutputMatrix.im);
	return 0;
}

double Registration::WeightCalculation(double offset)
{
	double weight = 0;
	if (offset > 0)
	{
		offset = offset;
	}
	else
	{
		offset = -offset;
	}

	if (offset < 1.0)
	{
		weight = 1.0 - 2.0 * offset * offset + offset * offset * offset;
	}
	else if (offset >= 1.0 && offset <= 2.0)
	{
		weight = 4.0 - 8.0 * offset + 5.0 * offset * offset - offset * offset * offset;
	}
	else
	{
		weight = 0.0;
	}

	return weight;
}

int Registration::registration_subpixel(ComplexMat& Master, ComplexMat& Slave, int blocksize, int interp_times)
{
	if (Master.GetRows() < 1 ||
		Master.GetCols() < 1 ||
		Master.GetRows() != Slave.GetRows() ||
		Master.GetCols() != Slave.GetCols() ||
		blocksize < 1 ||
		interp_times < 1)
	{
		fprintf(stderr, "%s\n\n", "registration_subpixel(): input check failed!\n\n");
		return -1;
	}

	int nsubr = Master.GetRows() / blocksize; //子块行数
	int nsubc = Master.GetCols() / blocksize; //子块列数
	int nsub = nsubr * nsubc;  //子块总数
	if (nsubc < 1 || nsubr < 1)
	{
		fprintf(stderr, "%s\n\n", "registration_subpixel: subblockszie, nsubc < 1 || nsubr < 1");
		return -1;
	}
	Mat sub_r_offset = Mat::zeros(nsub, 1, CV_64F); //行亚像素偏移量
	Mat sub_c_offset = Mat::zeros(nsub, 1, CV_64F); //列亚像素偏移量

	Mat m = Mat::zeros(nsub, 1, CV_32S); //子块中心行坐标
	Mat n = Mat::zeros(nsub, 1, CV_32S); //子块中心列坐标
	Mat indx = Mat::zeros(nsub, 1, CV_32S); //索引


	int count = 0;
	int ret;
	Mat coherence;
	Utils util;
	volatile bool parallel_flag = true;
#pragma omp parallel for schedule(guided) \
	private(ret)
	for (int i = 0; i < nsubc; i++)
	{
		if (!parallel_flag) continue;
		
		for (int j = 0; j < nsubr; j++)
		{
			if (!parallel_flag) continue;
			ComplexMat temp_in_slave, temp_out_slave, temp_in_master, temp_out_master;
			int offset_row, offset_col;
			//辅图像子块插值
			Slave.re(Range(j * blocksize, (j + 1) * blocksize), Range(i * blocksize, (i + 1) * blocksize)).copyTo(temp_in_slave.re);
			Slave.im(Range(j * blocksize, (j + 1) * blocksize), Range(i * blocksize, (i + 1) * blocksize)).copyTo(temp_in_slave.im);
			//util.cvmat2bin("E:\\zgb1\\InSAR\\InSAR\\bin\\re.bin", temp_in_slave.re);
			ret = interp_paddingzero(temp_in_slave, temp_out_slave, interp_times);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}

			//主图像子块插值
			Master.re(Range(j * blocksize, (j + 1) * blocksize), Range(i * blocksize, (i + 1) * blocksize)).copyTo(temp_in_master.re);
			Master.im(Range(j * blocksize, (j + 1) * blocksize), Range(i * blocksize, (i + 1) * blocksize)).copyTo(temp_in_master.im);
			ret = interp_paddingzero(temp_in_master, temp_out_master, interp_times);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}

			//实相关函数求取亚像素偏移量
			ret = real_coherent(temp_out_master, temp_out_slave, &offset_row, &offset_col);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			double offset_row_sub, offset_col_sub, mean_coh;
			//ret = util.real_coherence(temp_out_master, temp_out_slave, coherence);
			//if (ret < 0)
			//{
			//	parallel_flag = false;
			//	continue;
			//}
			//mean_coh = mean(coherence)[0];
			offset_row_sub = double(offset_row) / double(interp_times);
			offset_col_sub = double(offset_col) / double(interp_times);
			sub_r_offset.at<double>(i * nsubr + j, 0) = offset_row_sub;
			sub_c_offset.at<double>(i * nsubr + j, 0) = offset_col_sub;

			//子块中心坐标
			n.at<int>(i * nsubr + j, 0) = blocksize / 2 + i * blocksize;
			m.at<int>(i * nsubr + j, 0) = blocksize / 2 + j * blocksize;

			ret = interp_cubic(temp_in_slave, temp_in_slave, sub_r_offset.at<double>(i * nsubr + j, 0), sub_c_offset.at<double>(i * nsubr + j, 0));//子辅图像插值
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			ret = util.real_coherence(temp_in_master, temp_in_slave, coherence);
			if (ret < 0)
			{
				parallel_flag = false;
				continue;
			}
			if (mean(coherence)[0] > 0.4)
			{
				indx.at<int>(i * nsubr + j, 0) = 1;
			}

		}
	}
	if (parallel_check(parallel_flag, "registration_subpixel()", parallel_error_head)) return -1;
	int NoneZero = countNonZero(indx);//非零元素个数
	count = 0;
	if (NoneZero == 0)
	{
		fprintf(stderr, "registration_subpixel(): NoneZero == 0\n\n");
		return -1;
	}
	Mat sub_r_offset_sifted = Mat::zeros(NoneZero, 1, CV_64F); //筛选后行亚像素偏移量
	Mat sub_c_offset_sifted = Mat::zeros(NoneZero, 1, CV_64F); //筛选后列亚像素偏移量

	Mat m_sifted = Mat::zeros(NoneZero, 1, CV_32S); //筛选后子块中心行坐标
	Mat n_sifted = Mat::zeros(NoneZero, 1, CV_32S); //筛选后子块中心列坐标


	for (int k = 0; k < nsub; k++)
	{
		if (indx.at<int>(k, 0) > 0)
		{

			sub_r_offset_sifted.at<double>(count, 0) = sub_r_offset.at<double>(k, 0);
			sub_c_offset_sifted.at<double>(count, 0) = sub_c_offset.at<double>(k, 0);
			m_sifted.at<int>(count, 0) = m.at<int>(k, 0);
			n_sifted.at<int>(count, 0) = n.at<int>(k, 0);

			count++;
		}
	}


	Mat para;
	ret = all_subpixel_move(m_sifted, n_sifted, sub_r_offset_sifted, sub_c_offset_sifted, para);//拟合辅图像偏移量
	//测试
	//cout << para << "\n";
	//
	if (return_check(ret, "all_subpixel_move(*, *, *, *)", error_head)) return -1;
	ret = interp_cubic(Slave, Slave, para);
	if (return_check(ret, "interp_cubic(*, *, *)", error_head)) return -1;
	return 0;
}

int Registration::coregistration_subpixel(ComplexMat& master, ComplexMat& slave, int blocksize, int interp_times, int* offset_row,
	int* offset_col)
{
	if (master.isempty() ||
		slave.isempty() ||
		//blocksize * 5 > (slave.GetCols() < slave.GetRows() ? slave.GetCols() : slave.GetRows()) ||
		blocksize < 8||interp_times < 1 ||
		master.type() != slave.type() ||
		(master.type() != CV_64F && master.type() != CV_32F && master.type() != CV_16S)
		)
	{
		fprintf(stderr, "coregistration_pixel(): input check failed!\n");
		return -1;
	}

	//粗配准
	ComplexMat slave_r, master_small, slave_small;
	slave_r = master;
	slave_r.re = 0.0; slave_r.im = 0.0;
	int nr0 = master.GetRows() > 10000 ? 10000 : master.GetRows();
	int nc0 = master.GetCols() > 10000 ? 10000 : master.GetCols();
	nr0 = slave.GetRows() > nr0 ? nr0 : slave.GetRows();
	nc0 = slave.GetCols() > nc0 ? nc0 : slave.GetCols();

	master_small = master(cv::Range(0, nr0), cv::Range(0, nc0));
	slave_small = slave(cv::Range(0, nr0), cv::Range(0, nc0));
	master_small.convertTo(master_small, CV_64F);
	slave_small.convertTo(slave_small, CV_64F);
	int offset_rows_pre, offset_cols_pre;
	real_coherent(master_small, slave_small, &offset_rows_pre, &offset_cols_pre);

	int start_r, start_r2, end_r, end_r2, start_c, end_c, start_c2, end_c2;
	if (offset_rows_pre > 0)
	{
		start_r = offset_rows_pre;
		start_r2 = 0;
		end_r = (slave.GetRows() - offset_rows_pre) > master.GetRows() ? (offset_rows_pre + master.GetRows()) : slave.GetRows();
		end_r2 = start_r2 + (end_r - start_r);
	}
	else
	{
		//start_r = (slave.GetRows() + offset_rows_pre) > master.GetRows() ? (slave.GetRows() + offset_rows_pre - master.GetRows()) : 0;
		//start_r2 = (slave.GetRows() + offset_rows_pre) > master.GetRows() ? 0 : (master.GetRows() - slave.GetRows() - offset_rows_pre);
		//end_r = slave.GetRows() + offset_rows_pre;
		//end_r2 = master.GetRows();

		start_r = 0;
		start_r2 = -offset_rows_pre;
		end_r = (slave.GetRows() - offset_rows_pre) > master.GetRows() ? (master.GetRows() + offset_rows_pre) : slave.GetRows();
		end_r2 = (slave.GetRows() - offset_rows_pre) > master.GetRows() ? master.GetRows() : slave.GetRows() - offset_rows_pre;
	}
	if (offset_cols_pre > 0)
	{
		start_c = offset_cols_pre;
		start_c2 = 0;
		end_c = (slave.GetCols() - offset_cols_pre) > master.GetCols() ? (offset_cols_pre + master.GetCols()) : slave.GetCols();
		end_c2 = start_c2 + (end_c - start_c);
	}
	else
	{
		//start_c = (slave.GetCols() + offset_cols_pre) > master.GetCols() ? (slave.GetCols() + offset_cols_pre - master.GetCols()) : 0;
		//start_c2 = (slave.GetCols() + offset_cols_pre) > master.GetCols() ? 0 : (master.GetCols() - slave.GetCols() - offset_cols_pre);
		//end_c = slave.GetCols() + offset_cols_pre;
		//end_c2 = master.GetCols();

		start_c = 0;
		start_c2 = -offset_cols_pre;
		end_c = (slave.GetCols() - offset_cols_pre) > master.GetCols() ? (master.GetCols() + offset_cols_pre) : slave.GetCols();
		end_c2 = (slave.GetCols() - offset_cols_pre) > master.GetCols() ? master.GetCols() : slave.GetCols() - offset_cols_pre;
	}
	slave.re(cv::Range(start_r, end_r), cv::Range(start_c, end_c)).copyTo(slave_r.re(cv::Range(start_r2, end_r2), cv::Range(start_c2, end_c2)));
	slave.im(cv::Range(start_r, end_r), cv::Range(start_c, end_c)).copyTo(slave_r.im(cv::Range(start_r2, end_r2), cv::Range(start_c2, end_c2)));
	slave = slave_r;

	/*---------------------------------------*/
	/*              求取偏移量矩阵           */
	/*---------------------------------------*/
	interp_times = interp_times > 32 ? 32 : interp_times;//限定最多32倍插值
	Utils util;
	int m = (master.GetRows()) / blocksize;
	int n = (master.GetCols()) / blocksize;
	if (m * n < 10)
	{
		fprintf(stderr, "coregistration_pixel(): try smaller blocksize!\n");
		return -1;
	}
	Mat offset_r = Mat::zeros(m, n, CV_64F); Mat offset_c = Mat::zeros(m, n, CV_64F);
	Mat offset_coord_row = Mat::zeros(m, n, CV_64F); 
	Mat offset_coord_col = Mat::zeros(m, n, CV_64F);
	Mat sentinel0 = Mat::zeros(m, n, CV_64F);
	//子块中心坐标
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			offset_coord_row.at<double>(i, j) = ((double)blocksize) / 2 * (double)(2 * i + 1);
			offset_coord_col.at<double>(i, j) = ((double)blocksize) / 2 * (double)(2 * j + 1);
		}
	}
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < m; i++)
	{
		ComplexMat master_sub, slave_sub, master_sub_interp, slave_sub_interp, master1, slave1;
		Mat amplitude_slave, sign, coh1;
		int offset_row, offset_col;
		double mean_coh, coh_thresh = 0.05;
		for (int j = 0; j < n; j++)
		{
			//计算相关系数判断是否是有效数据
			master.re(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(master1.re);
			master.im(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(master1.im);
			if (master1.type() != CV_64F) master1.convertTo(master1, CV_64F);
			slave.re(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(slave1.re);
			slave.im(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(slave1.im);
			if (slave1.type() != CV_64F) slave1.convertTo(slave1, CV_64F);
			registration_pixel(master1, slave1);
			util.complex_coherence(master1, slave1, 7, 7, coh1);
			mean_coh = cv::mean(coh1)[0];
			if (mean_coh < coh_thresh)
			{
				sentinel0.at<double>(i, j) = 1.0;
				continue;
			}
			//主图像子块插值
			master.re(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(master_sub.re);
			master.im(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(master_sub.im);
			if (master_sub.type() != CV_64F) master_sub.convertTo(master_sub, CV_64F);
			interp_paddingzero(master_sub, master_sub_interp, interp_times);
			//辅图像子块插值
			slave.re(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(slave_sub.re);
			slave.im(Range(i * blocksize, (i + 1) * blocksize), Range(j * blocksize, (j + 1) * blocksize)).copyTo(slave_sub.im);
			if (slave_sub.type() != CV_64F) slave_sub.convertTo(slave_sub, CV_64F);
			amplitude_slave = slave_sub.re;
			int count_zero = 0;
			for (int ii = 0; ii < amplitude_slave.rows; ii++)
			{
				for (int jj = 0; jj < amplitude_slave.cols; jj++)
				{
					if (fabs(amplitude_slave.at<double>(ii, jj)) < 0.0000001) count_zero++;
				}
			}
			//sign = amplitude_slave < 0.0000001;
			int thresh = blocksize * blocksize / 4;
			if (count_zero > thresh)
			{
				sentinel0.at<double>(i, j) = 1.0;
				continue;
			}
			interp_paddingzero(slave_sub, slave_sub_interp, interp_times);
			//求取偏移量
			real_coherent(master_sub_interp, slave_sub_interp, &offset_row, &offset_col);
			offset_r.at<double>(i, j) = (double)offset_row / (double)interp_times;
			offset_c.at<double>(i, j) = (double)offset_col / (double)interp_times;
			
		}
	}

	/*---------------------------------------*/
	/*    拟合偏移量（将坐标做归一化处理）   */
	/*---------------------------------------*/

	/*
	* 拟合公式为 offser_row/offser_col = a0 + a1*x + a2*y
	*/
	
	////剔除outliers
	Mat sentinel = Mat::zeros(m, n, CV_64F);
	int ix, iy, count = 0, c = 0; double delta, thresh = 2.0;
	//for (int i = 0; i < m; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		count = 0;
	//		//上
	//		ix = j; 
	//		iy = i - 1; iy = iy < 0 ? 0 : iy;
	//		delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
	//		delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
	//		if (fabs(delta) >= thresh) count++;
	//		//下
	//		ix = j;
	//		iy = i + 1; iy = iy > m - 1 ? m - 1 : iy;
	//		delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
	//		delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
	//		if (fabs(delta) >= thresh) count++;
	//		//左
	//		ix = j - 1; ix = ix < 0 ? 0 : ix;
	//		iy = i; 
	//		delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
	//		delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
	//		if (fabs(delta) >= thresh) count++;
	//		//右
	//		ix = j + 1; ix = ix > n - 1 ? n - 1 : ix;
	//		iy = i;
	//		delta = fabs(offset_c.at<double>(i, j) - offset_c.at<double>(iy, ix));
	//		delta += fabs(offset_r.at<double>(i, j) - offset_r.at<double>(iy, ix));
	//		if (fabs(delta) >= thresh) count++;

	//		if (count > 2) { sentinel.at<double>(i, j) = 1.0; }
	//	}
	//}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (sentinel.at<double>(i, j) > 0.5 || sentinel0.at<double>(i, j) > 0.5) c++;
		}
	}
	Mat offset_c_0, offset_r_0, offset_coord_row_0 , offset_coord_col_0;
	offset_c_0 = Mat::zeros(m * n - c, 1, CV_64F);
	offset_r_0 = Mat::zeros(m * n - c, 1, CV_64F);
	offset_coord_row_0 = Mat::zeros(m * n - c, 1, CV_64F);
	offset_coord_col_0 = Mat::zeros(m * n - c, 1, CV_64F);
	count = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (sentinel.at<double>(i, j) < 0.5 && sentinel0.at<double>(i, j) < 0.5)
			{
				offset_r_0.at<double>(count, 0) = offset_r.at<double>(i, j);
				offset_c_0.at<double>(count, 0) = offset_c.at<double>(i, j);
				offset_coord_row_0.at<double>(count, 0) = offset_coord_row.at<double>(i, j);
				offset_coord_col_0.at<double>(count, 0) = offset_coord_col.at<double>(i, j);
				count++;
			}
		}
	}


	offset_c = offset_c_0;
	offset_r = offset_r_0;
	offset_coord_row = offset_coord_row_0;
	offset_coord_col = offset_coord_col_0;
	m = 1; n = count;
	if (count < 11)
	{
		fprintf(stderr, "coregistration_pixel(): insufficient valide sub blocks!\n");
		return -1;
	}
	double offset_x = (double)master.GetCols() / 2;
	double offset_y = (double)master.GetRows() / 2;
	double scale_x = (double)master.GetCols();
	double scale_y = (double)master.GetRows();
	offset_coord_row -= offset_y;
	offset_coord_col -= offset_x;
	offset_coord_row /= scale_y;
	offset_coord_col /= scale_x;
	Mat A = Mat::ones(m * n, 3, CV_64F);
	Mat temp, A_t;
	offset_coord_col.copyTo(A(Range(0, m * n), Range(1, 2)));

	offset_coord_row.copyTo(A(Range(0, m * n), Range(2, 3)));



	cv::transpose(A, A_t);

	Mat b_r, b_c, coef_r, coef_c, error_r, error_c, b_t, a, a_t;
	
	A.copyTo(a);
	cv::transpose(a, a_t);
	offset_r.copyTo(b_r);
	b_r = A_t * b_r;

	offset_c.copyTo(b_c);
	b_c = A_t * b_c;

	A = A_t * A;

	double rms1 = -1.0; double rms2 = -1.0;
	Mat eye = Mat::zeros(m * n, m * n, CV_64F);
	for (int i = 0; i < m * n; i++)
	{
		eye.at<double>(i, i) = 1.0;
	}
	if (cv::invert(A, error_r, cv::DECOMP_LU) > 0)
	{
		cv::transpose(offset_r, b_t);
		error_r = b_t * (eye - a * error_r * a_t) * offset_r;
		rms1 = sqrt(error_r.at<double>(0, 0) / double(m * n));
	}
	if (cv::invert(A, error_c, cv::DECOMP_LU) > 0)
	{
		cv::transpose(offset_c, b_t);
		error_c = b_t * (eye - a * error_c * a_t) * offset_c;
		rms2 = sqrt(error_c.at<double>(0, 0) / double(m * n));
	}
	if (!cv::solve(A, b_r, coef_r, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "coregistration_subpixel(): matrix defficiency!\n");
		return -1;
	}
	if (!cv::solve(A, b_c, coef_c, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "coregistration_subpixel(): matrix defficiency!\n");
		return -1;
	}

	/*---------------------------------------*/
	/*    双线性插值获取重采样后的辅图像     */
	/*---------------------------------------*/

	Mat tt(1, 3, CV_64F);
	tt.at<double>(0, 0) = 1.0;
	tt.at<double>(0, 1) = (0.0 - offset_x) / scale_x;
	tt.at<double>(0, 2) = (0.0 - offset_y) / scale_y;
	if (offset_row) *offset_row = sum(tt * coef_r)[0];
	if (offset_col) *offset_col = sum(tt * coef_c)[0];


	int rows = master.GetRows(); int cols = master.GetCols();
	ComplexMat slave_tmp;
	int master_type = master.type();
	slave_tmp = master;
	int rows_slave = slave.GetRows(); int cols_slave = slave.GetCols();
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		double x, y, ii, jj; Mat tmp(1, 3, CV_64F); Mat result;
		int mm, nn, mm1, nn1;
		double offset_rows, offset_cols, upper, lower;
		for (int j = 0; j < cols; j++)
		{
			jj = (double)j;
			ii = (double)i;
			x = (jj - offset_x) / scale_x;
			y = (ii - offset_y) / scale_y;
			tmp.at<double>(0, 0) = 1.0;
			tmp.at<double>(0, 1) = x;
			tmp.at<double>(0, 2) = y;

			result = tmp * coef_r;
			offset_rows = result.at<double>(0, 0);
			result = tmp * coef_c;
			offset_cols = result.at<double>(0, 0);

			ii += offset_rows;
			jj += offset_cols;
			
			mm = (int)floor(ii); nn = (int)floor(jj);
			if (mm < 0 || nn < 0 || mm > rows_slave - 1 || nn > cols_slave - 1)
			{
				if (master_type == CV_64F)
				{
					slave_tmp.re.at<double>(i, j) = 0.0;
					slave_tmp.im.at<double>(i, j) = 0.0;
				}
				else if (master_type == CV_32F)
				{
					slave_tmp.re.at<float>(i, j) = 0.0;
					slave_tmp.im.at<float>(i, j) = 0.0;
				}
				else
				{
					slave_tmp.re.at<short>(i, j) = 0.0;
					slave_tmp.im.at<short>(i, j) = 0.0;
				}
			}
			else
			{
				mm1 = mm + 1; nn1 = nn + 1;
				mm1 = mm1 >= rows_slave - 1 ? rows_slave - 1 : mm1;
				nn1 = nn1 >= cols_slave - 1 ? cols_slave - 1 : nn1;
				if (master_type == CV_64F)
				{
					//实部插值
					upper = slave.re.at<double>(mm, nn) + (slave.re.at<double>(mm, nn1) - slave.re.at<double>(mm, nn)) * (jj - (double)nn);
					lower = slave.re.at<double>(mm1, nn) + (slave.re.at<double>(mm1, nn1) - slave.re.at<double>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.re.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = slave.im.at<double>(mm, nn) + (slave.im.at<double>(mm, nn1) - slave.im.at<double>(mm, nn)) * (jj - (double)nn);
					lower = slave.im.at<double>(mm1, nn) + (slave.im.at<double>(mm1, nn1) - slave.im.at<double>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.im.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
				else if (master_type == CV_32F)
				{
					//实部插值
					upper = slave.re.at<float>(mm, nn) + (slave.re.at<float>(mm, nn1) - slave.re.at<float>(mm, nn)) * (jj - (double)nn);
					lower = slave.re.at<float>(mm1, nn) + (slave.re.at<float>(mm1, nn1) - slave.re.at<float>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.re.at<float>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = slave.im.at<float>(mm, nn) + (slave.im.at<float>(mm, nn1) - slave.im.at<float>(mm, nn)) * (jj - (double)nn);
					lower = slave.im.at<float>(mm1, nn) + (slave.im.at<float>(mm1, nn1) - slave.im.at<float>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.im.at<float>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
				else
				{
					//实部插值
					upper = (double)slave.re.at<short>(mm, nn) + double(slave.re.at<short>(mm, nn1) - slave.re.at<short>(mm, nn)) * (jj - (double)nn);
					lower = (double)slave.re.at<short>(mm1, nn) + double(slave.re.at<short>(mm1, nn1) - slave.re.at<short>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.re.at<short>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = (double)slave.im.at<short>(mm, nn) + double(slave.im.at<short>(mm, nn1) - slave.im.at<short>(mm, nn)) * (jj - (double)nn);
					lower = (double)slave.im.at<short>(mm1, nn) + double(slave.im.at<short>(mm1, nn1) - slave.im.at<short>(mm1, nn)) * (jj - (double)nn);
					slave_tmp.im.at<short>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
			}

		}
	}
	slave = slave_tmp;
	return 0;
}

int Registration::every_subpixel_move(int i, int j, Mat& coefficient, double* offset_row, double* offset_col)
{
	if (i < 1 || j < 1 || coefficient.cols < 1 || coefficient.rows < 12 || coefficient.type() != CV_64F)
	{
		fprintf(stderr, "every_subpixel_move(): input check failed!\n\n");
		return -1;
	}
	double a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2;
	a1 = coefficient.at<double>(0, 0);
	b1 = coefficient.at<double>(1, 0);
	c1 = coefficient.at<double>(2, 0);
	d1 = coefficient.at<double>(3, 0);
	e1 = coefficient.at<double>(4, 0);
	f1 = coefficient.at<double>(5, 0);
	a2 = coefficient.at<double>(6, 0);
	b2 = coefficient.at<double>(7, 0);
	c2 = coefficient.at<double>(8, 0);
	d2 = coefficient.at<double>(9, 0);
	e2 = coefficient.at<double>(10, 0);
	f2 = coefficient.at<double>(11, 0);

	*offset_row = a1 + b1 * i + c1 * j + d1 * i * i + e1 * j * j + f1 * i * j;
	*offset_col = a2 + b2 * i + c2 * j + d2 * i * i + e2 * j * j + f2 * i * j;
	return 0;
}

int Registration::all_subpixel_move(Mat& Coordinate_x, Mat& Coordinate_y, Mat& offset_row, Mat& offset_col, Mat& para)
{
	if (Coordinate_x.rows < 1 ||
		Coordinate_x.cols < 1 ||
		Coordinate_x.rows != Coordinate_y.rows ||
		Coordinate_x.cols != Coordinate_y.cols ||
		offset_row.rows != offset_col.rows ||
		offset_row.cols != offset_col.cols)
	{
		fprintf(stderr, "all_subpixel_move(): input size/type check failed!\n\n");
		return -1;
	}
	int N = Coordinate_x.rows;
	Coordinate_x.convertTo(Coordinate_x, CV_64F);
	Coordinate_y.convertTo(Coordinate_y, CV_64F);
	offset_row.convertTo(offset_row, CV_64F);
	offset_col.convertTo(offset_col, CV_64F);

	Mat connect_h_x[] = { Mat::ones(N, 1, CV_64F), Coordinate_x, Coordinate_y, Coordinate_x.mul(Coordinate_x), Coordinate_y.mul(Coordinate_y),
		Coordinate_y.mul(Coordinate_x) };
	Mat matrix, matrix_t;
	hconcat(connect_h_x, 6, matrix);
	Mat para1;
	transpose(matrix, matrix_t);
	if (!solve(matrix_t * matrix, matrix_t * offset_row, para1, DECOMP_LU))
	{
		fprintf(stderr, "all_subpixel_move(): cant' solve least square problem!\n");
		return -1;
	}

	Mat connect_h_y[] = { Mat::ones(N, 1, CV_64F), Coordinate_x, Coordinate_y, Coordinate_x.mul(Coordinate_x), Coordinate_y.mul(Coordinate_y),
		Coordinate_y.mul(Coordinate_x) };
	hconcat(connect_h_y, 6, matrix);
	transpose(matrix, matrix_t);
	Mat para2;
	if (!solve(matrix_t * matrix, matrix_t * offset_col, para2, DECOMP_LU))
	{
		fprintf(stderr, "all_subpixel_move(): cant' solve least square problem!\n");
		return -1;
	}
	Mat connect[] = { para1, para2 };
	vconcat(connect, 2, para);
	return 0;
}

int Registration::gcps_sift(int rows, int cols, int move_rows, int move_cols, Mat& gcps)
{
	if (fabs(move_rows) > rows ||
		fabs(move_cols) > cols ||
		rows < 1 ||
		cols < 1 ||
		gcps.cols != 5 ||
		gcps.type() != CV_64F ||
		gcps.channels() != 1||
		gcps.rows < 3)
	{
		fprintf(stderr, "gcps_sift(): input check failed!\n\n");
		return -1;
	}
	Mat gcps_tmp, GCP;
	gcps(Range(0, 2), Range(0, 5)).copyTo(GCP);
	gcps.copyTo(gcps_tmp);
	Mat index = Mat::zeros(gcps.rows, 1, CV_64F);
	int count = 0;
	//////////////////////////Scenario 1///////////////////////////////
	if (move_rows > 0 && move_cols > 0)
	{
		for (int i = 0; i < gcps_tmp.rows; i++)
		{
			if (gcps_tmp.at<double>(i, 0) <= double(rows - move_rows) &&
				gcps_tmp.at<double>(i, 1) <= double(cols - move_cols))
			{
				index.at<double>(i, 0) = 1.0;
				count++;
			}
		}
		if (count == 0)
		{
			fprintf(stderr, "all ground control points have been sifted out!\n\n");
			return -1;
		}
		gcps_tmp = Mat::zeros(count, gcps.cols, CV_64F);
		count = 0;
		for (int i = 0; i < gcps.rows; i++)
		{
			if (index.at<double>(i, 0) > 0.5)
			{
				gcps(Range(i, i + 1), Range(0, gcps.cols)).copyTo(gcps_tmp(Range(count, count + 1), Range(0, gcps.cols)));
				count++;
			}
		}
	}
	//////////////////////////Scenario 2///////////////////////////////
	if (move_rows > 0 && move_cols <= 0)
	{
		for (int i = 0; i < gcps_tmp.rows; i++)
		{
			if (gcps_tmp.at<double>(i, 0) <= double(rows - move_rows) &&
				gcps_tmp.at<double>(i, 1) >= double(1 - move_cols))
			{
				index.at<double>(i, 0) = 1.0;
				gcps_tmp.at<double>(i, 1) = gcps_tmp.at<double>(i, 1) + move_cols;
				count++;
			}
		}
		if (count == 0)
		{
			fprintf(stderr, "all ground control points have been sifted out!\n\n");
			return -1;
		}
		Mat gcps_tmp1 = Mat::zeros(count, gcps.cols, CV_64F);
		count = 0;
		for (int i = 0; i < gcps.rows; i++)
		{
			if (index.at<double>(i, 0) > 0.5)
			{
				gcps_tmp(Range(i, i + 1), Range(0, gcps.cols)).copyTo(gcps_tmp1(Range(count, count + 1), Range(0, gcps.cols)));
				count++;
			}
		}
		gcps = gcps_tmp1;
	}
	//////////////////////////Scenario 3///////////////////////////////
	if (move_rows <= 0 && move_cols > 0)
	{
		for (int i = 0; i < gcps_tmp.rows; i++)
		{
			if (gcps_tmp.at<double>(i, 0) >= double(1 - move_rows) &&
				gcps_tmp.at<double>(i, 1) <= double(cols - move_cols))
			{
				index.at<double>(i, 0) = 1.0;
				gcps_tmp.at<double>(i, 0) = gcps_tmp.at<double>(i, 0) + move_rows;
				count++;
			}
		}
		if (count == 0)
		{
			fprintf(stderr, "all ground control points have been sifted out!\n\n");
			return -1;
		}
		Mat gcps_tmp1 = Mat::zeros(count, gcps.cols, CV_64F);
		count = 0;
		for (int i = 0; i < gcps.rows; i++)
		{
			if (index.at<double>(i, 0) > 0.5)
			{
				gcps_tmp(Range(i, i + 1), Range(0, gcps.cols)).copyTo(gcps_tmp1(Range(count, count + 1), Range(0, gcps.cols)));
				count++;
			}
		}
		gcps = gcps_tmp1;
	}
	//////////////////////////Scenario 4///////////////////////////////
	if (move_rows <= 0 && move_cols <= 0)
	{
		for (int i = 0; i < gcps_tmp.rows; i++)
		{
			if (gcps_tmp.at<double>(i, 0) >= double(1 - move_rows) &&
				gcps_tmp.at<double>(i, 1) >= double(1 - move_cols))
			{
				index.at<double>(i, 0) = 1.0;
				gcps_tmp.at<double>(i, 0) = gcps_tmp.at<double>(i, 0) + move_rows;
				gcps_tmp.at<double>(i, 1) = gcps_tmp.at<double>(i, 1) + move_cols;
				count++;
			}
		}
		if (count == 0)
		{
			fprintf(stderr, "all ground control points have been sifted out!\n\n");
			return -1;
		}
		Mat gcps_tmp1 = Mat::zeros(count, gcps.cols, CV_64F);
		count = 0;
		for (int i = 0; i < gcps.rows; i++)
		{
			if (index.at<double>(i, 0) > 0.5)
			{
				gcps_tmp(Range(i, i + 1), Range(0, gcps.cols)).copyTo(gcps_tmp1(Range(count, count + 1), Range(0, gcps.cols)));
				count++;
			}
		}
		gcps = gcps_tmp1;
	}
	Mat GCPS(2 + gcps.rows, 5, CV_64F);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			GCPS.at<double>(i, j) = GCP.at<double>(i, j);
		}
	}
	for (int i = 2; i < GCPS.rows; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			GCPS.at<double>(i, j) = gcps.at<double>(i - 2, j);
		}
	}
	GCPS.copyTo(gcps);
	return 0;
}

int Registration::getDEMRgAzPos(
	Mat& DEM,
	Mat& stateVector, 
	Mat& rangePos,
	Mat& azimuthPos, 
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
	double lon_spacing,
	double lat_spacing
)
{
	if (DEM.empty() ||
		DEM.type() != CV_16S ||
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
		fprintf(stderr, "getDEMRgAzPos(): input check failed!\n");
		return -1;
	}
	//初始化轨道类
	orbitStateVectors stateVectors(stateVector, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;

	int DEM_rows = DEM.rows; int DEM_cols = DEM.cols;
	rangePos.create(DEM_rows, DEM_cols, CV_64F);
	azimuthPos.create(DEM_rows, DEM_cols, CV_64F);
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
				rangePos.at<double>(i, j) = -1.0;
				azimuthPos.at<double>(i, j) = -1.0;
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
				rangePos.at<double>(i, j) = -1.0;
				azimuthPos.at<double>(i, j) = -1.0;
			}
			else
			{
				rangePos.at<double>(i, j) = rangeIndex;
				azimuthPos.at<double>(i, j) = azimuthIndex;
			}
		}
	}
	return 0;
}

int Registration::fitSlaveOffset(Mat& slaveOffset, Mat& masterRange,
	Mat& masterAzimuth, double* a0, double* a1, double* a2)
{
	if (slaveOffset.empty() ||
		slaveOffset.type() != CV_64F || 
		masterRange.type() != CV_64F ||
		masterAzimuth.type() != CV_64F ||
		masterRange.rows != slaveOffset.rows ||
		masterRange.cols != slaveOffset.cols ||
		masterAzimuth.rows != slaveOffset.rows ||
		masterAzimuth.cols != slaveOffset.cols
		)
	{
		fprintf(stderr, "fitSlaveOffset(): input check failed!\n");
		return -1;
	}
	int count = 0, nr, nc;
	nr = slaveOffset.rows;
	nc = slaveOffset.cols;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(slaveOffset.at<double>(i, j) - invalidOffset) > 0.0001) count++;
		}
	}
	if (count < 4)
	{
		fprintf(stderr, "fitSlaveOffset(): not enough valid offsetpoints!\n");
		return -1;
	}
	Mat offset(count, 1, CV_64F);
	Mat range(count, 1, CV_64F);
	Mat azimuth(count, 1, CV_64F);
	Mat A = Mat::ones(count, 3, CV_64F);
	count = 0;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(slaveOffset.at<double>(i, j) - invalidOffset) > 0.0001)
			{
				offset.at<double>(count, 0) = slaveOffset.at<double>(i, j);
				range.at<double>(count, 0) = masterRange.at<double>(i, j);
				azimuth.at<double>(count++, 0) = masterAzimuth.at<double>(i, j);
			}
		}
	}
	range.copyTo(A(cv::Range(0, count), cv::Range(1, 2)));
	azimuth.copyTo(A(cv::Range(0, count), cv::Range(2, 3)));
	Mat A_t, b, coef;
	cv::transpose(A, A_t);
	A = A_t * A;
	b = A_t * offset;
	if (!cv::solve(A, b, coef, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "fitSlaveOffset(): matrix defficiency!\n");
		return -1;
	}
	if (a0) *a0 = coef.at<double>(0, 0);
	if (a1) *a1 = coef.at<double>(1, 0);
	if (a2) *a2 = coef.at<double>(2, 0);
	return 0;
}

int Registration::computeSlaveOffset(
	Mat& masterRange,
	Mat& masterAzimuth, 
	Mat& slaveRange,
	Mat& slaveAzimuth,
	Mat& slaveAzimuthOffset, 
	Mat& slaveRangeOffset
)
{
	if (masterRange.rows != masterAzimuth.rows ||
		masterRange.cols != masterAzimuth.cols ||
		masterRange.rows != slaveRange.rows ||
		masterRange.cols != slaveRange.cols ||
		masterRange.rows != slaveAzimuth.rows ||
		masterRange.cols != slaveAzimuth.cols
		)
	{
		fprintf(stderr, "computeSlaveOffset(): input check failed!\n");
		return -1;
	}
	slaveAzimuthOffset.create(masterRange.size(), CV_64F);
	slaveRangeOffset.create(masterRange.size(), CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < slaveAzimuthOffset.rows; i++)
	{
		for (int j = 0; j < slaveAzimuthOffset.cols; j++)
		{
			//方位向偏移量计算
			if (masterAzimuth.at<double>(i, j) < -0.5 || slaveAzimuth.at<double>(i, j) < -0.5)
			{
				slaveAzimuthOffset.at<double>(i, j) = invalidOffset;
			}
			else
			{
				slaveAzimuthOffset.at<double>(i, j) = slaveAzimuth.at<double>(i, j) - masterAzimuth.at<double>(i, j);
			}
			//距离向偏移量计算
			if (masterRange.at<double>(i, j) < -0.5 || slaveRange.at<double>(i, j) < -0.5)
			{
				slaveRangeOffset.at<double>(i, j) = invalidOffset;
			}
			else
			{
				slaveRangeOffset.at<double>(i, j) = slaveRange.at<double>(i, j) - masterRange.at<double>(i, j);
			}
		}
	}
	return 0;
}

int Registration::performBilinearResampling(
	ComplexMat& slave, 
	int dstHeight,
	int dstWidth, 
	double a0Rg, double a1Rg, double a2Rg, 
	double a0Az, double a1Az, double a2Az,
	int* offset_row,
	int* offset_col
)
{
	if (slave.isempty() || dstHeight < 2 || dstWidth < 2 ||
		(slave.type() != CV_16S && slave.type() != CV_64F && slave.type() != CV_32F)
		)
	{
		fprintf(stderr, "performBilinearResampling(): input check failed!\n");
		return -1;
	}
	ComplexMat slcResampled;
	int type = slave.type();
	if (type == CV_16S)
	{
		slcResampled.re.create(dstHeight, dstWidth, CV_16S);
		slcResampled.im.create(dstHeight, dstWidth, CV_16S);
	}
	else if (type == CV_32F)
	{
		slcResampled.re.create(dstHeight, dstWidth, CV_32F);
		slcResampled.im.create(dstHeight, dstWidth, CV_32F);
	}
	else {
		slcResampled.re.create(dstHeight, dstWidth, CV_64F);
		slcResampled.im.create(dstHeight, dstWidth, CV_64F);
	}
	
	Mat coef_r(3, 1, CV_64F), coef_c(3, 1, CV_64F);
	coef_r.at<double>(0, 0) = a0Az;
	coef_r.at<double>(1, 0) = a1Az;
	coef_r.at<double>(2, 0) = a2Az;
	coef_c.at<double>(0, 0) = a0Rg;
	coef_c.at<double>(1, 0) = a1Rg;
	coef_c.at<double>(2, 0) = a2Rg;
	if (offset_row && offset_col)
	{
		Mat tt(1, 3, CV_64F);
		tt.at<double>(0, 0) = 1.0;
		tt.at<double>(0, 1) = 0.0;
		tt.at<double>(0, 2) = 0.0;

		*offset_row = sum(tt * coef_r)[0];
		*offset_col = sum(tt * coef_c)[0];
	}
	int rows = dstHeight; int cols = dstWidth;
	int cols_slave = slave.GetCols(); int rows_slave = slave.GetRows();
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		double x, y, ii, jj; Mat tmp(1, 3, CV_64F); Mat result;
		int mm, nn, mm1, nn1;
		double offset_rows, offset_cols, upper, lower;
		for (int j = 0; j < cols; j++)
		{
			jj = (double)j;
			ii = (double)i;
			tmp.at<double>(0, 0) = 1.0;
			tmp.at<double>(0, 1) = jj;
			tmp.at<double>(0, 2) = ii;

			result = tmp * coef_r;
			offset_rows = result.at<double>(0, 0);
			result = tmp * coef_c;
			offset_cols = result.at<double>(0, 0);

			ii += offset_rows;
			jj += offset_cols;

			mm = (int)floor(ii); nn = (int)floor(jj);
			if (mm < 0 || nn < 0 || mm > rows_slave - 1 || nn > cols_slave - 1)
			{
				if (type == CV_16S)
				{
					slcResampled.re.at<short>(i, j) = 0;
					slcResampled.im.at<short>(i, j) = 0;
				}
				else if (type == CV_32F)
				{
					slcResampled.re.at<float>(i, j) = 0.0;
					slcResampled.im.at<float>(i, j) = 0.0;
				}
				else
				{
					slcResampled.re.at<double>(i, j) = 0.0;
					slcResampled.im.at<double>(i, j) = 0.0;
				}
				
			}
			else
			{
				mm1 = mm + 1; nn1 = nn + 1;
				mm1 = mm1 >= rows_slave - 1 ? rows_slave - 1 : mm1;
				nn1 = nn1 >= cols_slave - 1 ? cols_slave - 1 : nn1;
				if (type == CV_16S)
				{
					//实部插值
					upper = slave.re.at<short>(mm, nn) + (slave.re.at<short>(mm, nn1) - slave.re.at<short>(mm, nn)) * (jj - (double)nn);
					lower = slave.re.at<short>(mm1, nn) + (slave.re.at<short>(mm1, nn1) - slave.re.at<short>(mm1, nn)) * (jj - (double)nn);
					slcResampled.re.at<short>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = slave.im.at<short>(mm, nn) + (slave.im.at<short>(mm, nn1) - slave.im.at<short>(mm, nn)) * (jj - (double)nn);
					lower = slave.im.at<short>(mm1, nn) + (slave.im.at<short>(mm1, nn1) - slave.im.at<short>(mm1, nn)) * (jj - (double)nn);
					slcResampled.im.at<short>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
				else if (type == CV_32F)
				{
					//实部插值
					upper = slave.re.at<float>(mm, nn) + (slave.re.at<float>(mm, nn1) - slave.re.at<float>(mm, nn)) * (jj - (double)nn);
					lower = slave.re.at<float>(mm1, nn) + (slave.re.at<float>(mm1, nn1) - slave.re.at<float>(mm1, nn)) * (jj - (double)nn);
					slcResampled.re.at<float>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = slave.im.at<float>(mm, nn) + (slave.im.at<float>(mm, nn1) - slave.im.at<float>(mm, nn)) * (jj - (double)nn);
					lower = slave.im.at<float>(mm1, nn) + (slave.im.at<float>(mm1, nn1) - slave.im.at<float>(mm1, nn)) * (jj - (double)nn);
					slcResampled.im.at<float>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
				else
				{
					//实部插值
					upper = slave.re.at<double>(mm, nn) + (slave.re.at<double>(mm, nn1) - slave.re.at<double>(mm, nn)) * (jj - (double)nn);
					lower = slave.re.at<double>(mm1, nn) + (slave.re.at<double>(mm1, nn1) - slave.re.at<double>(mm1, nn)) * (jj - (double)nn);
					slcResampled.re.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
					//虚部插值
					upper = slave.im.at<double>(mm, nn) + (slave.im.at<double>(mm, nn1) - slave.im.at<double>(mm, nn)) * (jj - (double)nn);
					lower = slave.im.at<double>(mm1, nn) + (slave.im.at<double>(mm1, nn1) - slave.im.at<double>(mm1, nn)) * (jj - (double)nn);
					slcResampled.im.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				}
				
			}

		}
	}
	slave = slcResampled;
	return 0;
}


