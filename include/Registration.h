#pragma once
#ifndef __REGISTRATION__H__
#define __REGISTRATION__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"
#include"..\include\Utils.h"



class InSAR_API Registration
{
public:
	Registration();
	~Registration();
	/*求取两幅辅图像的实相关函数
	 参数1 主图像（复）
	 参数2 辅图像（复）
	 参数3 行偏移量（返回值）
	 参数4 列偏移量（返回值）
	*/
	int real_coherent(ComplexMat& Master, ComplexMat& Slave, int* offset_row, int* offset_col);
	/*2D FFTSHIFT(原地操作)*/
	int fftshift2(Mat& matrix);
	/*2D FFT
	 参数1 输入矩阵
	 参数2 输出矩阵
	*/
	int fft2(Mat& Src, Mat& Dst);
	/*像元级配准（原地操作）
	 参数1 主图像（复）（输入值/返回值）
	 参数2 辅图像（复）（输入值/返回值）
	 参数3 行偏移量（返回值）
	 参数4 列偏移量（返回值）
	*/
	int registration_pixel(ComplexMat& Master, ComplexMat& Slave, int* move_r = NULL, int* move_c = NULL);
	/*频域补零插值
	 参数1 输入矩阵（复）
	 参数2 输出矩阵（复）
	 参数3 插值倍数（大于1， 且是2的n次幂）
	*/
	int interp_paddingzero(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, int interp_times);
	/*立方插值
	 参数1 输入矩阵
	 参数2 输出矩阵
	 参数3 行偏移量
	 参数4 列偏移量
	*/
	int interp_cubic(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, double offset_row, double offset_col);
	/*立方插值
	参数1 输入矩阵
	参数2 输出矩阵
	参数3 拟合系数
	*/
	int interp_cubic(ComplexMat& InputMatrix, ComplexMat& OutputMatrix, Mat& Coefficient);
	/*计算每个像素的偏移量
	 参数1 像素行号
	 参数2 像素列号
	 参数3 拟合系数
	 参数4 行偏移量
	 参数5 列偏移量
	*/
	int every_subpixel_move(int i, int j, Mat& coefficient, double* offset_row, double* offset_col);
	/*计算卷积核权重*/
	double WeightCalculation(double offset);
	/*亚像素级配准
	 参数1 主图像（复）
	 参数2 辅图像（复）
	 参数3 子块大小（大于1）
	 参数4 插值倍数（大于1）
	*/
	int registration_subpixel(ComplexMat& Master, ComplexMat& Slave, int blocksize, int interp_times);
	/** @brief 精配准
	
	@param master              主图像
	@param slave               辅图像
	@param blocksize           主图像参考块分块大小（blocksize×blocksize，blocksize为2的n次幂）
	@param interp_times        插值倍数(InSAR要求至少8倍插值)
	*/
	int coregistration_subpixel(
		ComplexMat& master,
		ComplexMat& slave,
		int blocksize,
		int interp_times
	);
	/*拟合像素偏移量
	 参数1 行序列号
	 参数2 列序列号
	 参数3 行偏移量
	 参数4 列偏移量
	 参数5 拟合系数（返回值）
	*/
	int all_subpixel_move(Mat& Coordinate_x, Mat& Coordinate_y, Mat& offset_row, Mat& offset_col, Mat& para);
	/*根据粗配准偏移量筛选控制点
	* 参数1 原始图像行数
	* 参数2 原始图像列数
	* 参数3 行偏移
	* 参数4 列偏移
	* 参数5 控制点信息
	*/
	int gcps_sift(int rows, int cols, int move_rows, int move_cols, Mat& gcps);


private:
	char error_head[256];
	char parallel_error_head[256];

};




#endif // !__REGISTRATION__H__
