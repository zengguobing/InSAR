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
	/** @brief 精配准（支持16位整型和64位浮点型输入）
	
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

	/*@brief 根据DEM和卫星轨道等辅助参数计算DEM点在SAR图像中的位置
	* @param DEM                          DEM（short型矩阵）
	* @param stateVector                  卫星轨道数据（未插值）
	* @param rangePos                     DEM点在SAR图像中距离向坐标（列，double型矩阵，返回值）
	* @param azimuthPos                   DEM点在SAR图像中方位向坐标（行，double型矩阵，返回值）
	* @param lon_upperleft                84坐标系DEM左上角经度
	* @param lat_upperleft                84坐标系DEM左上角纬度
	* @param offset_row                   SAR图像在原场景中的行偏移量
	* @param offset_col                   SAR图像在原场景中的列偏移量
	* @param sceneHeight                  SAR图像场景高度
	* @param sceneWidth                   SAR图像场景宽度
	* @param prf                          SAR卫星雷达脉冲重复频率
	* @param rangeSpacing                 距离向采样间隔（m）
	* @param wavelength                   波长
	* @param nearRangeTime                最近斜距时间
	* @param acquisitionStartTime         方位向采样开始时间
	* @param acquisitionStopTime          方位向采样结束时间
	* @param lon_spacing                  84坐标系DEM经度采样间隔（°）
	* @param lat_spacing                  84坐标系DEM纬度采样间隔（°）
	* @return 成功返回0，否则返回-1
	*/
	int getDEMRgAzPos(
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
	);
	/*@brief 拟合辅图像偏移（1阶拟合，offset = a0 + a1 * x + a2 * y）
	* @param slaveOffset                           偏移量
	* @param masterRange                           DEM点在主图中的距离向坐标（列数，double型矩阵）
	* @param masterAzimuth                         DEM点在主图中的方位向坐标（行数，double型矩阵）
	* @param a0                                    拟合系数
	* @param a1                                    拟合系数
	* @param a2                                    拟合系数
	* @return 成功返回0，否则返回-1
	*/
	int fitSlaveOffset(
		Mat& slaveOffset,
		Mat& masterRange,
		Mat& masterAzimuth,
		double* a0,
		double* a1,
		double* a2
	);
	/*@brief 计算辅图像偏移
	* @param masterRange                           DEM点在主图中的距离向坐标（列数）
	* @param masterAzimuth                         DEM点在主图中的方位向坐标（行数）
	* @param slaveRange                            DEM点在辅图中的距离向坐标（列数）
	* @param slaveAzimuth                          DEM点在辅图中的方位向坐标（行数）
	* @param slaveAzimuthOffset                    辅图像方位向偏移（返回值）
	* @param slaveRangeOffset                      辅图像距离向偏移（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int computeSlaveOffset(
		Mat& masterRange,
		Mat& masterAzimuth,
		Mat& slaveRange,
		Mat& slaveAzimuth,
		Mat& slaveAzimuthOffset,
		Mat& slaveRangeOffset
	);
	/*@brief 复图像双线性插值重采样（inplace，原地操作）
	* @param slc                                   待重采样图像（原地操作）
	* @param dstHeight                             重采样图像高度
	* @param dstWidth                              重采样图像宽度
	* @param a0Rg                                  距离向偏移拟合系数
	* @param a1Rg                                  距离向偏移拟合系数
	* @param a2Rg                                  距离向偏移拟合系数
	* @param a0Az                                  方位向偏移拟合系数
	* @param a1Az                                  方位向偏移拟合系数
	* @param a2Az                                  方位向偏移拟合系数
	* @param offset_row                            辅图像左上角相对于主图像的行偏移量（返回值）
	* @param offset_col                            辅图像左上角相对于主图像的列偏移量（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int performBilinearResampling(
		ComplexMat& slc,
		int dstHeight,
		int dstWidth,
		double a0Rg, double a1Rg, double a2Rg,
		double a0Az, double a1Az, double a2Az,
		int* offset_row = NULL,
		int* offset_col = NULL
	);
private:
	char error_head[256];
	char parallel_error_head[256];
	double invalidOffset = -9999.0;

};




#endif // !__REGISTRATION__H__
