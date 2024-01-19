#pragma once
#ifndef __FILTER__H__
#define __FILTER__H__
#include"..\include\Package.h"
#include"..\include\Utils.h"




class InSAR_API Filter
{
public:
	Filter();
	~Filter();
	/*2D CZT变换
	 参数1 输入矩阵
	 参数2 2D CZT变换结果（返回值）
	 参数3 CZT变换参数
	 参数4 CZT变换参数
	 参数5 CZT变换参数
	 参数6 CZT变换参数
	*/
	int czt2(Mat& src, Mat& dst, int M, int N, double theta0, double phi0);
	/*均值滤波器（原地操作）
	 参数1 输入（返回值）
	 参数2 滤波窗口大小（必须是奇数）
	*/
	int meanfilter(Mat& Src, int WndSize);
	/*2D FFTSHIFT(原地操作)*/
	int fftshift2(Mat& matrix);
	/*斜坡自适应滤波器
	 参数1 待滤波输入（干涉相位）
	 参数2 滤波结果（滤波后干涉相位）
	 参数3 滤波窗口大小（必须是奇数）
	 参数4 均值滤波窗口大小（必须是奇数）
	*/
	int slope_adaptive_filter(Mat& phase, Mat& phase_filtered, int wndsize_filter, int wndsize_prefilter);
	/** @brief 深度学习滤波（路径不要有中文）
	 
	 @param filter_dl_path              滤波程序路径
	 @param tmp_path                    中间文件保存路径
	 @param dl_model_file               深度学习模型文件
	 @param phase                       待滤波相位
	 @param phase_filtered              滤波后相位（返回值）
	*/
	int filter_dl(
		const char* filter_dl_path,
		const char* tmp_path,
		const char* dl_model_file,
		Mat& phase,
		Mat& phase_filtered
	);
	/*经典Goldstein滤波
	* 参数1 待滤波相位
	* 参数2 滤波后相位
	* 参数3 滤波器参数
	* 参数4 傅里叶变换窗口大小
	* 参数5 补零窗口大小
	*/
	int Goldstein_filter(
		Mat& phase,
		Mat& phase_filter,
		double alpha,
		int n_win,
		int n_pad
	);
	/*经典Goldstein滤波（并行）
	* 参数1 待滤波相位
	* 参数2 滤波后相位
	* 参数3 滤波器参数
	* 参数4 傅里叶变换窗口大小
	* 参数5 补零窗口大小
	*/
	int Goldstein_filter_parallel(
		Mat& phase,
		Mat& phase_filter,
		double alpha,
		int n_win,
		int n_pad
	);
	// 按二维高斯函数实现高斯滤波
	int GaussianFilter(Mat& src, Mat& dst, Mat window);
	int GenerateGaussMask(Mat& Mask, int window_height, int win_width, double sigma);
private:
	char error_head[256];
	char parallel_error_head[256];
};


#endif // !__FILTER__H__
