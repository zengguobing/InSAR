#pragma once
#ifndef __DEM__H__
#define __DEM__H__
#include"..\include\Package.h"
#include"..\include\Deflat.h"
#ifdef _DEBUG
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "Deflat_d.lib")
#else
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "Deflat.lib")
#endif // _DEBUG


class InSAR_API Dem
{
public:
	Dem();
	~Dem();
	/*相位高程转换
	* 参数1 解缠相位
	* 参数2 平地相位
	* 参数3 高程值（返回值）
	* 参数4 主星辅助参数
	* 参数5 辅星辅助参数
	* 参数6 主星轨道参数
	* 参数7 辅星轨道参数
	* 参数8 多普勒中心频率
	* 参数9 地面控制点
	* 参数10 配准结果（包括配准偏移量和配准前图像尺寸）
	* 参数11 主星成像间隔时间
	* 参数12 辅星成像间隔时间
	* 参数13 干涉多视倍数
	* 参数14 收发模式（1：单发单收，2：单发双收）
	* 参数15 牛顿迭代次数（默认30次）
	*/
	int phase2dem_newton_iter(
		Mat unwrapped_phase,
		Mat flat_phase,
		Mat& dem,
		Mat auxi_m,
		Mat auxi_s,
		Mat orbit_m,
		Mat orbit_s,
		Mat doppler_frequency,
		Mat gcps,
		Mat regis_out,
		double delta_m,
		double delta_s,
		int multilook_times,
		int mode, 
		int iters
	);

private:
	char error_head[256];
	char parallel_error_head[256];

};

Dem::Dem()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "DEM_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "DEM_DLL_ERROR: error happens when using parallel computing in function: ");
}

Dem::~Dem()
{
}

#endif // !__DEM__H__
