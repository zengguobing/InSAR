#pragma once
#ifndef __DEM__H__
#define __DEM__H__
#include"..\include\Package.h"
#include"..\include\Deflat.h"



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
	/** @brief 牛顿迭代法反演高程
	
	@param unwrapped_phase_file                            解缠相位h5文件
	@param dem                                             高程反演结果（返回值）
	@param project_path                                    工程路径
	@param iter_times                                      迭代次数
	@param mode                                            收发方式（1自发自收（默认），2单发双收）
	@return  成功返回0，否则返回-1
	*/
	int dem_newton_iter(
		const char* unwrapped_phase_file,
		Mat& dem,
		const char* project_path,
		int iter_times,
		int mode = 1
	);

	/** @brief 牛顿迭代法反演高程（测试版）
	@param unwrapped_phase_file                            解缠相位h5文件
	@param dem                                             高程反演结果（返回值）
	@param project_path                                    工程路径
	@param iter_times                                      迭代次数
	@param mode                                            收发方式（1自发自收（默认），2单发双收）
	@return  成功返回0，否则返回-1
	*/
	int dem_newton_iter_test(
		const char* unwrapped_phase_file,
		Mat& dem,
		const char* project_path,
		int iter_times,
		int mode = 1
	);

	/** @brief 牛顿迭代法反演高程（测试版）
	@param unwrapped_phase_file                            解缠相位h5文件
	@param dem                                             高程反演结果（高程返回值）
	@param lon                                             高程反演结果（经度返回值）
	@param lat                                             高程反演结果（纬度返回值）
	@param dem_x                                           高程反演结果（X坐标返回值）
	@param dem_y                                           高程反演结果（Y坐标返回值）
	@param dem_z                                           高程反演结果（Z坐标返回值）
	@param error_llh                                       控制点经纬高误差（返回值, n_gcps×3）
	@param error_xyz                                       控制点XYZ误差（返回值，n_gcps×3）
	@param project_path                                    工程路径
	@param iter_times                                      迭代次数
	@param mode                                            收发方式（1自发自收（默认），2单发双收）
	@return  成功返回0，否则返回-1
	*/
	int dem_newton_iter_14(
		const char* unwrapped_phase_file,
		Mat& dem,
		Mat& lon,
		Mat& lat,
		Mat& dem_x,
		Mat& dem_y,
		Mat& dem_z,
		Mat& error_llh,
		Mat& error_xyz,
		const char* project_path,
		int iter_times,
		int mode = 1
	);

	/** @brief 牛顿迭代法反演高程（双频乒乓模式）
	@param unwrapped_phase_file                            解缠相位h5文件
	@param dem                                             高程反演结果（高程返回值）
	@param lon                                             高程反演结果（经度返回值）
	@param lat                                             高程反演结果（纬度返回值）
	@param dem_x                                           高程反演结果（X坐标返回值）
	@param dem_y                                           高程反演结果（Y坐标返回值）
	@param dem_z                                           高程反演结果（Z坐标返回值）
	@param error_llh                                       控制点经纬高误差（返回值, n_gcps×3）
	@param error_xyz                                       控制点XYZ误差（返回值，n_gcps×3）
	@param project_path                                    工程路径
	@param iter_times                                      迭代次数
	@param mode                                            收发方式（1自发自收（默认），2单发双收）
	@return  成功返回0，否则返回-1
	*/
	int dem_newton_iter_14_dualfreqpingpong(
		const char* unwrapped_phase_file,
		Mat& dem,
		Mat& lon,
		Mat& lat,
		Mat& dem_x,
		Mat& dem_y,
		Mat& dem_z,
		Mat& error_llh,
		Mat& error_xyz,
		const char* project_path,
		int iter_times,
		int mode = 1
	);

private:
	char error_head[256];
	char parallel_error_head[256];

};


#endif // !__DEM__H__
