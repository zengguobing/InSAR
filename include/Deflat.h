#pragma once
#ifndef __DEFLAT__H__
#define __DEFLAT__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"
#include"..\include\Utils.h"



class InSAR_API Deflat
{
public:
	Deflat();
	~Deflat();
	/*卫星轨道拟合
	 参数1 轨道数据
	 参数2 拟合系数（返回值）
	*/
	int Orbit_Polyfit(Mat& Orbit, Mat& coef);
	/*获取地面点的卫星成像方位向时刻
	 参数1 卫星轨道中心方位向时刻
	 参数2 轨道拟合系数
	 参数3 地面点坐标
	 参数4 卫星轨道起始方位向时刻
	 参数5 地面点的卫星成像方位向时刻（返回值）
	*/
	int get_satellite_aztime_NEWTON(double center_time, Mat& coef, Mat pos_xyz, double start_time, double* aztime);
	/*获取给定方位向时刻的卫星位置
	 参数1 给定方位向时刻
	 参数2 轨道拟合系数
	 参数3 卫星位置（返回值）
	*/
	int get_xyz(double aztime, Mat& coef, Mat& pos_xyz);
	/*获取给定方位向时刻的卫星速度
	参数1 给定方位向时刻
	参数2 轨道拟合系数
	参数3 卫星速度（返回值）
	*/
	int get_vel(double aztime, Mat& coef, Mat& vel_xyz);
	/*获取给定方位向时刻的卫星加速度
	参数1 给定方位向时刻
	参数2 轨道拟合系数
	参数3 卫星加速度（返回值）
	*/
	int get_acc(double aztime, Mat& coef, Mat& acc_xyz);
	/*去平地函数
	 参数1 待去平地干涉相位
	 参数2 去平地干涉相位
	 参数3 平地干涉相位
	 参数4 辅助参数（雷达载频等参数）
	 参数5 去平地控制点
	 参数6 主星轨道
	 参数7 辅星轨道
	 参数8 收发模式（单发单收1，单发双收2）
	 参数9 干涉多视倍数（>= 1）
	*/
	int deflat(	
		Mat& phase,
		Mat& phase_deflat,
		Mat& flat_phase,
		Mat auxi,
		Mat gcps,
		Mat orbit_main,
		Mat orbit_slave,
		int mode,
		int multilook_times
	);
	/** @brief 去平地
	
	@param stateVec1         轨道数据1（主星轨道，未插值）
	@param stateVec2         轨道数据2（主星轨道，未插值）
	@param lon_coef          坐标转换系数矩阵（图像坐标-->经度）
	@param lat_coef          坐标转换系数矩阵（图像坐标-->纬度）
	@param phase             输入相位
	@param offset_row        主图像左上角在原始图像中的行偏移量
	@param offset_col        主图像左上角在原始图像中的列偏移量
	@param height            平地高度
	@param time_interval1    轨道1插值时间间隔（1/prf）
	@param time_interval2    轨道2插值时间间隔（1/prf）
	@param mode              收发方式(1：自发自收，2：单发双收)
	@param wavelength        波长
	@param phase_deflated    去平地相位（返回值）
	@param flat_phase        平地相位（返回值，未缠绕）
	*/
	int deflat(
		const Mat& stateVec1,
		const Mat& stateVec2,
		const Mat& lon_coef,
		const Mat& lat_coef,
		const Mat& phase,
		int offset_row,
		int offset_col,
		double height,
		double time_interval1,
		double time_interval2,
		int mode,
		double wavelength,
		Mat& phase_deflated,
		Mat& flat_phase
	);
	/** @brief 利用外部DEM数据去除地形相位

	@param phase                 输入相位
	@param dem                   外部DEM数据
	@param dem_range_lon         外部DEM经度范围（1×2矩阵，第一个最小经度，第二个最大经度）
	@param dem_range_lat         外部DEM纬度范围（1×2矩阵，第一个最小经度，第二个最大经度）
	@param stateVec1             主星轨道数据（未插值）
	@param stateVec2             辅星轨道数据（未插值）
	@param lon_coef              坐标转换系数矩阵（图像坐标-->经度）
	@param lat_coef              坐标转换系数矩阵（图像坐标-->纬度）
	@param inc_coef              下视角拟合系数（1×11矩阵）
	@param offset_row            主图像左上角在原始图像中的行偏移量
	@param offset_col            主图像左上角在原始图像中的列偏移量
	@param interp_interval1      主星轨道插值时间间隔（1/prf）
	@param interp_interval2      辅星轨道插值时间间隔（1/prf）
	@param mode                  收发方式(1：自发自收，2：单发双收)
	@param wavelength            波长
	@param B_effect              垂直基线
	@param phase_detopo          去地形相位（输出值）
	*/
	int topo_removal(
		const Mat& phase,
		const Mat& dem,
		const Mat& dem_range_lon,
		const Mat& dem_range_lat,
		const Mat& stateVec1,
		const Mat& stateVec2,
		const Mat& lon_coef,
		const Mat& lat_coef,
		const Mat& inc_coef,
		int offset_row,
		int offset_col,
		double interp_interval1,
		double interp_interval2,
		int mode,
		double wavelength,
		double B_effect,
		Mat& phase_detopo
	);
private:
	char error_head[256];
	char parallel_error_head[256];

};


#endif // !__DEFLAT__H__
