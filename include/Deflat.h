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
	@param stateVec2         轨道数据2（辅星轨道，未插值）
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
	@param flat_phase_coef   平地相位拟合系数（返回值）
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
		Mat& flat_phase_coef
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
	/*@brief 利用SRTM高程数据和卫星轨道参数模拟地形相位
	* @param topography_phase                       地形相位（返回值）
	* @param statevector1                           主星轨道
	* @param statevector2                           辅星轨道
	* @param lon_cofficient                         坐标转换系数矩阵（图像坐标-->经度）
	* @param lat_cofficient                         坐标转换系数矩阵（图像坐标-->纬度）
	* @param inc_cofficient                         坐标转换系数矩阵（图像坐标-->下视角）
	* @param prf1                                   主星prf
	* @param prf2                                   辅星prf
	* @param sceneHeight                            主图行数
	* @param sceneWidth                             主图列数
	* @param offset_row                             主图在原图像中的行偏移
	* @param offset_col                             主图在原图像中的列偏移
	* @param nearRangeTime                          主图最近斜距时间（s）
	* @param rangeSpacing                           主图斜距向采样间隔（m）
	* @param wavelength                             波长
	* @param acquisition_start_time                 主图原图像方位向采样开始时间
	* @param acquisition_stop_time                  主图原图像方位向采样结束时间
	* @param DEMpath                                下载的SRTM高程文件保存路径
	* @param interp_times                           SRTM高程插值倍数（默认为20）
	* @return 成功返回0，否则返回-1
	*/
	int topography_simulation(
		Mat& topography_phase,
		Mat& statevector1,
		Mat& statevector2,
		Mat& lon_cofficient,
		Mat& lat_cofficient,
		Mat& inc_cofficient,
		double prf1, double prf2,
		int sceneHeight, int sceneWidth,
		int offset_row, int offset_col,
		double nearRangeTime, double rangeSpacing, double wavelength,
		double acquisition_start_time, double acquisition_stop_time,
		const char* DEMpath,
		int interp_times = 20
	);
	/*@brief 将WGS84坐标DEM投影到相应的SAR坐标系中
	* @param DEM84                        84坐标系DEM（short型矩阵）
	* @param mappedDEM                    投影DEM（返回值,short型矩阵）
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
	* @param stateVector                  卫星轨道数据（未插值）
	* @param interp_times                 84坐标系DEM插值倍数（默认值为10）
	* @param lon_spacing                  84坐标系DEM经度采样间隔（°）
	* @param lat_spacing                  84坐标系DEM纬度采样间隔（°）
	* @return 成功返回0，否则返回-1
	*/
	int demMapping(
		Mat& DEM84,
		Mat& mappedDEM,
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
		int interp_times = 10,
		double lon_spacing = 5.0 / 6000.0,
		double lat_spacing = 5.0 / 6000.0
	);
	/*@brief 将WGS84坐标DEM投影到相应的SAR坐标系中（投影经纬度也返回）
	* @param DEM84                        84坐标系DEM（short型矩阵）
	* @param mappedDEM                    投影DEM（返回值,short型矩阵）
	* @param mappedLat                    投影纬度坐标返回值，double型矩阵
	* @param mappedLon                    投影经度坐标（返回值，double型矩阵）
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
	* @param stateVector                  卫星轨道数据（未插值）
	* @param interp_times                 84坐标系DEM插值倍数（默认值为10）
	* @param lon_spacing                  84坐标系DEM经度采样间隔（°）
	* @param lat_spacing                  84坐标系DEM纬度采样间隔（°）
	* @return 成功返回0，否则返回-1
	*/
	int demMapping(
		Mat& DEM84,
		Mat& mappedDEM,
		Mat& mappedLat,
		Mat& mappedLon,
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
		int interp_times = 10,
		double lon_spacing = 5.0 / 6000.0,
		double lat_spacing = 5.0 / 6000.0
	);
	/*@brief 去除配准SLC图像中的参考地形相位（包括平地相位）
	* @param slc_deramped                           去参考相位后SLC图像（返回值）
	* @param mappedDEM                              配准主图像坐标系DEM
	* @param mappedLat                              DEM纬度坐标
	* @param mappedLon                              DEM经度坐标
	* @param slcH5File                              配准SLC图像h5文件
	* @return 成功返回0，否则返回-1
	*/
	int SLC_deramp(
		ComplexMat& slc,
		Mat& mappedDEM,
		Mat& mappedLat,
		Mat& mappedLon,
		const char* slcH5File
	);
	/*@brief 去除配准SAR图像数据堆栈的参考地形相位（包括平地相位）
	* @param SLCH5Files                             配准SAR图像数据堆栈h5文件
	* @param reference                              参考主图像序号
	* @param demPath                                DEM下载保存路径
	* @param outSLCH5Files                          去参考相位后SAR图像数据堆栈h5文件
	* @return 成功返回0，否则返回-1
	*/
	int SLCs_deramp(
		vector<string>& SLCH5Files,
		int reference,
		const char* demPath,
		vector<string>& outSLCH5Files
	);
	/*@brief 根据投影至SAR坐标系的DEM和卫星系统参数模拟地形相位
	* @param mappedDEM                              投影至SAR坐标系的DEM（short型）
	* @param topography_phase                       地形相位（返回值）
	* @param inc_coefficient                        下视角拟合系数
	* @param B_effect                               垂直基线长度
	* @param nearRangeTime                          最近斜距时间
	* @param offset_row                             SAR图像在原SAR图像中的行偏移量
	* @param offset_col                             SAR图像在原SAR图像中的列偏移量
	* @param wavelength                             波长
	* @param rangeSpacing                           距离向采样间隔
	* @return 成功返回0，否则返回-1
	*/
	int topography_phase_simulation(
		Mat& mappedDEM,
		Mat& topography_phase,
		Mat& inc_coefficient,
		double B_effect,
		double nearRangeTime,
		int offset_row,
		int offset_col,
		double wavelength,
		double rangeSpacing
	);
	/*@brief 计算所需DEM的地理边界
	* @param lat_coefficient                        地理坐标转换系数（行列-->纬度）
	* @param lon_coefficient                        地理坐标转换系数（行列-->经度）
	* @param sceneHeight                            场景高度
	* @param sceneWidth                             场景宽度
	* @param offset_row                             场景在原图像中的行偏移量
	* @param offset_col                             场景在原图像中的列偏移量
	* @param lonMax                                 最大经度（返回值）
	* @param latMax                                 最大纬度（返回值）
	* @param lonMin                                 最小经度（返回值）
	* @param latMin                                 最小纬度（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int computeImageGeoBoundry(
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
	);
	/*@brief 根据地理边界信息计算所需下载的SRTM高程文件名
	* @param lonMin                       最小经度
	* @param lonMax                       最大经度
	* @param latMin                       最小纬度
	* @param latMax                       最大纬度
	* @param name                         文件名
	* @return 成功返回0，否则返回-1
	*/
	static int getSRTMFileName(
		double lonMin,
		double lonMax,
		double latMin,
		double latMax,
		vector<string>& name
	);
	/*@brief 根据地理边界信息获取SRTM高程
	* @param filepath                     下载的SRTM高程文件保存路径
	* @param DEM_out                      DEM数据（返回值，short型）
	* @param lonUpperLeft                 左上角经度（返回值）
	* @param latUpperLeft                 左上角纬度（返回值）
	* @param lonMin                       最小经度
	* @param lonMax                       最大经度
	* @param latMin                       最小纬度
	* @param latMax                       最大纬度
	* @return 成功返回0，否则返回-1
	*/
	int getSRTMDEM(
		const char* filepath,
		Mat& DEM_out,
		double* lonUpperLeft,
		double* latUpperLeft,
		double lonMin,
		double lonMax,
		double latMin,
		double latMax
	);
	/*@brief 下载SRTM高程数据
	* @param name                         文件名
	* @return 成功返回0，否则返回-1
	*/
	int downloadSRTM(const char* name);

private:
	char error_head[256];
	char parallel_error_head[256];
	/*DEM数据（高程为short型数据）*/
	Mat rawDEM;
	/*DEM数据行数*/
	int rows;
	/*DEM数据列数*/
	int cols;
	///*左上角经度*/
	//double lonUpperLeft;
	///*左上角纬度*/
	//double latUpperLeft;
	///*右下角经度*/
	//double lonLowerRight;
	///*右下角纬度*/
	//double latLowerRight;
	/*DEM文件保存路径*/
	string DEMPath;
	/*纬度采样间隔（默认为5.0/6000.0）*/
	double latSpacing;
	/*经度采样间隔（默认为5.0/6000.0）*/
	double lonSpacing;
	/*SRTM全球高程url*/
	string SRTMURL = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/";
};


#endif // !__DEFLAT__H__
