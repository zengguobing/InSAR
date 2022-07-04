#pragma once
#ifndef __SLC_SIMULATOR__H__
#define __SLC_SIMULATOR__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"


/*----------------------------------------*/
/*          单视复图像仿真类              */
/*----------------------------------------*/
class InSAR_API SLC_simulator
{
public:
	SLC_simulator();
	~SLC_simulator();
	/*@brief 根据入射角计算后向散射系数
	* @param incidenceAngle                             入射角
	* @param sigma                                      后向散射系数(返回值)
	* @return 成功返回0，否则返回-1
	*/
	int reflectivity(Mat& incidenceAngle, Mat& sigma);
	/*@brief 计算后向散射系数
	* @param DEM_x                                      DEM转换坐标（X）
	* @param DEM_y                                      DEM转换坐标（Y）
	* @param DEM_z                                      DEM转换坐标（Z）
	* @param satellitePos_x                             DEM对应的卫星成像点坐标（X）
	* @param satellitePos_y                             DEM对应的卫星成像点坐标（Y）
	* @param satellitePos_z                             DEM对应的卫星成像点坐标（Z）
	* @param sigma                                      后向散射系数（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int computeIncidenceAngle(
		Mat& DEM_x,
		Mat& DEM_y,
		Mat& DEM_z,
		Mat& satellitePos_x,
		Mat& satellitePos_y,
		Mat& satellitePos_z,
		Mat& sigma
	);
	/*@brief 根据轨道和场景DEM以及成像参数生成单视复图像
	* @param stateVec                                  轨道数据
	* @param dem                                       场景DEM（来自SRTM）
	* @param lon_upperleft                             DEM左上角经度
	* @param lat_upperleft                             DEM左上角维度
	* @param sceneHeight                               仿真SLC场景高度
	* @param sceneWidth                                仿真SLC场景宽度
	* @param nearRange                                 最近斜距
	* @param prf                                       脉冲重复频率
	* @param wavelength                                波长
	* @param rangeSpacing                              距离向采样间隔（m）
	* @param azimuthSpacing                            方位向采样间隔（m）
	* @param acquisitionStartTime                      方位向起始时间
	* @param acquisitionStopTime                       方位向结束时间
	* @param SNR                                       热噪声信噪比（dB）
	* @param slc                                       仿真SLC图像（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int generateSLC(
		Mat& stateVec,
		Mat& dem,
		double lon_upperleft,
		double lat_upperleft,
		int sceneHeight,
		int sceneWidth,
		double nearRange,
		double prf,
		double wavelength,
		double rangeSpacing,
		double azimuthSpacing,
		double acquisitionStartTime,
		double acquisitionStopTime,
		double SNR,
		ComplexMat& slc
	);
private:
	char error_head[256];

};




#endif // !__SLC_SIMULATOR__H__