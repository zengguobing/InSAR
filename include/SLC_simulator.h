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
	/*@brief 根据轨道和场景DEM以及成像参数生成单视复图像（单发单收）
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
	* @param GCP                                       控制点信息（返回值，n×5，行/列/经/纬/高）
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
		ComplexMat& slc,
		Mat& GCP
	);
	/*@brief 根据轨道和场景DEM以及成像参数生成单视复图像（单发双收或者单发单收，生成2幅图）
	* @param stateVec1                                 轨道数据1（主星，信号发射轨道+信号接收轨道）
	* @param stateVec2                                 轨道数据2（辅星，信号接收轨道）
	* @param dem                                       场景DEM（来自SRTM）
	* @param lon_upperleft                             DEM左上角经度
	* @param lat_upperleft                             DEM左上角维度
	* @param sceneHeight1                              主图仿真SLC场景高度
	* @param sceneWidth1                               主图仿真SLC场景宽度
	* @param sceneHeight2                              辅图仿真SLC场景高度
	* @param sceneWidth2                               辅图仿真SLC场景宽度
	* @param nearRange1                                主图最近斜距
	* @param nearRange2                                辅图最近斜距
	* @param prf                                       脉冲重复频率
	* @param wavelength                                波长
	* @param rangeSpacing                              距离向采样间隔（m）
	* @param azimuthSpacing                            方位向采样间隔（m）
	* @param acquisitionStartTime1                     主图方位向起始时间
	* @param acquisitionStopTime1                      主图方位向结束时间
	* @param acquisitionStartTime2                     辅图方位向起始时间
	* @param acquisitionStopTime2                      辅图方位向结束时间
	* @param SNR                                       热噪声信噪比（dB）
	* @param slc1                                      主星仿真SAR图像（返回值）
	* @param slc2                                      辅星仿真SAR图像（返回值）
	* @param GCP                                       控制点信息（返回值，n×5，行/列/经/纬/高）
	* @param mode                                      收发模式（1：单发单收，2：单发双收，默认为单发单收）
	* @return 成功返回0，否则返回-1
	*/
	int generateSLC(
		Mat& stateVec1,
		Mat& stateVec2,
		Mat& dem,
		double lon_upperleft,
		double lat_upperleft,
		int sceneHeight1,
		int sceneWidth1,
		int sceneHeight2,
		int sceneWidth2,
		double nearRange1,
		double nearRange2,
		double prf,
		double wavelength,
		double rangeSpacing,
		double azimuthSpacing,
		double acquisitionStartTime1,
		double acquisitionStopTime1,
		double acquisitionStartTime2,
		double acquisitionStopTime2,
		double SNR,
		ComplexMat& slc1,
		ComplexMat& slc2,
		Mat& GCP,
		int mode = 1
	);
	/*@brief 根据轨道和场景DEM以及成像参数生成单视复图像（单发双收，乒乓模式，生成4幅图）
	* @param stateVec1                                 轨道数据1（主星，信号发射轨道+信号接收轨道）
	* @param stateVec2                                 轨道数据2（辅星，信号接收轨道）
	* @param dem                                       场景DEM（来自SRTM）
	* @param lon_upperleft                             DEM左上角经度
	* @param lat_upperleft                             DEM左上角维度
	* @param sceneHeight1                              主图仿真SLC场景高度
	* @param sceneWidth1                               主图仿真SLC场景宽度
	* @param sceneHeight2                              辅图仿真SLC场景高度
	* @param sceneWidth2                               辅图仿真SLC场景宽度
	* @param nearRange1                                主图最近斜距
	* @param nearRange2                                辅图最近斜距
	* @param prf                                       脉冲重复频率
	* @param wavelength                                波长
	* @param rangeSpacing                              距离向采样间隔（m）
	* @param azimuthSpacing                            方位向采样间隔（m）
	* @param acquisitionStartTime1                     主图方位向起始时间
	* @param acquisitionStopTime1                      主图方位向结束时间
	* @param acquisitionStartTime2                     辅图方位向起始时间
	* @param acquisitionStopTime2                      辅图方位向结束时间
	* @param SNR                                       热噪声信噪比（dB）
	* @param slc1                                      主星自发自收（返回值）
	* @param slc2                                      主星收辅星（返回值）
	* @param slc3                                      辅星自发自收（返回值）
	* @param slc4                                      辅星收主星（返回值）
	* @param GCP                                       控制点信息（返回值，n×5，行/列/经/纬/高）
	* @return 成功返回0，否则返回-1
	*/
	int generateSLC(
		Mat& stateVec1,
		Mat& stateVec2,
		Mat& dem,
		double lon_upperleft,
		double lat_upperleft,
		int sceneHeight1,
		int sceneWidth1,
		int sceneHeight2,
		int sceneWidth2,
		double nearRange1,
		double nearRange2,
		double prf,
		double wavelength,
		double rangeSpacing,
		double azimuthSpacing,
		double acquisitionStartTime1,
		double acquisitionStopTime1,
		double acquisitionStartTime2,
		double acquisitionStopTime2,
		double SNR,
		ComplexMat& slc1,
		ComplexMat& slc2,
		ComplexMat& slc3,
		ComplexMat& slc4,
		Mat& GCP
	);
private:
	char error_head[256];

};




#endif // !__SLC_SIMULATOR__H__