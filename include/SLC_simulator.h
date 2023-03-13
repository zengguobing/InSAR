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
	/*@brief 根据轨道和场景DEM以及成像参数参考斜距（单发双收，乒乓模式）
	* @param stateVec1                                 轨道数据1（主星，信号发射轨道+信号接收轨道）
	* @param stateVec2                                 轨道数据2（辅星，信号接收轨道）
	* @param dem                                       场景DEM（来自SRTM）
	* @param lon_upperleft                             DEM左上角经度
	* @param lat_upperleft                             DEM左上角维度
	* @param sceneHeight1                              主图仿真SLC场景高度
	* @param sceneWidth1                               主图仿真SLC场景宽度
	* @param nearRange1                                主图最近斜距
	* @param prf                                       脉冲重复频率
	* @param wavelength                                波长
	* @param rangeSpacing                              距离向采样间隔（m）
	* @param azimuthSpacing                            方位向采样间隔（m）
	* @param acquisitionStartTime1                     主图方位向起始时间
	* @param acquisitionStopTime1                      主图方位向结束时间
	* @param acquisitionStartTime2                     辅图方位向起始时间
	* @param acquisitionStopTime2                      辅图方位向结束时间
	* @param R1                                        主星斜距
	* @param R2                                        辅星斜距
	* @return 成功返回0，否则返回-1
	*/
	int generateSlantrange(
		Mat& stateVec1,
		Mat& stateVec2,
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
		double acquisitionStartTime1,
		double acquisitionStopTime1,
		double acquisitionStartTime2,
		double acquisitionStopTime2,
		Mat& R1,
		Mat& R2
	);
	/*@brief 乒乓模式去参考相位
	* @param mappedDEM                              配准主图像坐标系DEM
	* @param mappedLat                              DEM纬度坐标
	* @param mappedLon                              DEM经度坐标
	* @param slcH5File1                             配准SLC图像h5文件（主星发主星收）
	* @param slcH5file2                             配准SLC图像h5文件（辅星发主星收）
	* @param slcH5file3                             配准SLC图像h5文件（辅星发辅星收）
	* @param slcH5file4                             配准SLC图像h5文件（主星发辅星收）
	* @param slcH5File1_out                         配准SLC图像h5文件（主星发主星收，去参考相位返回值）
	* @param slcH5file2_out                         配准SLC图像h5文件（辅星发主星收，去参考相位返回值）
	* @param slcH5file3_out                         配准SLC图像h5文件（辅星发辅星收，去参考相位返回值）
	* @param slcH5file4_out                         配准SLC图像h5文件（主星发辅星收，去参考相位返回值）
	* @return 成功返回0，否则返回-1
	*/
	int SLC_deramp(
		Mat& mappedDEM,
		Mat& mappedLat,
		Mat& mappedLon,
		const char* slcH5File1,
		const char* slcH5File2,
		const char* slcH5File3,
		const char* slcH5File4,
		const char* slcH5File1_out,
		const char* slcH5File2_out,
		const char* slcH5File3_out,
		const char* slcH5File4_out
	);
	/*@brief 去参考相位
	* @param master_index                           主图序列号（1-based）
	* @param mappedDEM                              配准主图像坐标系DEM
	* @param mappedLat                              DEM纬度坐标
	* @param mappedLon                              DEM经度坐标
	* @param mode                                   收发模式（1：单发单收，2：单发双收，3：乒乓模式，4：双频乒乓模式）
	* @param slcH5FilesList                         配准h5文件数组
	* @param slcH5FileListOut                       去参考后h5文件数组
	* @return 成功返回0，否则返回-1
	*/
	int SLC_deramp_14(
		vector<string>& slcH5FilesList,
		vector<string>& slcH5FilesListOut,
		int master_index,
		Mat& mappedDEM,
		Mat& mappedLat,
		Mat& mappedLon,
		int mode
	);
	/*@brief 乒乓模式重新加入参考相位
	* @param mappedDEM                              配准主图像坐标系DEM
	* @param mappedLat                              DEM纬度坐标
	* @param mappedLon                              DEM经度坐标
	* @param slcH5File1                             配准SLC图像h5文件（主星发主星收）
	* @param slcH5file2                             配准SLC图像h5文件（辅星发主星收）
	* @param slcH5file3                             配准SLC图像h5文件（辅星发辅星收）
	* @param slcH5file4                             配准SLC图像h5文件（主星发辅星收）
	* @param slcH5File1_out                         配准SLC图像h5文件（主星发主星收，加参考相位返回值）
	* @param slcH5file2_out                         配准SLC图像h5文件（辅星发主星收，加参考相位返回值）
	* @param slcH5file3_out                         配准SLC图像h5文件（辅星发辅星收，加参考相位返回值）
	* @param slcH5file4_out                         配准SLC图像h5文件（主星发辅星收，加参考相位返回值）
	* @return 成功返回0，否则返回-1
	*/
	int SLC_reramp(
		Mat& mappedDEM,
		Mat& mappedLat,
		Mat& mappedLon,
		const char* slcH5File1,
		const char* slcH5File2,
		const char* slcH5File3,
		const char* slcH5File4,
		const char* slcH5File1_out,
		const char* slcH5File2_out,
		const char* slcH5File3_out,
		const char* slcH5File4_out
	);
	/*@brief 乒乓模式干涉相位估计
	* @param estimation_wndsize                     估计窗口大小
	* @param slcH5File1                             配准SLC图像h5文件（主星发主星收）
	* @param slcH5file2                             配准SLC图像h5文件（辅星发主星收）
	* @param slcH5file3                             配准SLC图像h5文件（辅星发辅星收）
	* @param slcH5file4                             配准SLC图像h5文件（主星发辅星收）
	* @param slcH5File1_out                         配准SLC图像h5文件（主星发主星收，返回值）
	* @param slcH5file2_out                         配准SLC图像h5文件（辅星发主星收，返回值）
	* @param slcH5file3_out                         配准SLC图像h5文件（辅星发辅星收，返回值）
	* @param slcH5file4_out                         配准SLC图像h5文件（主星发辅星收，返回值）
	* @return 成功返回0，否则返回-1
	*/
	int MB_phase_estimation(
		int estimation_wndsize,
		const char* slcH5File1,
		const char* slcH5File2,
		const char* slcH5File3,
		const char* slcH5File4,
		const char* slcH5File1_out,
		const char* slcH5File2_out,
		const char* slcH5File3_out,
		const char* slcH5File4_out
	);
	/*@brief 双频乒乓模式干涉相位估计
	* @param estimation_wndsize                     估计窗口大小
	* @param slcH5File1                             配准SLC图像h5文件（频率1，主星发主星收）
	* @param slcH5file2                             配准SLC图像h5文件（频率1，辅星发主星收）
	* @param slcH5file3                             配准SLC图像h5文件（频率1，辅星发辅星收）
	* @param slcH5file4                             配准SLC图像h5文件（频率1，主星发辅星收）
	* @param slcH5File5                             配准SLC图像h5文件（频率2，主星发主星收）
	* @param slcH5file6                             配准SLC图像h5文件（频率2，辅星发主星收）
	* @param slcH5file7                             配准SLC图像h5文件（频率2，辅星发辅星收）
	* @param slcH5file8                             配准SLC图像h5文件（频率2，主星发辅星收）
	* @param slcH5File1_out                         配准SLC图像h5文件（频率1，主星发主星收，返回值）
	* @param slcH5file2_out                         配准SLC图像h5文件（频率1，辅星发主星收，返回值）
	* @param slcH5file3_out                         配准SLC图像h5文件（频率1，辅星发辅星收，返回值）
	* @param slcH5file4_out                         配准SLC图像h5文件（频率1，主星发辅星收，返回值）
	* @param slcH5File5_out                         配准SLC图像h5文件（频率2，主星发主星收，返回值）
	* @param slcH5file6_out                         配准SLC图像h5文件（频率2，辅星发主星收，返回值）
	* @param slcH5file7_out                         配准SLC图像h5文件（频率2，辅星发辅星收，返回值）
	* @param slcH5file8_out                         配准SLC图像h5文件（频率2，主星发辅星收，返回值）
	* @return 成功返回0，否则返回-1
	*/
	int MB_phase_estimation(
		int estimation_wndsize,
		const char* slcH5File1,
		const char* slcH5File2,
		const char* slcH5File3,
		const char* slcH5File4,
		const char* slcH5File5,
		const char* slcH5File6,
		const char* slcH5File7,
		const char* slcH5File8,
		const char* slcH5File1_out,
		const char* slcH5File2_out,
		const char* slcH5File3_out,
		const char* slcH5File4_out,
		const char* slcH5File5_out,
		const char* slcH5File6_out,
		const char* slcH5File7_out,
		const char* slcH5File8_out
	);
	/*@brief 乒乓模式双频无模糊最大似然相位估计
	* @param phase_reference                参考相位（用于最大似然搜索的区间确定）
	* @param wrapped_phase_low              低频相位
	* @param wrapped_phase_high             高频相位
	* @param outphase                       融合相位
	* @param lat_coef                       纬度转换系数
	* @param lon_coef                       经度转换系数
	* @param wavelength_low                 低频波长
	* @param wavelength_high                高频波长
	* @param nearRangeTime                  主图最近斜距时间（s）
	* @param rangeSpacing                   主图距离向采样间隔（m）
	* @param azimuthSpacing                 主图方位向采样间隔（m）
	* @param offset_row                     主图像在原图中的行偏移
	* @param offset_col                     主图像在原图中的列偏移
	* @param start1                         主图拍摄起始时间 
	* @param end1                           主图拍摄结束时间
	* @param start2                         辅图拍摄起始时间
	* @param end2                           辅图拍摄结束时间
	* @param statevec1                      主星轨道
	* @param statevec2                      辅星轨道
	* @param prf                            脉冲重复频率
	* @param demPath                        SRTM DEM路径
	* @return 成功返回0，否则返回-1
	*/
	int pingpong_MLE(
		Mat& phase_reference,
		Mat& wrapped_phase_low,
		Mat& wrapped_phase_high,
		Mat& outphase,
		Mat& lat_coef,
		Mat& lon_coef,
		double wavelength_low,
		double wavelength_high,
		double nearRangeTime,
		double rangeSpacing,
		double azimuthSpacing,
		int offset_row,
		int offset_col,
		double start1,
		double end1,
		double start2,
		double end2,
		Mat& statevec1,
		Mat& statevec2,
		double prf,
		string demPath
	);
private:
	char error_head[256];

};




#endif // !__SLC_SIMULATOR__H__