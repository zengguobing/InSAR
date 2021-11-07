#pragma once
#ifndef _SARCOMM_H_INCLUDED
#define _SARCOMM_H_INCLUDED
#include "globalparam.h"
//#include "../../../Include/REGISStruct.h"

struct Para
{
	double AziDis;      //图像方位向采样间隔
	double TSpace_R;	//图像距离向采样间隔
	double AziStartTime;	//图像方位向起始慢时刻
	double RanStartTime;	//图像距离向起始快时刻
	double frequency;	//发射信号中心载频
};

typedef enum {
	COMPLEX_IMAGE,		/* complex image */
	PHASE_IMAGE,		/* phase image */
	DIFF_PHASE_IMAGE //Differential phase image
} Esip_LayerType;

typedef enum {
	CU4,		//无符号复数,实虚部各4位,先实部后虚部
	CU5,		//无符号复数,实虚部各5位,先实部后虚部
	CU8,		//无符号复数,实虚部各8位,先实部后虚部
	CU16,		//无符号复数,实虚部各16位,先实部后虚部
	CU32,		//无符号复数,实虚部各32位,先实部后虚部
	CS4,		//有符号复数,实虚部各4位,先实部后虚部
	CS5,		//有符号复数,实虚部各5位,先实部后虚部
	CS8,		//有符号复数,实虚部各8位,先实部后虚部
	CS16,		//有符号复数,实虚部各16位,先实部后虚部
	CS32,		//有符号复数,实虚部各32位,先实部后虚部
	S16,		//有符号16位实数
	S32		    //有符号32位实数
} Esip_PixelType;

typedef enum {
	HH,
	HV,
	VH,
	VV
} Esip_Polarization; //极化方式

typedef enum {
	STRIPE,	
	BEAMING
} Esip_SARMode;

typedef enum {
	RIGHT,	
	LEFT
} Esip_Look; //左右视类型

/* -------------------------------------------------------------------- */
/* Class date     时间结构体                                            */                                            
/* -------------------------------------------------------------------- */
#pragma pack (2)
typedef struct {
	long year;			//年
	long month;			//月
	long date;			//日
	long hour;			//时
	long minute;		//分
	double second;      //秒
} Esip_ClassDate;

/* -------------------------------------------------------------------- */
/* Class beam_position 波位参数结构体                                   */                                            
/* -------------------------------------------------------------------- */
typedef struct {
	double prf;			//脉冲重复频率/Hz
	double lookAngle;	//波束中心下视角/rad
	double spreadRatio;	//波束距离向展宽系数k2/无单位
	double pulseWidth;	//发射信号脉宽/s
	double bandWidth;	//发射信号带宽/Hz
	double peakPower;	//发射信号峰值功率/W
	double sampleFrequency;	//接收信号采样频率/Hz
	double Da;	           //天线方位向尺寸/m ，方位向波束宽度为0.886*lamda/Da
	double Dr;	           //天线距离向尺寸/m ，距离向波束宽度为k2*0.886*lamda/Dr
	double spotAngle;	 //聚束波位的前后扫描范围/rad，常用值+-0.6度/工作时长
	double slipFactor;   //滑动因子
} Esip_Classbeam_position;

/* -------------------------------------------------------------------- */
/* Class sar_pos雷达天线相位中心测量数据结构体（地心固连坐标系）        */                                            
/* -------------------------------------------------------------------- */
typedef struct {
	double gpsTime;			//GPS时间/s
	double positionX;	//位置x/m
	double positionY;	//位置y/m
	double positionZ;	//位置z/m
	double velocityX;	//速度vx/ m/s
	double velocityY;	//速度vy/ m/s
	double velocityZ;	//速度vz/ m/s
} Esip_ClassSar_pos;

/* -------------------------------------------------------------------- */
/* Class sar_attitude卫星姿态角测量数据        */                                            
/* -------------------------------------------------------------------- */
typedef struct {
	double gpsTime;			//GPS时间/s
	double rollAngle;	//横滚角
	double pitchAngle;	//俯仰角
	double yawAngle;	//航偏角
} Esip_ClassSar_attitude;

/* -------------------------------------------------------------------- */
/* Class tag 标志点数据结构体        */                                            
/* -------------------------------------------------------------------- */
typedef struct {
	long no;			 //标志点序号
	double pixelAzimuth; //主图像中的方位向象素值/Pixel
	double pixelRange;	 //主图像中的距离向向象素值/Pixel
	double coordinateX;	 //地固坐标x/m
	double coordinateY;	 //地固坐标y/m
	double coordinateZ;	 //地固坐标z/m
} Esip_ClassTag;

/* -------------------------------------------------------------------- */
/*      Structure definitions from eprj.h, with some type               */
/*      simplifications.                                                */
/* -------------------------------------------------------------------- */
typedef struct {
	double x;			/* coordinate x-value */
	double y;			/* coordinate y-value */
} Eprj_Coordinate;


/* -------------------------------------------------------------------- */
/*      Esip_HeaderTag structure               */
/* -------------------------------------------------------------------- */
typedef struct {
	char label[16];			/* label:ESIP_HEADER_TAG */
	GInt32 headerPos;		/* point to Esip_file node */
} Esip_HeaderTag;

/* -------------------------------------------------------------------- */
/*      Esip_HeaderTag structure               */
/* -------------------------------------------------------------------- */
#pragma pack (2)
typedef struct {
	short version;			//version
	short type;		//file type 
	GInt32 freeListPtr;    //the points to list of freed blocks
	GInt32 rootEntryPtr;
	GUInt16 entryHeaderLength;
	GUInt32 dictionaryPtr;
} Esip_File;

/* -------------------------------------------------------------------- */
/*      Esip_Entry structure  头部信息             */
/* -------------------------------------------------------------------- */
#pragma pack (1)
typedef struct {
	GUInt32 nextPtr;	  //next node
	GUInt32 prevPtr;	  //previous node 
	GUInt32 parentPtr;    //parent node
	GUInt32 childPtr;     //child node
	GUInt32 dataPtr;     //data pointer
	GUInt32 datasize;    //data size
	char name[64];   //node name
	char type[32];   //node type
	ModTime modeTime;
} Esip_Entry;


/* -------------------------------------------------------------------- */
/*      Img_Layer structure               */
/* -------------------------------------------------------------------- */
typedef struct {
	double longitudeSceneLocation;			//经度观测区域
	double latitudeSceneLocation;       	//维度观测区域 
	long   azimuthWidth;    //方位向像素点数/像素
	long   rangHeight;  //距离向像素点数/像素
	Esip_LayerType layerType;
	int pixelType;
	long blockWidth;
	long blockheight;
} Esip_ImgLayer;

/* -------------------------------------------------------------------- */
/*      Esip_ClassBlock structure               */
/* -------------------------------------------------------------------- */
typedef struct {
	GUInt32 offset;			//points to the byte location in the file
	long size;  //The number of bytes in the block
} Esip_ClassBlock;

/* -------------------------------------------------------------------- */
/*      ImgBlock structure               */
/* -------------------------------------------------------------------- */
typedef struct {
	long   numvirtualblocks;    //The number of blocks in this layer
	long   numobjectsperblock;  //The number of pixels represented by one block
	Esip_ClassBlock Block[32];
	ModTime modeTime;
} Esip_ImgBlock;

/* -------------------------------------------------------------------- */
/*      LayerInfo structure               */
/* -------------------------------------------------------------------- */
typedef struct {
	short   dataSimulation;    //数据仿真方式:仿真数据1;真实数据2
	double  azimuthSample;     //图像方位向采样间隔/秒
	double  rangeSample;       //图像距离向采样间隔/秒
	double  azimuthStartTime;  //图像方位向起始慢时刻
	double  rangeStartTime;    //图像距离向起始快时刻
} Esip_LayerInfo;

/* -------------------------------------------------------------------- */
/*      SensorInfo 传感器信息节点               */
/* -------------------------------------------------------------------- */
typedef struct {
	Esip_Polarization polarization; //极化方式
	Esip_SARMode sarMode;          //工作模式
	Esip_ClassDate TimeStart;      //回波起始时间
	long   numberOfBeamPosition;  //条带模式波位号
	Esip_Classbeam_position transmitter; //发射星波位参数
	Esip_Classbeam_position receiver; //发接收星波位参数
	double  carrierFrequency;       //发射信号中心载频/Hz
	short signModulationFrequency; //发射信号调频率极性:＋1；－1
	double squintTransmitter; //发射天线斜视角/rad
	double squintReceiver;    //接收天线斜视角/rad
	Esip_Look look;           //左右视类型
	double  resolutionAzimuth;  //理论方位分辨率
	double  resolutionRange;  //理论斜距分辨率
	//Esip_ClassDate timeStartTransmitter; //发射雷达相位中心数据仿真起始时间
	//long sampleNumberTransmitter; //发射雷达相位中心数据数
	//double sampleIntervalTransmitter; //发射雷达相位中心数据采样间隔
	//Esip_ClassSar_pos *antennaPosTransmitter; //发射雷达相位中心数据
	//Esip_ClassDate timeStartReceiver; //接收雷达相位中心数据仿真起始时间
	//long sampleNumberReceiver; //接收雷达相位中心数据数
	//double sampleIntervalReceiver; //接收雷达相位中心数据采样间隔
	//Esip_ClassSar_pos *antennaPosReceiver; //接收雷达相位中心数据
	//long numofDopplerCenterFrequency; //参考多普勒中心频率数据数目
	//double *dopplerCenterFrequency; //参考多普勒中心频率数据
} Esip_SensorInfo;

/* -------------------------------------------------------------------- */
/*      OrbitInfo structure     轨道信息节点          */
/* -------------------------------------------------------------------- */
typedef struct {
	double  earthRotationSpeed;     //地球自转速度
	double  AverGroundVelocity;     //平均地速Vg/m/s
	char  satelliteTransmitter[64]; //发射卫星名
	char  satelliteReceiver[64];    //接收卫星名
	Esip_ClassSar_pos  positionTransmitter;  //发射卫星的星历数据
	Esip_ClassSar_pos  positionReceiver;    //发射卫星的星历数据
	Esip_ClassSar_attitude attitudeTransmitter; //发射卫星姿态测量数据
	Esip_ClassSar_attitude attitudeReceiver; //接收卫星姿态测量数据
} Esip_OrbitInfo;

/* -------------------------------------------------------------------- */
/*      SceneInfo structure     地面信息节点          */
/* -------------------------------------------------------------------- */
typedef struct {
	long  numOfTag;     //地面控制点数
	Esip_ClassTag *dataPointTag; //标志点数据
	char fileNameDEM[64]; //粗DEM数据文件名
} Esip_SceneInfo;


/* -------------------------------------------------------------------- */
//新增附加数据2018.5.30      
/*Esip_Additional_Data  附加数据结构体
/* -------------------------------------------------------------------- */
#if 0
typedef struct {
	Esip_ClassDate  sendRadarClassData;  //发射雷达相位中心数据仿真起始时间
	long  sendRadarXWZXDataNum;          //发射雷达相位中心数据数
	double sendRadarXWZXDataInteral;     //发射雷达相位中心数据采样间隔
	Esip_ClassSar_pos *sendRadarSarPos;  //发射雷达相位中心数据
	Esip_ClassDate  recvRadarClassData;  //接收雷达相位中心数据仿真起始时间
	long  recvRadarXWZXDataNum;          //接收雷达相位中心数据数
	double recvRadarXWZXDataInteral;     //接收雷达相位中心数据采样间隔
	Esip_ClassSar_pos *recvRadarSarPos;  //接收雷达相位中心数据
	long dopplerCenterFrepNum;
	double *usedRangeDopplerCenterFrep;
	long gndCtrlPointNum;
	Esip_ClassTag *tagPointData; //标志点数据
	char fileNameDEM[64]; //粗DEM数据文件名
} Esip_Additional_Data;
#endif
struct Esip_Additional_Data
{
	double AziDis;      //图像方位向采样间隔(m)
	double TSpace_R;	//图像距离向采样间隔(m)
	double AziStartTime;	//图像方位向起始慢时刻
	double RanStartTime;	//图像距离向起始快时刻
	double frequency;	//发射信号中心载频

	Esip_ClassDate  sendRadarClassData;  //发射雷达相位中心数据仿真起始时间
	long  sendRadarXWZXDataNum;          //发射雷达相位中心数据数
	double sendRadarXWZXDataInteral;     //发射雷达相位中心数据采样间隔
	double *sendRadarSarPos = NULL;  //发射雷达相位中心数据
	Esip_ClassDate  recvRadarClassData;  //接收雷达相位中心数据仿真起始时间
	long  recvRadarXWZXDataNum;          //接收雷达相位中心数据数
	double recvRadarXWZXDataInteral;     //接收雷达相位中心数据采样间隔
	double *recvRadarSarPos;  //接收雷达相位中心数据
	long dopplerCenterFrepNum;
	double *usedRangeDopplerCenterFrep;
	long gndCtrlPointNum;
	double *tagPointData; //标志点数据
	char fileNameDEM[512]; //粗DEM数据文件名

} ;

/* -------------------------------------------------------------------- */
/*      ProcessInfo structure     处理信息节点          */
/* -------------------------------------------------------------------- */
struct Esip_ProcessInfo {
	short  numOfSource;     //产生该数据的源数据个数，0为不在本文档
	char *dataSource[64];   //产生该数据的源数据所在的节点名称
	short filerSizeAzimuth; //滤波方位向窗口大小/Pixels
	short filerSizeRange;   //滤波距离向窗口大小/Pixels
	double flatFrequency;   //去平地的平地频率/1/pixel
	short TimePhaseUnwrapping; //相位展开次数
	char cmments[1024];     //对该数据处理者、方法或者过程的描述
} ;

struct ExtImgHeaderInfo {
	long  fileKey;        //密钥,识别文件类型
	long reserved1;      //预留，待扩展
	char filenameTransmitter[64];   //发射雷达相位中心状态数据文件名
	char filenameReceiver[64];   //接收雷达相位中心状态数据文件名
	char filenameTagData[64];   //标志点数据文件名
	char filenameEchoesData[64]; //回波数据文件名
	char filenameTransmitterOfMeasure[64];   //发射雷达相位中心状态测量数据文件名
	char filenameReceiverOfMeasure[64];   //发射雷达相位中心状态测量数据文件名
	char filenameSendSatelliteAttitude[64];   //发射卫星姿态测量文件名
	char filenameRecvSatelliteAttitude[64];   //接收卫星姿态测量文件名
	char nameSendSatellite[64];   //发射卫星名
	char nameRecvSatellite[64];   //接收卫星名
	char dataSimulation[64]; //直接将数据copy这个buf中
	long workmode;   //工作模式
	long noBeamPosition; //波位号，暂时不用
	Esip_ClassDate timeStartEhcoes; //回波仿真起始时间
	long reserved2;  //预留2
	Esip_Classbeam_position transmitter; //发射星波位参数
	Esip_Classbeam_position receiver; //发接收星波位参数
	double  carrierFrequency;  //发射信号中心载频
	long  signModulationFrequency; //发射信号调频率极性
	long reserved3;  //预留3
	double squintTransmitter; //发射天线斜视角/rad
	double squintReceiver;    //接收天线斜视角/rad
	long look;           //左右视类型
	long reserved4;  //预留4
	long   azimuthLength;    //方位向像素点数/像素 //(即nr,对应图像高度) 
	long   rangeLength;     //距离向像素点数/像素 //(即nc,对应图像宽度)
	double  azimuthSample;     //图像方位向采样间隔/秒
	double  rangeSample;       //图像距离向采样间隔/秒
	double  azimuthStartTime;  //图像方位向起始慢时刻
	double  rangeStartTime;    //图像距离向起始快时刻
	double avgVG;  //平均地速
	double  resolutionAzimuth;  //理论方位分辨率
	double  resolutionRange;  //理论斜距分辨率
} ;  //扩展图像文件结构

struct ExtImgHandle_t {
	FILE	*fp;
	//char	*pszPath;
	char    *pszFilename;
	ExtImgHeaderInfo *pszExtHeaderInfo;
	GUInt32	nDataPos;
	GUInt32 nBytesPerPixel;

};

struct HFAInfo_t {
	FILE	*fp;
	char	*pszPath;
	char    *pszFilename; /* sans path */

	GUInt32 nEndOfFile;
	GUInt32 nHeaderPos;
	GUInt32	nRootPos;
	GUInt32	nDictionaryPos;

	GInt16	nEntryHeaderLength;
	GInt32	nVersion;

	int         bTreeDirty;
	//HFAEntry	*poRoot;
	Esip_Entry *poRoot;

	//HFADictionary *poDictionary;
	char	*pszDictionary;

	int		nXSize;
	int		nYSize;

	int		nBands;
	//HFABand	**papoBand;

	void	*pMapInfo;
	void    *pDatum;
	void     *pProParameters;

} ;

//typedef struct demPara
//{
//	//界面参数
//	int m;          //迭代次数，从界面获取
//	//数据
//	int nr;             //解缠相位的行数
//	int nc;             //解缠相位的列数
//	int tag_M;              //主图像地面控制点数
//	Para para_M;       //主图像参数
//	Para para_S;	   //辅图像参数
//}demPara;

typedef HFAInfo_t *HFAHandle;

#endif /* ndef _SARCOMM_H_INCLUDED */