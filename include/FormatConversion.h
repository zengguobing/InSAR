#pragma once
#ifndef __FORMATCONVERSION__H__
#include<string>
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"
#include"hdf5.h"
#include"..\include\tinyxml.h"
#define Big2Little32(A) ((uint32_t)(A&0xff000000)>>24|(uint32_t)(A&0x00ff0000)>>8 | (uint32_t)(A&0x0000ff00)<<8|(uint32_t)(A&0x000000ff)<<24)
#define Big2Little16(A) ((uint16_t)(A&0xff00)>>8 | (uint16_t)(A&0x00ff)<<8)
/*********************************************************/
/***************   XML文件参数读写类库    ****************/
/*********************************************************/

class InSAR_API XMLFile
{
public:
	XMLFile();
	~XMLFile();

	/** @brief 创建新的工程文件
	
	@param project_path       工程路径
	@param project_name       工程名
	@param project_version    工程文件版本
	*/
	int XMLFile_creat_new_project(
		const char* project_path,
		const char* project_name,
		const char* project_version
	);
	/** @brief 添加导入TSX节点

	@param node_name      图像名
	@param node_path      图像路径
	*/
	int XMLFile_add_origin(
		const char* node_name,
		const char* node_path
	);
	/** @brief 添加裁剪图像节点

	@param datanode_node  裁剪图像节点名
	@param node_name      裁剪图像名
	@param node_path      裁剪图像路径
	@param Row_offset     行偏移量
	@param Col_offset     列偏移量
	@param lon            中心经度
	@param lat            中心纬度
	@param width          裁剪宽度
	@param height         裁剪高度
	*/
	int XMLFile_add_cut(
		const char* datenode_name,
		const char* node_name,
		const char* node_path,
		int Row_offset,
		int Col_offset,
		double lon, double lat,
		double width, double height
	);
	/** @brief 添加配准图像节点

	@param datanode_node  配准图像节点名
	@param node_name      配准图像名
	@param node_path      配准图像路径
	@param Row_offset     行偏移量
	@param Col_offset     列偏移量
	@param master_index   主图像序号
	@param interp_times   插值倍数（2的n次幂）
	@param block_size     子块尺寸（2的n次幂）
	@param temporal_baseline 时间基线估计
	@param B_effect       垂直基线估计
	@param B_parallel     水平基线估计
	*/
	int XMLFile_add_regis(
		const char* datanode_name,
		const char* node_name,
		const char* node_path,
		int Row_offset,
		int Col_offset,
		int master_index, int interp_times, int block_size,
		const char* temporal_baseline, const char* B_effect, const char* B_parallel
	);
	/** @brief 返回字符串

	@param n			输入整数值
	*/
	string int2str(int n);
	/** @brief 返回整型

	@param s			输入字符串
	*/
	int str2int(const string& s);
	/** @brief 保存XML

	@param save_path    保存路径
	*/
	int XMLFile_save(
		const char* save_path
	);




	/*XML文件加载
	* 参数1：待加载文件名
	*/
	int XMLFile_load(const char* xmlFileName);
	/*
	* 按名称查找节点
	* 参数1：根节点（返回值）
	*/
	int get_root(TiXmlElement*& root);
	/*
	* 按节点查询子节点个数
	* 参数1：节点名
	* 参数2：子节点数量（返回值）
	*/
	int get_children_count(TiXmlElement* pRoot, int* count);
	/*
	* 按名称查找节点
	* 参数1：节点名
	* 参数2：节点名
	* 参数3：节点指针（返回值）
	*/
	int _find_node(TiXmlElement* pRoot, const char* node_name, TiXmlElement*& pnode);
	/*
	* 按名称查找节点
	* 参数1：节点名
	* 参数2：节点指针（返回值）
	*/
	int find_node(const char* node_name, TiXmlElement*& pnode);
	/*
	* 按名称及属性值查找节点
	* 参数1：根节点（遍历用，非文件根节点）
	* 参数2：节点名
	* 参数3：节点属性名（关键属性）
	* 参数4：节点属性值
	* 参数5：节点指针（返回值）
	*/
	int find_node_with_attribute(
		TiXmlElement* pRoot,
		const char* node_name,
		const char* attribute_name,
		const char* attribute_value,
		TiXmlElement*& pnode);
	/*
	* 从XML文件中读出字符串参数
	* 参数1：参数名（节点名）
	* 参数2：参数值（输出）
	*/
	int get_str_para(const char* node_name, string& value);
	/*
	* 从XML文件中读出double类型参数
	* 参数1：参数名（节点名）
	* 参数2：参数值（输出）
	*/
	int get_double_para(const char* node_name, double* value);
	/*
	* 从XML文件中读出整型参数
	* 参数1：参数名（节点名）
	* 参数2：参数值（输出）
	*/
	int get_int_para(const char* node_name, int* value);



	/*************************************************/
	/**********   Sensor-specific functions **********/
	/*************************************************/


	/*
	* 从TerraSAR-X GEOREF.xml文件中读出控制点参数
	* 参数1：控制点参数（输出值，N×6矩阵，每列分别为：经，纬，高，行，列，下视角）
	*/
	int get_gcps_from_TSX(Mat& gcps);
	/*
	* 从TerraSAR-X的主XML文件中读出轨道数据
	* 参数1：轨道参数（输出值，N×7矩阵，每列分别是：GPS时间，位置，速度）
	*/
	int get_stateVec_from_TSX(Mat& stateVec);
	/*
	* 从TerraSAR-X主xml文件中读出多普勒中心频率参数估计系数
	* 参数1：多普勒中心频率（输出值， N×(N_order+2)矩阵，N为多普勒中心估计数，N_order为斜距向多普勒中心拟合阶数, 另外两个是零阶系数和参考点距离向时间）
	*/
	int get_dopplerCentroid_from_TSX(Mat& doppler);


	/*
	* 从sentinel1卫星数据xml文件中读出地面控制点数据
	* 参数1：控制点数据（输出值，N×6矩阵，每列分别为：经，纬，高，行，列，下视角）
	*/
	int get_gcps_from_sentinel(Mat& gcps);
	/*
	* 从sentinel1卫星数据xml文件中读出多普勒中心频率参数
	* 参数1：多普勒中心频率（输出值， N×(N_order+2)矩阵，N为多普勒中心估计数，N_order为斜距向多普勒中心拟合阶数, 另外两个是零阶系数和参考点距离向时间）
	*/
	int get_dopplerCentroid_from_sentinel(Mat& doppler);
	/*
	* 从sentinel1卫星数据xml文件中读出轨道参数
	* 参数1：轨道参数（输出值，N×7矩阵，每列分别是：GPS时间，位置，速度）
	*/
	int get_stateVec_from_sentinel(Mat& stateVec);
	

private:
	char m_xmlFileName[2048];
	TiXmlDocument doc;
	int data_node_count;
	char error_head[256];


	
};



/**************************************************/
/*********           格式转换类库        **********/
/**************************************************/

class InSAR_API FormatConversion
{
public:
	FormatConversion();
	~FormatConversion();

	/*
	* 功能：将字符串格式的UTC时间转换为GPS时间
	* 参数1：UTC时间
	* 参数2：GPS时间
	*/
	int utc2gps(const char* utc_time, double* gps_time);
	/** @brief 创建新的h5文件（若文件已存在则覆盖）
	
	@param filename     文件名
	*/
	int creat_new_h5(const char* filename);
	/*
	* 功能：向h5文件中写入实数矩阵（干涉相位，相干系数等）,input_array可以是16位整型或者double/float型
	* 参数1：文件名
	* 参数2：dataset名
	* 参数3 待写入矩阵
	*/
	int write_array_to_h5(const char* filename, const char* dataset_name, const Mat& input_array);
	/*
	* 功能：从h5文件中读出实数矩阵(读出类型为double型、16位整型或者32位整型)
	* 参数1：文件名
	* 参数2：dataset名
	* 参数3：输出矩阵
	*/
	int read_array_from_h5(const char* filename, const char* dataset_name, Mat& out_array);
	/*
	* 功能：从h5文件中读取矩阵数据子集
	* 参数1：h5文件名
	* 参数2：dataset名
	* 参数3：行偏移量(从0开始)
	* 参数4：列偏移量（从0开始）
	* 参数5：子集行数
	* 参数6：子集列数
	* 参数7：输出子集矩阵
	*/
	int read_subarray_from_h5(
		const char* filename,
		const char* dataset_name,
		int offset_row,
		int offset_col,
		int rows_subarray,
		int cols_subarray,
		Mat& out_array
	);
	/** @brief 向已有的H5文件指定数据集矩阵中的指定位置写入子矩阵
	
	@param h5_filename                   h5文件名
	@param dataset_name                  数据集名称
	@param subarray                      子矩阵数据
	@param offset_row                    行偏移量（从0开始）
	@param offset_col                    列偏移量（从0开始）
	@param rows_subarray                 子集（矩阵）行数
	@param cols_subarray                 子集（矩阵）列数
	@return 成功返回0， 否则返回-1
	*/
	int write_subarray_to_h5(
		const char* h5_filename,
		const char* dataset_name,
		Mat& subarray,
		int offset_row,
		int offset_col,
		int rows_subarray,
		int cols_subarray
	);
	/*
	* 功能：向h5文件写入字符串参数
	* 参数1：文件名
	* 参数2：dataset名
	* 参数3：待写入字符串参数
	*/
	int write_str_to_h5(const char* filename, const char* dataset_name, const char* str);
	/*
	* 功能：从h5文件中读出字符串参数
	* 参数1：文件名
	* 参数2：dataset名
	* 参数3：输出字符串参数
	*/
	int read_str_from_h5(const char* filename, const char* dataset_name, string& string);
	/*
	* 功能：向h5文件写入复图像数据（SLC），如果已经存在则不写入。
	* 参数1：文件名
	* 参数2：复数据
	*/
	int write_slc_to_h5(const char* filename, const ComplexMat& slc);
	/*
	* 功能：从h5文件中读出slc数据
	* 参数1：文件名
	* 参数2：输出slc
	*/
	int read_slc_from_h5(const char* filename, ComplexMat& slc);


	/*------------------------------------------------*/
	/*            TerraSAR-X产品数据导入工具          */
	/*------------------------------------------------*/

	/*
	* 从TerraSAR-X卫星的.cos数据中读出slc数据(不改变类型，仍然是16位整型)
	* 参数1：.cos文件名
	* 参数2：复数据矩阵（输出值）
	*/
	int read_slc_from_TSXcos(const char* filename, ComplexMat& slc);
	/** @brief 将TerraSAR-X卫星数据格式转换为自定义的h5格式
	
	@param cosar_filename                     TerraSAR-X .cos文件名
	@param xml_filename                       TerraSAR-X 主xml文件名
	@param GEOREF_filename                    TerraSAR-X GEOREF.xml文件名
	@param dst_h5_filename                    目标h5文件（若文件已经存在则覆盖）
	@return 成功返回0，否则返回-1
	*/
	int TSX2h5(
		const char* cosar_filename,
		const char* xml_filename,
		const char* GEOREF_filename,
		const char* dst_h5_filename
	);
	/** @brief 将TerraSAR-X卫星数据格式转换为自定义的h5格式
	
	@param xml_filename                       TerraSAR-X 主xml文件名
	@param dst_h5_filename                    目标h5文件
	@return 成功返回0，否则返回-1
	*/
	int TSX2h5(
		const char* xml_filename,
		const char* dst_h5_filename
	);





	/*------------------------------------------------*/
	/*            Sentinel1产品数据导入工具           */
	/*------------------------------------------------*/

	/*
	* 从sentinel1卫星精密轨道数据文件中读出精密轨道数据
	* 参数1：精密轨道数据文件
	* 参数2：粗轨道数据起始时间
	* 参数3：粗轨道数据结束时间
	* 参数4：目标h5文件
	*/
	int read_POD(const char* POD_filename, double start_time, double stop_time, const char* dst_h5_filename);
	/*
	* 功能：从sentinel1卫星数据中读出slc数据（读出数据类型为16位整型）
	* 参数1：sentinel1卫星数据文件名
	* 参数2：xml参数文件
	* 参数3：复矩阵（读出的slc数据）
	*/
	int read_slc_from_Sentinel(const char* filename, const char* xml_filename, ComplexMat& slc);
	/*
	* 功能：sentinel1数据deburst(数据类型为16位整型)
	* 参数1：用于deburst的xml参数文件
	* 参数2：待处理slc数据（原地操作）
	* 参数3：无效行积累标志数据(CV_32S型)
	*/
	int sentinel_deburst(const char* xml_filename, ComplexMat& slc, Mat& sentinel);
	/*
	* 将sentinel1卫星数据格式转换为自定义的h5格式
	* 参数1：tiff格式文件名（储存SLC图像）
	* 参数2：xml文件名
	* 参数3：目标h5文件名
	* 参数4：精密轨道数据文件
	*/
	int sentinel2h5(
		const char* tiff_filename,
		const char* xml_filename,
		const char* dst_h5_filename,
		const char* POD_file = NULL
	);



	/*------------------------------------------------*/
	/*              ALOS1/2产品数据导入工具               */
	/*------------------------------------------------*/

	/*
	* 从ALOS1/2 CEOS格式Level-1.1产品中读取slc数据(float型)
	* 参数1：图像文件
	* 参数2：slc数据
	*/
	int read_slc_from_ALOS(const char* img_file, ComplexMat& slc);
	/*
	* 从ALOS1/2 CEOS格式Level-1.1产品中读取卫星轨道数据
	* 参数1：ALOS LED文件
	* 参数2：轨道数据
	*/
	int read_stateVec_from_ALOS(const char* LED_file, Mat& stateVec);
	/*
	* 从ALOS1/2 CEOS格式Level-1.1产品中读取经纬坐标与图像坐标之间的转换关系系数
	* 参数1：ALOS LED文件
	* 参数2：图像坐标（行、列）与经度之间的转换关系（行列-->经度）
	* 参数3：图像坐标（行、列）与纬度之间的转换关系（行列-->纬度）
	* 参数4：经纬度与行坐标之间的转换关系（经纬度-->行数）
	* 参数5：经纬度与列坐标之间的转换关系（经纬度-->列数）
	*/
	int read_conversion_coefficient_from_ALOS(
		const char* LED_file,
		Mat& lon_coefficient,
		Mat& lat_coefficient,
		Mat& row_coefficient,
		Mat& col_coefficient
	);
	/*
	* 将ALOS1/2 CEOS格式Level-1.1产品数据读出到自定义的hdf5文件中
	* 参数1：ALOS IMG文件
	* 参数2：ALOS LED文件
	* 参数3：自定义h5文件
	*/
	int ALOS2h5(const char* IMG_file, const char* LED_file, const char* dst_h5);
	/** @brief 将原h5文件中的参数信息拷贝到另一个h5中
	
    @param Input_file        原始h5文件
	@param Output_file       输出h5文件
	*/
	int Copy_para_from_h5_2_h5(const char* Input_file, const char* Output_file);

private:
	char error_head[256];
	char parallel_error_head[256];

};

















#endif // !__FORMATCONVERSION__H__
