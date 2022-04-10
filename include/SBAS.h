#pragma once
#ifndef __SBAS__H__
#define __SBAS__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"


/*----------------------------------------*/
/*          小基线集三角网络节点          */
/*----------------------------------------*/
class InSAR_API SBAS_node
{
public:
	/*默认构造函数*/
	SBAS_node();
	/*拷贝构造函数*/
	SBAS_node(const SBAS_node& cp);
	/*初始化构造函数1*/
	SBAS_node(int num_neigh_edge);
	~SBAS_node();
	/*赋值函数（深拷贝）*/
	SBAS_node operator = (const SBAS_node& src);


	/*垂直基线(m)*/
	double B_spatial;
	/*时间基线(day)*/
	double B_temporal;
	/*邻接边数量*/
	int num_neigh_edges;
	/*邻接边序号*/
	int* neigh_edges;
	/*节点相位*/
	double phase;
	/*节点横坐标*/
	double x;
	/*节点纵坐标*/
	double y;
	/*是否已解缠(默认未解缠)*/
	bool b_unwrapped;
	/*高程残差*/
	double epsilon_height;
	/*线性形变速率*/
	double deformation_vel;
private:

	
};

/*----------------------------------------*/
/*        小基线集三角网络边结构体        */
/*----------------------------------------*/
struct SBAS_edge
{
	/*是否为边界(默认为否)*/
	bool isBoundry;
	/*边序号*/
	int num;
	/*端点1序号*/
	int end1;
	/*端点2序号*/
	int end2;

	/*积分增益*/
	double gain;
	/*相位梯度（定义为大序号端点 - 小序号端点）*/
	double phase_gradient;
	/*线性形变速度梯度(定义为大序号端点 - 小序号端点)*/
	double delta_deformation_vel;
	/*高程残差梯度(定义为大序号端点 - 小序号端点)*/
	double delta_epsilon_height;


	/*默认构造函数*/
	SBAS_edge()
	{
		isBoundry = false;
		end1 = end2 = num = -1;
		phase_gradient = delta_deformation_vel = delta_epsilon_height = gain = 0.0;
	}
	/*拷贝构造函数*/
	SBAS_edge(const SBAS_edge& cp)
	{
		isBoundry = cp.isBoundry;
		end1 = cp.end1;
		end2 = cp.end2;
		num = cp.num;
		phase_gradient = cp.phase_gradient;
		delta_deformation_vel = cp.delta_deformation_vel;
		delta_epsilon_height = cp.delta_epsilon_height;
		gain = cp.gain;
	}
	/*赋值（深拷贝）*/
	SBAS_edge operator = (const SBAS_edge& cp)
	{
		isBoundry = cp.isBoundry;
		end1 = cp.end1;
		end2 = cp.end2;
		num = cp.num;
		phase_gradient = cp.phase_gradient;
		delta_deformation_vel = cp.delta_deformation_vel;
		delta_epsilon_height = cp.delta_epsilon_height;
		gain = cp.gain;
		return *this;
	}
};

/*----------------------------------------*/
/*      小基线集三角网络三角形结构体      */
/*----------------------------------------*/
struct SBAS_triangle
{
	/*三角形序号*/
	int num;
	/*点1*/
	int p1;
	/*点2*/
	int p2;
	/*点3*/
	int p3;
	/*三角形残差值*/
	double residue;
	/*相邻三角形序号1*/
	int neigh1;
	/*相邻三角形序号2*/
	int neigh2;
	/*相邻三角形序号3*/
	int neigh3;
	/*边1（从1开始）*/
	int edge1;
	/*边2（从1开始）*/
	int edge2;
	/*边3（从1开始）*/
	int edge3;

	/*默认构造函数*/
	SBAS_triangle()
	{
		num = p1 = p2 = p3 = neigh1 = neigh2 = neigh3 = edge1 = edge2 = edge3 = 0;
		residue = 0.0;
	}
	/*拷贝构造函数*/
	SBAS_triangle(const SBAS_triangle& cp)
	{
		this->edge1 = cp.edge1;
		this->edge2 = cp.edge2;
		this->edge3 = cp.edge3;
		this->neigh1 = cp.neigh1;
		this->neigh2 = cp.neigh2;
		this->neigh3 = cp.neigh3;
		this->num = cp.num;
		this->p1 = cp.p1; this->p2 = cp.p2; this->p3 = cp.p3;
		this->residue = cp.residue;
	}
	/*赋值(深拷贝)*/
	SBAS_triangle operator= (const SBAS_triangle& cp)
	{
		this->edge1 = cp.edge1;
		this->edge2 = cp.edge2;
		this->edge3 = cp.edge3;
		this->neigh1 = cp.neigh1;
		this->neigh2 = cp.neigh2;
		this->neigh3 = cp.neigh3;
		this->num = cp.num;
		this->p1 = cp.p1; this->p2 = cp.p2; this->p3 = cp.p3;
		this->residue = cp.residue;
		return *this;
	}
};

/*----------------------------------------*/
/*             小基线集方法类             */
/*----------------------------------------*/
class InSAR_API SBAS
{
public:
	SBAS();
	~SBAS();

	/*@brief 写入时空基线三角网络节点
	* @param nodeFile                   节点文件
	* @param B_temporal                 时间基线（1×n，单位：day）
	* @param B_effect                   空间基线（1×n，单位：m）
	* @return 成功返回0，否则返回-1
	*/
	int write_spatialTemporal_node(
		const char* nodeFile,
		Mat& B_temporal,
		Mat& B_effect
	);
	/*@brief 设置时空基线三角网络节点的时空基线值和节点坐标
	* @param nodes                     三角网络节点数组
	* @param B_temporal                 时间基线（1×n，单位：day）
	* @param B_effect                   空间基线（1×n，单位：m）
	* @return 成功返回0，否则返回-1
	*/
	int set_spatialTemporalBaseline(
		vector<SBAS_node>& nodes,
		Mat& B_temporal,
		Mat& B_effect
	);
	/** @brief 从.edge文件读取Delaunay三角网的边信息
	* @param edge_file               .edge文件
	* @param num_nodes               节点数
	* @param edges                   Delaunay三角网边数组（返回值）
	* @param node_neighbours         每个节点的邻接边数（返回值）
	* @return  成功返回0， 否则返回-1
	*/
	int read_edges(
		const char* edge_file,
		int num_nodes,
		vector<SBAS_edge>& edges,
		vector<int>& node_neighbours
	);
	/*@brief 初始化SBAS_node（只初始化三角网络关系）
	* @param node                   SBAS_node节点数组
	* @param edges                  SBAS_edge边数组
	* @param node_neighbours        每个节点的邻接边数
	* @return  成功返回0， 否则返回-1
	*/
	int init_SBAS_node(
		vector<SBAS_node>& node,
		vector<SBAS_edge>& edges,
		vector<int>& node_neighbours
	);
	/*@brief 初始化SBAS_triangle（只初始化三角网络关系）
	* @param ele_file            .ele文件
	* @param neigh_file          .neigh文件
	* @param triangles           SBAS_triangle三角形结构体数组
	* @param edges               SBAS_edge三角形边结构体数组
	* @param nodes               SBAS_node三角形节点数组
	* @return 成功返回0，否则返回-1
	*/
	int init_SBAS_triangle(
		const char* ele_file,
		const char* neigh_file,
		vector<SBAS_triangle>& triangles,
		vector<SBAS_edge>& edges,
		vector<SBAS_node>& nodes
	);
	/*@brief 计算时空基线三角网络残差点（残差积分方向为逆时针）
	* @param nodes                时空基线三角网络节点数组
	* @param edges                时空基线三角网络边数组
	* @param triangles            时空基线三角网络三角形数组
	* @return 成功返回0，否则返回-1
	*/
	int compute_spatialTemporal_residue(
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangles
	);







	/*@brief 根据时空基线三角网生成差分干涉相位（主图像时间>辅图像时间）
	* @param SLCH5Files          SLC图像文件
	* @param edges               SBAS_edge三角形边结构体数组
	* @param nodes               SBAS_node三角形节点数组
	* @param multilook_az        方位向多视倍数
	* @param multilook_rg        距离向多视倍数
	* @param ifgSavePath         差分相位h5文件保存路径
	* @return 成功返回0，否则返回-1
	*/
	int generate_interferograms(
		vector<string>& SLCH5Files,
		vector<SBAS_edge>& edges,
		vector<SBAS_node>& nodes,
		int multilook_az,
		int multilook_rg,
		const char* ifgSavePath
	);
	/*@brief 估计差分干涉相位数据堆栈相关系数并生成高相干掩膜
	* @param phaseFiles               差分干涉相位数据文件
	* @param wndsize_rg               相关系数估计距离向窗口大小（奇数）
	* @param wndsize_az               相关系数估计方位向窗口大小（奇数）
	* @param coherence_thresh         高相干阈值
	* @param mask                     高相干掩膜（返回值，int型）
	* @return 成功返回0，否则返回-1
	*/
	int generate_high_coherence_mask(
		vector<string>& phaseFiles,
		int wndsize_rg,
		int wndsize_az,
		double coherence_thresh,
		Mat& mask
	);
	/*@brief 根据高相干掩膜矩阵写入三角网节点文件
	* @param mask                     高相干掩膜矩阵(int型)
	* @param nodeFile                 三角网络节点
	* @return 成功返回0，否则返回-1
	*/
	int write_high_coherence_node(
		Mat& mask,
		const char* nodeFile
	);
	/*@brief 根据高相干掩膜矩阵设置三角网络节点的坐标值
	* @param mask                     高相干掩膜矩阵（int型）
	* @param nodes                    高相干三角网络节点数组
	* @return 成功返回0，否则返回-1
	*/
	int set_high_coherence_node_coordinate(
		Mat& mask,
		vector<SBAS_node>& nodes
	);
	/*@brief 根据高相干掩膜矩阵设置三角网络节点相位值
	* @param mask                     高相干掩膜矩阵（int型）
	* @param nodes                    高相干三角网络节点数组
	* @param phase                    干涉相位
	* @return 成功返回0，否则返回-1
	*/
	int set_high_coherence_node_phase(
		Mat& mask,
		vector<SBAS_node>& nodes,
		Mat& phase
	);

private:
	char error_head[256];
};




#endif // !__SBAS__H__
