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
	/*权值*/
	double weight;


	/*默认构造函数*/
	SBAS_edge()
	{
		isBoundry = false;
		end1 = end2 = num = -1;
		phase_gradient = delta_deformation_vel = delta_epsilon_height = gain = 0.0;
		weight = 1.0;
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
		weight = cp.weight;
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
		weight = cp.weight;
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
	/*@brief 将最小费用流问题写入DIMACS文件，准备求解（时空基线三角网）
	* @param DIMACS_file          DIMACS文件
	* @param nodes                时空基线三角网络节点数组
	* @param edges                时空基线三角网络边数组
	* @param triangle             时空基线三角网络三角形数组
	* @return 成功返回0，否则返回-1
	*/
	int writeDIMACS_temporal(
		const char* DIMACS_file,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangle
	);
	/*@brief 根据时空基线三角网生成差分干涉相位（主图像时间>辅图像时间）,生成的差分相位与时空基线三角网络边数组对应
	* @param SLCH5Files          SLC图像文件
	* @param edges               时空基线三角网络边结构体数组
	* @param nodes               时空基线三角网络节点数组
	* @param multilook_az        方位向多视倍数
	* @param multilook_rg        距离向多视倍数
	* @param ifgSavePath         差分相位h5文件保存路径
	* @param b_save_images       是否保存为图片（默认为否）
	* @return 成功返回0，否则返回-1
	*/
	int generate_interferograms(
		vector<string>& SLCH5Files,
		vector<SBAS_edge>& edges,
		vector<SBAS_node>& nodes,
		int multilook_az,
		int multilook_rg,
		const char* ifgSavePath,
		bool b_save_images = false
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
	/*@brief 根据高相干掩膜矩阵设置三角网络节点相位值，并计算相邻节点之间的相位梯度（梯度定义为大序号-小序号）
	* @param mask                     高相干掩膜矩阵（int型）
	* @param nodes                    高相干三角网络节点数组
	* @param edges                    高相干三角网络边数组
	* @param phase                    干涉相位
	* @return 成功返回0，否则返回-1
	*/
	int set_high_coherence_node_phase(
		Mat& mask,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		Mat& phase
	);
	
	/*@brief 将最小费用流问题写入DIMACS文件，准备求解（高相干点三角网）
	* @param DIMACS_file          DIMACS文件
	* @param nodes                高相干点三角网络节点数组
	* @param edges                高相干点三角网络边数组
	* @param triangles            高相干点三角网络三角形数组
	* @return 成功返回0，否则返回-1
	*/
	int writeDIMACS_spatial(
		const char* DIMACS_file,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangles
	);
	/*@brief 从最小费用流求解器解算文件中读取结果(更新梯度)
	* @param MCF_solution_file    最小费用流求解器解算文件
	* @param nodes                三角网络节点数组
	* @param edges                三角网络边数组
	* @param triangles            三角网络三角形数组
	* @param obj_value            最优目标值
	* @param flowcount            总流动值
	* @return 成功返回0，否则返回-1
	*/
	int readDIMACS(
		const char* MCF_solution_file,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangles,
		double* obj_value,
		double* flowcount = NULL
	);
	/*@brief 将差分干涉相位堆栈高相干点之间的相位梯度信息（与高相干点三角网络的边对应）保存在h5文件中
	* @param phaseFiles           差分干涉相位数据堆栈文件
	* @param mask                 高相干点掩膜
	* @param nodes                高相干三角网络节点数组
	* @param edges                高相干三角网络边数组
	* @param dstH5File            保存h5文件
	* @return 成功返回0，否则返回-1
	*/
	int saveGradientStack(
		vector<string>& phaseFiles,
		Mat& mask,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		const char* dstH5File
	);
	/*@brief 估计差分干涉相位数据堆栈相关系数并生成高相干掩膜
	* @param phaseFiles               差分干涉相位数据堆栈文件
	* @param wndsize_rg               相关系数估计距离向窗口大小（奇数）
	* @param wndsize_az               相关系数估计方位向窗口大小（奇数）
	* @param coherence_thresh         高相干点相关系数阈值（0~1）
	* @param count_thresh             掩膜筛选阈值（0~1，若某点相关系数大于阈值的图幅数大于count_thresh×图幅数，则该点为高相干点）
	* @param mask                     高相干掩膜（返回值，int型）
	* @return 成功返回0，否则返回-1
	*/
	int generate_high_coherence_mask(
		vector<string>& phaseFiles,
		int wndsize_rg,
		int wndsize_az,
		double coherence_thresh,
		double count_thresh,
		Mat& mask
	);
	/*@brief 洪水淹没法高相干点积分解缠（三角网络边梯度已经用MCF求解过）
	* @param nodes                   高相干点三角网络节点
	* @param edges                   高相干点三角网络边
	* @param start                   解缠起始节点
	* @param b_zero_start            解缠起始点值是否设置为0
	* @return 成功返回0，否则返回-1
	*/
	int floodFillUnwrap(
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		int start,
		bool b_zero_start = false
	);
	/*@brief 用相干系数设置高相干点三角网络边的权重
	* @param coherence               相干系数
	* @param nodes                   高相干点三角网络节点
	* @param edges                   高相干点三角网络边
	* @return 成功返回0，否则返回-1
	*/
	int set_weight_by_coherence(
		Mat& coherence,
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges
	);
	/*@brief 从三角网络节点中获取解缠相位
	* @param nodes                   三角网络节点
	* @param phase                   相位（返回值，inplace操作）
	* @return 成功返回0，否则返回-1
	*/
	int retrieve_unwrapped_phase(
		vector<SBAS_node>& nodes,
		Mat& phase
	);
	/*@brief 根据原始差分相位，计算高相干点三角网络残差点（残差积分方向为逆时针）
	* @param nodes                   高相干三角网络节点
	* @param edges                   高相干三角网络边
	* @param triangles               高相干三角网络三角形
	* @return 成功返回0，否则返回-1
	*/
	int compute_high_coherence_residue(
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangle
	);
	/*@brief 根据时间维解缠梯度，计算高相干点三角网络残差点（残差积分方向为逆时针）
	* @param nodes                   高相干三角网络节点
	* @param edges                   高相干三角网络边
	* @param triangles               高相干三角网络三角形
	* @return 成功返回0，否则返回-1
	*/
	int compute_high_coherence_residue_by_gradient(
		vector<SBAS_node>& nodes,
		vector<SBAS_edge>& edges,
		vector<SBAS_triangle>& triangle
	);
	/*@brief 检查三角网络中残差点数
	* @param triangles              三角网络三角形
	* @param num                    残差点数（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int residue_num(
		vector<SBAS_triangle>& triangles,
		int* num
	);
	/*@brief 根据时空基线分布和相应阈值得到干涉组合矩阵，并求出相应的时空基线
	* @param spatial                        空间基线(1×n，单位：m)
	* @param temporal                       时间基线（1×n，单位：day）
	* @param spatial_thresh                 空间基线阈值（m）
	* @param temporal_thresh                时间基线阈值（year）
	* @param formation_matrix               干涉组合矩阵（返回值，int型，n×n）
	* @param spatial_baseline               空间基线矩阵（与干涉组合矩阵对应，单位：m）
	* @param temporal_baseline              时间基线矩阵（与干涉组合矩阵对应，单位：year）
	* @return 成功返回0，否则返回-1
	*/
	int get_formation_matrix(
		Mat& spatial,
		Mat& temporal,
		double spatial_thresh,
		double temporal_thresh,
		Mat& formation_matrix,
		Mat& spatial_baseline,
		Mat& temporal_baseline
	);
	/*@brief 根据干涉组合矩阵生成差分干涉相位（主图像时间>辅图像时间）
	* @param SLCH5Files          SLC图像文件
	* @param formation_matrix    干涉相位组合矩阵
	* @param spatial_baseline    空间基线（与干涉组合矩阵对应）
	* @param temporal_baseline   时间基线（与干涉组合矩阵对应）
	* @param multilook_az        方位向多视倍数
	* @param multilook_rg        距离向多视倍数
	* @param ifgSavePath         差分相位h5文件保存路径
	* @param b_save_images       是否保存为图片（默认为否）
	* @param Goldstein_alpha     Goldstein滤波强度（默认为0.8）
	* @return 成功返回0，否则返回-1
	*/
	int generate_interferograms(
		vector<string>& SLCH5Files,
		Mat& formation_matrix,
		Mat& spatial_baseline,
		Mat& temporal_baseline,
		int multilook_az,
		int multilook_rg,
		const char* ifgSavePath,
		bool b_save_images = false,
		double alpha = 0.8
	);
	/*@brief 计算时间相关系数（temporal_coherence），评估时间序列估计效果
	* @param estimated_phase_series              时间序列差分相位估计结果(n×1)
	* @param phase_series                        原始差分相位(n×1)
	* @param temporal_coherence                  时间相关系数（返回值）
	* @return 成功返回0，否则返回-1
	*/
	int compute_temporal_coherence(
		Mat& estimated_phase_series,
		Mat& phase_series,
		double* temporal_coherence
	);
	/*@brief 基于同质像元识别的自适应多视干涉相位生成（分块读取、计算、储存）
	@param coregis_slc_files              配准并去斜后SAR图像数据堆栈（文件）
	@param ifgSavePath                    差分相位h5文件保存路径
	@param formation_matrix               SBAS干涉组合矩阵
	@param spatial_baseline               空间基线（与干涉组合矩阵对应）
	@param temporal_baseline              时间基线（与干涉组合矩阵对应）
	@param blocksize_row                  子块尺寸（行，必须大于同质检验搜索窗口半径）
	@param blocksize_col                  子块尺寸（列，必须大于同质检验搜索窗口半径）
	@param out_mask                       掩膜输出（int型，标记经过EVD法估计的像素点，与参数thresh_c1_to_c2有关）
	@param b_coh_est                      是否估计相关系数（默认是）
	@param homogeneous_test_wnd           同质检验搜索窗口大小（奇数，homogeneous_test_wnd×homogeneous_test_wnd， 默认为21×21）
	@param thresh_c1_to_c2                协方差矩阵第2特征值与第1特征值比值阈值（0-1之间，默认为0.7，小于阈值则进行EVD估计）
	@param b_normalize                    协方差矩阵是否归一化（默认是）
	@param b_save_images                  是否将干涉相位保存为图片（默认是）
	*/
	int adaptive_multilooking(
		vector<string>& coregis_slc_files,
		const char* ifgSavePath,
		Mat& formation_matrix,
		Mat& spatial_baseline,
		Mat& temporal_baseline,
		int blocksize_row,
		int blocksize_col,
		Mat& out_mask,
		bool b_coh_est = true,
		int homogeneous_test_wnd = 21,
		double thresh_c1_to_c2 = 0.7,
		bool b_normalize = true,
		bool b_save_images = true
	);
	/*@brief 轨道精炼重去平（一阶拟合）
	* @param unwrapped_phase                       解缠相位
	* @param mask                                  高相干点掩膜
	* @param coherence                             相关系数
	* @param coh_thresh                            选控制点相关系数阈值
	* @param reference                             参考点序号
	* @return 成功返回0，否则返回-1
	*/
	int refinement_and_reflattening(
		Mat& unwrapped_phase,
		Mat& mask,
		Mat& coherence,
		double coh_thresh,
		int reference
	);
private:
	char error_head[256];
};




#endif // !__SBAS__H__
