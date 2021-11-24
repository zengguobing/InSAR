#pragma once
#ifndef __UNWRAP__H__
#define __UNWRAP__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"
#include"..\include\Utils.h"



class InSAR_API Unwrap
{
public:
	Unwrap();
	~Unwrap();
	/*基于规则网络的最小费用流相位解缠算法
	  参数1：待解缠相位
	  参数2：解缠相位（返回值）
	  参数3：干涉相位相干系数
	  参数4：干涉相位残差点
	  参数5：最小费用流问题描述文件
	  参数6：最小费用流算法可执行文件路径
	*/
	int MCF(
		Mat& wrapped_phase,
		Mat& unwrapped_phase,
		Mat& coherence, Mat& residue,
		const char* MCF_problem_file,
		const char* MCF_EXE_PATH
	);
	/*基于不规则网络的最小费用流解缠算法
	* 参数1 待解缠相位
	* 参数2 解缠相位（返回值）
	* 参数3 mask
	* 参数4 Delaunay三角网络节点
	* 参数5 Delaunay三角网边
	* 参数6 Delaunay三角网边数量
	* 参数7 解缠起始节点
	* 参数8 是否绕过枝切线（默认为否）
	* 参数9 边长阈值（超过阈值则该边不参与解缠）
	*/
	int MCF(
		Mat& wrapped_phase,
		Mat& unwrapped_phase,
		Mat& mask,
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges,
		int start,
		bool pass = false,
		double thresh = 10000.0
	);
	/** @brief 基于Delaunay三角网络的最小费用流解缠算法
	
	@param wrapped_phase                           待解缠相位
	@param unwrapped_phase                         解缠相位（返回值）
	@param out_mask                                解缠结果掩膜（返回值）
	@param mask                                    需要解缠的像素点掩膜
	@param nodes                                   Delaunay三角网络节点
	@param edges                                   Delaunay三角网边
	@param start                                   解缠起始点
	@param pass                                    是否绕过枝切线（默认为否）
	@param thresh                                  边长阈值（超过阈值则该边不参与解缠）
	*/
	int MCF(
		const Mat& wrapped_phase,
		Mat& unwrapped_phase,
		Mat& out_mask,
		const Mat& mask,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges,
		int start,
		bool pass = false,
		double thresh = 10000.0
	);
	/*基于不规则网络的最小费用流解缠算法（用于第二次解缠，不能单独使用）
	* 参数1 已解缠相位（完成了第一次解缠的相位）
	* 参数2 第二次解缠的三角网络节点数组
	* 参数3 第二次解缠的三角网络边数组
	* 参数4 第二次解缠的三角网络边数量
	* 参数5 是否绕过枝切线（默认为否）
	* 参数6 边长阈值（超过阈值则该边不参与解缠）
	*/
	int MCF_second(
		Mat& unwrapped_phase,
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges,
		bool pass = false,
		double thresh = 1000000000.0
	);
	/** @brief 调用mcf.exe解Delaunay三角网络最小费用流问题
	
	@param MCF_problem_file                   最小费用流问题描述文件
	@param MCF_EXE_PATH                       最小费用流算法可执行文件路径
	@return 成功返回0，否则返回-1
	*/
	int mcf_delaunay(
		const char* MCF_problem_file,
		const char* MCF_EXE_PATH
	);
	/*结合质量图和最小费用流的解缠法（此法需要提前计算每条边的质量值）
	* 参数1 待解缠相位
	* 参数2 解缠相位（返回值）
	* 参数3 mask
	* 参数4 Delaunay三角网络节点
	* 参数5 Delaunay三角网边
	* 参数6 Delaunay三角网边数量
	* 参数7 解缠起始节点
	* 参数8 是否绕过枝切线（默认为否）
	* 参数9 边长阈值（超过阈值则该边不参与解缠）
	*/
	int QualityMap_MCF(
		Mat& wrapped_phase,
		Mat& unwrapped_phase,
		Mat& mask,
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges,
		int start,
		bool pass = false,
		double thresh = 10000.0
	);
	/** @brief 结合质量图和最小费用流的解缠法第一步（绕过残差点的质量图法）
	
	@param wrapped_phase                      待解缠相位
	@param unwrapped_phase                    解缠相位
	@param out_mask                           已解缠像元掩膜（返回值）
	@param nodes                              Delaunay三角网络节点
	@param edges                              Delaunay三角网络边
	@param distance_thresh                    边长阈值（超过阈值则该边不参与解缠）
	@param pass                               是否绕过残差点（默认是）
	@return 成功返回0，否则返回-1
	*/
	int _QualityGuided_MCF_1(
		const Mat& wrapped_phase,
		Mat& unwrapped_phase,
		Mat& out_mask,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges,
		double distance_thresh = 1.5,
		bool pass = true
	);
	/** @brief 结合质量图和最小费用流的解缠法第二步（第一步未解缠的使用最小费用流法解缠）
	
	@param unwrapped_phase                       已解缠相位（完成了第一次解缠的相位）
	@param nodes                                 未解缠像素点组成的Delaunay三角网络节点
	@param edges                                 未解缠像素点组成的Delaunay三角网络边
	@param distance_thresh                       边长阈值（超过该值不通过该边进行解缠处理）
	@return 成功返回0，否则返回-1
	*/
	int _QualityGuided_MCF_2(
		Mat& unwrapped_phase,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges,
		double distance_thresh = 1.5
	);
	/** @brief 结合质量图和最小费用流的解缠法
	 
	@param wrapped_phase                       待解缠相位
	@param unwrapped_phase                     解缠相位（返回值）
	@param coherence_thresh                    相关系数划分阈值（0~1之间，高相关像素采用质量图法解缠，低相关像素采用最小费用流法）
	@param distance_thresh                     边长阈值（超过此阈值不参与解缠和残差点计算）
	@param tmp_path                            中间结果保存路径
	@param EXE_path                            最小费用流求解器/Delaunay三角网生成器路径
	*/
	int QualityGuided_MCF(
		const Mat& wrapped_phase,
		Mat& unwrapped_phase,
		double coherence_thresh,
		double distance_thresh,
		const char* tmp_path,
		const char* EXE_path
	);
	/** @brief 统计费用流法解缠（SNAPHU）
	
	@param wrapped_phase_file                            缠绕相位文件（h5）
	@param unwrapped_phase                               解缠相位
	@param project_path                                  工程路径
	@param tmp_folder                                    中间结果保存路径
	@param exe_path                                      snaphu.exe路径
	@return 成功返回0，否则返回-1
	*/
	int snaphu(
		const char* wrapped_phase_file,
		Mat& unwrapped_phase,
		const char* project_path,
		const char* tmp_folder,
		const char* exe_path
	);

private:
	char error_head[256];
	char parallel_error_head[256];

};



#endif