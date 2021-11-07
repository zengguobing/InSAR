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
	/*基于不规则网络的最小费用流解缠算法（用于第二次解缠，不能单独使用）
	* 参数1 已解缠相位（完成了第一次解缠的相位）
	* 参数2 第二次解缠的三角网络节点数组
	* 参数3 第二次解缠的三角网络边数组
	* 参数4 第二次解缠的三角网络边数量
	* 参数5 解缠起始节点
	* 参数6 是否绕过枝切线（默认为否）
	* 参数7 边长阈值（超过阈值则该边不参与解缠）
	*/
	int MCF_second(
		Mat& unwrapped_phase,
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges,
		bool pass = false,
		double thresh = 1000000000.0
	);
	/*调用mcf.exe解Delaunay三角网络最小费用流问题
	* 参数1：最小费用流问题描述文件
	* 参数2：最小费用流算法可执行文件路径
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
private:
	char error_head[256];
	char parallel_error_head[256];

};

Unwrap::Unwrap()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "UNWRAP_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "UNWRAP_DLL_ERROR: error happens when using parallel computing in function: ");
}

Unwrap::~Unwrap()
{
}

#endif