#pragma once
#ifndef __UTILS__H__
#define __UTILS__H__
#include"..\include\Package.h"
#include"..\include\ComplexMat.h"
#include<fstream>
#include<iostream>
#include<string>
#include"..\include\sar_comm.h"
#define INPUTMAXSIZE 1024
#define PI 3.141592653589793238
#define VEL_C 299792458.0
/*********************************************************/
/*                Delaunay三角网 节点类                  */
/*********************************************************/
class InSAR_API tri_node
{
public:
	/*默认构造函数*/
	tri_node();
	/*拷贝构造函数*/
	tri_node(const tri_node& node);
	/*构造函数
	* 参数1 节点行数
	* 参数2 节点列数
	* 参数3 节点邻接边数
	* 参数4 节点相位
	*/
	tri_node(int, int, int, double);
	~tri_node();
	/*赋值函数（深拷贝赋值）*/
	tri_node operator = (const tri_node& src);
	/*获取节点相位
	* 参数1 相位指针（返回值）
	*/
	int get_phase(double* phi) const;
	/*获取节点行列坐标
	* 参数1 行序号
	* 参数2 列序号
	*/
	int get_pos(int* rows, int* cols) const;
	/*节点相位赋值
	* 参数1 输入相位
	*/
	int set_phase(double phi);
	/*获取邻接边指针
	* 参数1 指向邻接边指针的指针（返回值）
	* 参数2 邻接边个数指针（返回值）
	*/
	int get_neigh_ptr(long** ptr2ptr, int* num) const;
	/*改变解缠状态
	* 参数1 是否已经解缠
	*/
	int set_status(bool b_unwrapped);
	/*改变平衡状态
	* 参数1 是否属于残差平衡三角形
	*/
	int set_balance(bool b_balanced);
	/*打印邻接边序号
	* 
	*/
	int print_neighbour() const;
	/*获取邻接边个数
	* 参数1 邻接边个数指针
	*/
	int get_num_neigh(int* num_neigh) const;
	/*获取与另一节点的距离
	* 参数1 另一节点
	* 参数2 距离
	*/
	int get_distance(tri_node node, double* distance) const;
	/*获取解缠状态
	* 返回值（是否已解缠）
	*/
	bool get_status() const;
	/*获取平衡状态
	* 返回值（是否平衡）
	*/
	bool get_balance() const;
	/*返回是否节点属于残差三角形
	*/
	bool is_residue_node() const;
	/*设置节点是否属于残差节点
	*/
	int set_residue(bool b_res);
	/*获取形变速率*/
	double get_vel() const;
	/*获取高程误差*/
	double get_height() const;
	/*设置形变速率*/
	int set_vel(double vel);
	/*设置高程误差*/
	int set_height(double height);

private:

	/*****************InSAR处理变量*******************/

	/*是否已解缠(默认未解缠)*/
	bool b_unwrapped;
	/*是否属于残差节点*/
	bool b_residue;
	/*是否属于平衡三角形的顶点（默认为是），同时在PS-InSAR中充当是否节点被丢弃的标志(为true表示不被丢弃， 为false表示被丢弃)*/
	bool b_balanced;
	/*节点行数（起始值为0）*/
	int rows;
	/*节点列数（起始值为0）*/
	int cols;
	/*节点邻接边数*/
	int num_neigh_edges;
	/*节点相位*/
	double phase;
	/*节点邻接边序号*/
	long* neigh_edges;

	/*****************PS-InSAR处理变量*******************/
	
	/*形变速率*/
	double vel;
	/*高程误差*/
	double epsilon_height;
};

/*********************************************************/
/*             Delaunay三角网 三角形结构体               */
/*********************************************************/
struct triangle
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
	triangle()
	{
		num = p1 = p2 = p3 = neigh1 = neigh2 = neigh3 = edge1 = edge2 = edge3 = 0;
		residue = 0.0;
	}
	/*拷贝构造函数*/
	triangle(const triangle& cp)
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
	triangle operator= (const triangle& cp)
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


/*********************************************************/
/*             Delaunay三角网 三角形边结构体             */
/*********************************************************/
struct tri_edge
{
	/**********InSAR变量**********/

	/*积分增益（序号从小到大为正）*/
	double gain;
	/*相位质量（用于质量图法解缠）*/
	double quality;
	/*边序列号*/
	int num;
	/*端点1*/
	int end1;
	/*端点2*/
	int end2;
	/*残差边标志*/
	bool isResidueEdge;
	/*网络边界标志*/
	bool isBoundry;


	/**********PS_InSAR变量**********/

	/*线性形变速度差系数（4 * pi / lambda * Ti）*/
	//double coef_delta_vel;
	/*高程误差系数（4 * pi * bperp_i / lambda / R_i / sin_theta_i  ）*/
	//double coef_delta_height;
	/*线性形变速度差(定义为大坐标 - 小坐标)*/
	double delta_vel;
	/*高程误差(定义为大坐标 - 小坐标)*/
	double delta_height;
	/*模型相干系数*/
	double MC;
	/*端点相位差（相位差定义为：大序号端点减小序号端点）*/
	double phase_diff;

	/*默认构造函数*/
	tri_edge() {
		gain = 0.0;
		quality = 0.0;
		num = 0;
		end1 = 0; end2 = 0;
		isResidueEdge = false;
		isBoundry = false;
		delta_vel = 0.0;
		delta_height = 0.0; MC = 0.0; phase_diff = 0.0;
	}
	/*拷贝构造函数*/
	tri_edge(const tri_edge& cp)
	{
		gain = cp.gain;
		quality = cp.quality;
		num = cp.num;
		end1 = cp.end1; end2 = cp.end2;
		isResidueEdge = cp.isResidueEdge;
		isBoundry = cp.isBoundry;
		delta_vel = cp.delta_vel;
		delta_height = cp.delta_height; MC = cp.MC; phase_diff = cp.phase_diff;
	}
	/*赋值函数（深拷贝赋值）*/
	tri_edge operator = (const tri_edge& cp)
	{
		gain = cp.gain;
		quality = cp.quality;
		num = cp.num;
		end1 = cp.end1; end2 = cp.end2;
		isResidueEdge = cp.isResidueEdge;
		isBoundry = cp.isBoundry;
		delta_vel = cp.delta_vel;
		delta_height = cp.delta_height; MC = cp.MC; phase_diff = cp.phase_diff;
		return *this;
	}
};

/*********************************************************/
/*          Delaunay三角网 三角形边序列号结构体          */
/*********************************************************/

struct edge_index
{
	double quality;
	int num;
	edge_index() { num = 0; quality = 0.0; }
	friend bool operator < (struct edge_index a, struct edge_index b)
	{
		return a.quality > b.quality;
	}
	
};

/*-------------------------------------------------------*/
/*                   规则网格节点结构体                  */
/*-------------------------------------------------------*/

struct node_index
{
	/*节点行数（从0开始）*/
	int row;
	/*节点列数（从0开始）*/
	int col;
	/*默认构造函数*/
	node_index()
	{
		row = 0; col = 0;
	}
	/*拷贝构造函数*/
	node_index(const node_index& cp)
	{
		this->row = cp.row; this->col = cp.col;
	}
	/*赋值函数*/
	node_index operator = (const node_index& cp)
	{
		this->row = cp.row; this->col = cp.col;
		return *this;
	}
};


/*********************************************************/
/*               干涉SAR处理基本函数类库                 */
/*********************************************************/
class InSAR_API Utils
{
public:
	Utils();
	~Utils();
	/*计算矩阵梯度
	 参数1 源矩阵
	 参数2 行方向梯度（返回值）
	 参数3 列方向梯度（返回值）
	 参数4 是否补零使得梯度矩阵和源矩阵大小相同（默认补零）
	*/
	int diff(Mat& Src, Mat& diff_1, Mat& diff_2, bool same = true);
	/*计算干涉相位
	 参数1 主图像（复）
	 参数2 辅图像（复）
	 参数3 干涉相位（返回值）
	*/
	int generate_phase(const ComplexMat& Master, const ComplexMat& Slave, Mat& phase);

	/** @brief 最大似然相干估算器
	 
	@param master_image                       主图像（复）
	@param slave_image                        辅图像（复）
	@param coherence                          相干系数（返回值）
	@return 成功返回0，否则返回-1
	*/
	int real_coherence(ComplexMat& master_image, ComplexMat& slave_image, Mat& coherence);
	/** @brief 最大似然相干估算器（带估计窗口尺寸接口）
	
	@param master_image                       主图像（复）
	@param slave_image                        辅图像（复）
	@param est_wndsize_rg                     估计窗口距离向尺寸（奇数）
	@param est_wndsize_az                     估计窗口方位向尺寸（奇数）
	@param coherence                          相干系数（返回值）
	*/
	int real_coherence(
		const ComplexMat& master_image,
		const ComplexMat& slave_image,
		int est_wndsize_rg,
		int est_wndsize_az,
		Mat& coherence
	);
	/** @brief 频率无关相干估算器
	
	 @param master_image                        主图像（复）
	 @param slave_image                         辅图像（复）
	 @param coherence                           相干系数（返回值）
	 @return 成功返回0，否则返回-1
	*/
	int complex_coherence(ComplexMat& master_image, ComplexMat& slave_image, Mat& coherence);
	/** @brief 频率无关相干估算器（带估计窗口尺寸接口）
	
	@param master_image                         主图像
	@param slave_image                          辅图像
	@param est_wndsize_rg                       估计窗口距离向尺寸（奇数）
	@param est_wndsize_az                       估计窗口方位向尺寸（奇数）
	@param coherence                            相关系数（返回值）
	@return 成功返回0，否则返回-1
	*/
	int complex_coherence(
		const ComplexMat& master_image,
		const ComplexMat& slave_image,
		int est_wndsize_rg,
		int est_wndsize_az,
		Mat& coherence
	);
	/** @brief 根据干涉相位求相关系数
	@param phase                          输入相位
	@param coherence                      相关系数（返回值）
	@return 成功返回0，否则返回-1
	*/
	int phase_coherence(Mat& phase, Mat& coherence);
	/** @brief 根据干涉相位求相关系数（带估计窗口尺寸接口）
	
	@param phase                          输入相位
	@param est_wndsize_rg                 估计窗口距离向尺寸（奇数）
	@param est_wndsize_az                 估计窗口方位向尺寸（奇数）
	@param coherence                      相关系数（返回值）
	@return 成功返回0，否则返回-1
	*/
	int phase_coherence(
		const Mat& phase,
		int est_wndsize_rg,
		int est_wndsize_az,
		Mat& coherence
	);
	/*求解相位导数方差
	* 参数1 干涉相位
	* 参数2 相位导数方差（返回值）
	* 参数3 计算窗口大小（奇数）
	*/
	int phase_derivatives_variance(Mat& phase, Mat& phase_derivatives_variance, int wndsize = 3);
	/*最大可积距离
	* 参数1 原始相位
	* 参数2 最大可积距离（返回值）
	* 参数3 保守值（最大积分距离不能超过该值）
	*/
	int max_integrable_distance(Mat& phase, Mat& max_integrable_distance, double conservative_thresh = 20.0);
	/*FFTSHIFT
	 参数1 待fftshift的矩阵（原地进行fftshift操作）
	*/
	int fftshift(Mat& matrix);

	/*计算干涉相位图的残差值（点）
	 参数1 干涉相位
	 参数2 残差点矩阵（返回值）
	*/
	int residue(Mat& phase, Mat& residue);
	/*计算Delaunay三角网络的残差值（并且标注残差边和残差节点,便于解缠时避开）
	* 参数1 Delaunay三角网三角形结构体数组
	* 参数2 Delaunay三角网三角形数量
	* 参数3 Delaunay三角网节点数组
	* 参数4 Delaunay三角网边结构体数组
	* 参数5 Delaunay三角网边数量
	*/
	int residue(triangle* tri, int num_triangle, vector<tri_node>& nodes, tri_edge* edges, int num_edges);
	/** @brief 计算Delaunay三角网络的残差值（并且标注残差边和残差节点）
	
	@param triangle                              Delaunay三角网三角形结构体数组
	@param nodes                                 Delaunay三角网节点数组
	@param edges                                 Delaunay三角网边结构体数组
	@param distance_thresh                       边长度阈值（超过此阈值不参与残差点计算）
	@return 成功返回0，否则返回-1
	*/
	int residue(
		vector<triangle>& triangle,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges,
		double distance_thresh
	);
	/*计算mask（筛选高质量点）
	* 参数1 相关系数矩阵
	* 参数2 mask举矩阵（返回值）
	* 参数3 窗口半径
	* 参数4 阈值
	*/
	int gen_mask(Mat& coherence, Mat& mask, int wnd_size, double thresh);
	/*计算mask（筛选高质量点）
	* 参数1 相关系数矩阵
	* 参数2 相位导数方差
	* 参数3 mask举矩阵（返回值）
	* 参数4 窗口半径
	* 参数5 相关系数阈值
	* 参数6 相位导数方差阈值
	*/
	int gen_mask(
		Mat& coherence,
		Mat& phase_derivatives,
		Mat& mask, int wnd_size,
		double coh_thresh,
		double phase_derivative_thresh
	);
	/*根据设定阈值筛选残差点
	* 参数1 原始残差点矩阵
	* 参数2 筛选后残差点矩阵
	* 参数3 筛选阈值（大于0）
	* 参数4 残差点个数
	*/
	int residue_sift(Mat& residue_src, Mat& residue_dst, double thresh, long* num_residue);
	/*缠绕相位至（-pi,pi）
	 参数1 待缠绕相位
	 参数2 缠绕后的相位（返回值）
	*/
	int wrap(Mat& Src, Mat& Dst);

	/*按行或列累计积分
	 参数1 待积分数据
	 参数2 积分方向(dim = 1,按列计算 dim = 2,按行计算)
	*/
	int cumsum(Mat& phase, int dim);
	/*叉乘运算（三维）
	* 参数1 向量一(n * 3)
	* 参数2 向量二(n * 3)
	* 参数3 输出
	*/
	int cross(Mat& vec1, Mat& vec2, Mat& out);

	/*写入DIMACS文件（描述最小费用问题）
	 参数1 目标文件名
	 参数2 残差点矩阵
	 参数3 相干系数矩阵
	 参数4 残差点阈值(大于0)
	*/
	int write_DIMACS(const char* DIMACS_file_problem, Mat& residue, Mat& coherence, double thresh);
	/*写入DIMACS文件（描述最小费用问题，不规则三角网络）
	* 参数1 目标文件名
	* 参数2 Delaunay三角形结构体数组
	* 参数3 Delaunay三角形数量
	* 参数4 Delaunay三角网节点数组
	* 参数5 Delaunay三角网边结构体数组
	* 参数6 Delaunay三角网边数量
	* 参数7 每个节点的费用
	*/
	int write_DIMACS(
		const char* DIMACS_file_problem,
		triangle* tri,
		int num_triangle,
		vector<tri_node>& nodes,
		tri_edge* edges,
		long num_edges,
		Mat& cost
	);
	/** @brief 写入DIMACS文件（描述最小费用问题，Delaunay三角网络）
	
	@param DIMACS_file_problem                         目标DIMACS文件
	@param triange                                     Delaunay三角形结构体数组
	@param nodes                                       Delaunay三角网节点数组
	@param edges                                       Delaunay三角网边结构体数组
	@param cost                                        每个节点的费用
	@return 成功返回0，否则返回-1
	*/
	int write_DIMACS(
		const char* DIMACS_file_problem,
		vector<triangle>& triangle,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges,
		const Mat& cost
	);
	/*读取DIMACS文件（获取求解器求解结果）
	 参数1 最小费用流问题解文件
	 参数2 枝切路径1
	 参数3 枝切路径2
	 参数4 干涉相位图像行数
	 参数5 干涉相位图像列数
	*/
	int read_DIMACS(const char* DIMACS_file_solution, Mat& k1, Mat& k2, int rows, int cols);
	/*读取DIMACS文件（获取求解器求解结果）
	* 参数1 最小费用流问题解文件
	* 参数2 Delaunay三角网边结构体数组
	* 参数3 Delaunay三角网边数量
	* 参数4 Delaunay三角网节点数组
	* 参数5 Delaunay三角网三角形数组
	* 参数6 Delaunay三角网三角形数量
	*/
	int read_DIMACS(
		const char* DIMACS_file_solution,
		tri_edge* edges,
		int num_edges,
		vector<tri_node>& nodes,
		triangle* tri,
		int num_triangle
	);
	/** @brief 读取DIMACS文件（获取求解器求解结果）
	
	@param DIMACS_file_solution                         最小费用流问题解文件
	@param edges                                        Delaunay三角网边结构体数组
	@param nodes                                        Delaunay三角网节点数组
	@param triangle                                     Delaunay三角网三角形数组
	@param return 成功返回0，否则返回-1
	*/
	int read_DIMACS(
		const char* DIMACS_file_solution,
		vector<tri_edge>& edges,
		vector<tri_node>& nodes,
		vector<triangle>& triangle
	);
	/*将OpenCV Mat数据以二进制方式写入目标文件
	* 参数1 目标文件名
	* 参数2 待写入数据
	*/
	int cvmat2bin(const char* Dst_file, Mat& Src);
	/*从二进制文件中读数据，并将数据转换成OpenCV Mat格式
	* 参数1 二进制文件
	* 参数2 目标矩阵
	*/
	int bin2cvmat(const char* Src_file, Mat& Dst);
	/*InSAR多视处理（配准之后进行， 改变图像尺寸）
	* 参数1 主图像（SLC）
	* 参数2 辅图像（SLC）
	* 参数3 多视相位
	* 参数4 多视倍数（大于1）
	*/
	int multilook(ComplexMat& Master, ComplexMat& Slave, Mat& phase, int multilook_times);
	/** @brief InSAR多视处理（不改变图像尺寸）
	
	@param master_slc                    主图像
	@param slave_slc                     辅图像
	@param multilook_rg                  距离向多视倍数
	@param multilook_az                  方位向多视倍数
	@param multilooked_phase             多视相位
	*/
	int multilook(const ComplexMat& master, const ComplexMat& slave, int multilook_rg, int multilook_az, Mat& phase);
	/** @brief 将相位转换成cos和sin（实部和虚部）
	
	@param phase                     输入相位
	@param cos                       实部
	@param sin                       虚部
	@return 成功返回0，否则返回-1
	*/
	int phase2cos(const Mat& phase, Mat& cos, Mat& sin);
	/*84坐标系转经纬高坐标系
	* 参数1 84坐标系坐标
	* 参数2 经纬高坐标系坐标（度/度/米）
	*/
	int xyz2ell(Mat xyz, Mat& llh);
	/*经纬高坐标系转84坐标系
	* 参数1 经纬高坐标系坐标（纬度/经度/高度）
	* 参数2 84坐标系坐标
	*/
	int ell2xyz(Mat llh, Mat& xyz);


	/*
	* SAR图像文件读写操作工具
	*/

	/*从img文件中读出图像数据并转换成opencv mat
	* 参数1 img文件名
	* 参数2 cvmat
	*/
	int img2cvmat(const char* filename, Mat& cvmat);
	/*将cvmat数据写入到img文件中
	* 参数1 img文件名
	* 参数2 cvmat
	*/
	int cvmat2img(
		const char* filename, 
		Mat& cvmat,
		Esip_ImgLayer* imglayer,
		Esip_LayerInfo* layerinfo
	);
	/*将ComplexMat数据保存到img文件中(将SLC写入到img文件)
	* 参数1 目标img文件
	* 参数2 ComplexMat数据
	* 参数3 Esip_ImgLayer
	* 参数4 Esip_LayerInfo
	*/
	int complexmat2img(const char* filename, ComplexMat& complexmat, Esip_ImgLayer* imglayer,  Esip_LayerInfo* layerinfo);
	/*从img文件中读出复图像数据并转换成ComplexMat
	* 参数1 img文件名
	* 参数2 ComplexMat
	*/
	int img2complexmat(const char* filename, ComplexMat& complexmat);
	

	/*从img图像文件中读取复图像实部和虚部数据
	* 参数1 文件名
	* 参数2 图像宽度（列数）
	* 参数3 图像高度（行数）
	* 参数4 实部数据（返回值）
	* 参数5 虚部数据（返回值）
	*/
	int getRealImag(const char* filename, int width, int height, char* realBuf, char* imagBuf);
	/*从img文件中读取图像数据（干涉相位/差分干涉相位）
	* 参数1 源img文件
	* 参数2 读出数据块指针
	*/
	int getImagesData(const char* filename, char* outbuffer);
	/*获取img文件头部Flag（注意：此函数尽量不要单独使用，因为它没有关闭打开的文件）
	* 参数1 文件名
	* 参数2 Esip_HeaderTag结构体
	* 返回值 偏移后的文件指针
	*/
	FILE* getFileHeadFlag(const char* filename, Esip_HeaderTag* esip_header);
	/*获取根节点地址（注意：此函数尽量不要单独使用，因为它没有关闭打开的文件）
	* 参数1 文件指针
	* 参数2 文件头
	* 返回值 根节点地址
	*/
	GUInt32 getRootNodeptr(FILE* fp, GInt32 headerPtr);
	/*获取图层信息节点地址
	* 参数1 文件指针
	* 参数2 根节点地址
	* 返回值 图层信息节点地址
	*/
	GUInt32 getImgLayerNodeptr(FILE* fp, GInt32 rootPtr);
	/*获取图层信息子节点指针
	* 参数1 文件指针
	* 参数2 图层信息根节点地址
	* 参数3 图层信息地址
	* 参数4 图像块地址
	*/
	int getImgLayerChildNodeptr(FILE* fp, GInt32 imgLayerPtr, GInt32* layerInfoPtr, GInt32* imgBlockPtr);
	/*读取图像块数据
	* 参数1 文件指针
	* 参数2 图像块地址
	* 参数3 数据buffer
	*/
	int getImgBlockData(FILE* fp, GInt32 imgBlockPtr, char* dataBuf);
	/*获取img文件中图像的宽度/高度/图像类型/单像素占内存大小
	* 参数1 文件名
	* 参数2 宽度（返回值）
	* 参数3 高度（返回值）
	* 参数4 图像类型（SLC/phase/diff_phase）
	* 参数5 每个像素占用内存大小（byte）
	*/
	int getImgWidthAndHeight(const char* filename, GInt32* width, GInt32* height, Esip_LayerType* imgType, int* bytes_perpixel);
	/*从img文件中获取imglayer信息
	* 参数1 img文件路径
	* 参数2 imglayer结构体指针（返回值）
	*/
	int get_imglayer(const char* filename, Esip_ImgLayer* imglayer);
	/*从img文件中获取layerinfo信息
	* 参数1 img文件路径
	* 参数2 layerinfo结构体指针
	*/
	int get_layerinfo(const char* filename, Esip_LayerInfo* layerinfo);
	/*读取img文件中的控制点和轨道数据等辅助参数(若返回值小于0，则检查是否需要释放相应内存)
	* 参数1 img文件路径
	* 参数2 读出数据指针（返回值）
	*/
	int getAdditionalData(const char* filename, Esip_Additional_Data* esip_additional_data);
	/*向img文件写入控制点和轨道数据等辅助参数
	* 参数1 img文件路径
	* 参数2 写入数据指针
	*/
	int writeAdditionalData(const char* filename, Esip_Additional_Data* esip_additional_data);
	/*读取img文件中的控制点和轨道数据等辅助参数的地址
	* 
	*/
	int getAddtionalAddr(FILE* fp, GInt32 imgLayerPtr, long* AdditionalDataPtr);
	/*释放Esip_Additional_Data中的malloc分配的内存
	* 参数1 Esip_Additional_Data指针
	*/
	int free_AdditionalData(Esip_Additional_Data* ptr);
	/*初始化Esip_Additional_Data结构体里面的数据指针
	* 参数1 Esip_Additional_Data指针
	*/
	int init_AdditionalDataPtr(Esip_Additional_Data* ptr);
	/*初始化并填充Esip_Additional_Data结构体数据内容（从文本文件中读取）
	* 参数1 Esip_Additional_Data结构体
	* 参数2 目标文本文件名
	*/
	int init_AdditionalData(Esip_Additional_Data* addtional_data, const char* filename);
	/*获取当前时间
	* 参数1 时间结构体指针
	*/
	void gettimeofday(ModTime* tp);
	/*构造并写入EsipHeaderTag
	* 参数1 文件指针
	* 参数2 文件头地址
	*/
	int createEsipHeaderTag(FILE* fp, GInt32 headerPrt);
	/*写入文件头信息
	*/
	int createEsipFile(FILE* fp, int file_pos, Esip_File* esip_file);
	/*写入根节点*/
	int writeRootNode(FILE* fp, int nodePos, int* childPtr);
	/**/
	int createEsipEntry(FILE* fp, int nodePos, Esip_Entry* esip_entry);
	/*写入imglayer节点*/
	int writeImgLayerNode(FILE* fp, int nodePos, Esip_Entry* esip_entry, Esip_ImgLayer* esip_imglayer);
	/*写入imglayer数据*/
	int writeImgLayerData(FILE* fp, int nodePos, Esip_ImgLayer* esip_imgLayer);
	/**/
	int writeHfaNodeData(FILE* fp, int nodePos, void* buf, int buflen);
	/*向img文件写入实部和虚部数据（SLC）
	* 参数1 文件指针
	* 参数2 数据节点地址
	* 参数3 实部数据
	* 参数4 虚部数据
	* 参数5 SLC图像宽度
	* 参数6 SLC图像高度
	*/
	int writeRealAndImagData(FILE* fp, int nodePos, void* realBuf, void* imagBuf, int width, int height);
	/*构造Esip_ImgLayer结构体
	* 参数1 Esip_ImgLayer结构体指针
	* 参数2 图像宽度
	* 参数3 图像高度
	* 参数4 图像类型（单视复图像/干涉相位图/差分干涉相位图）
	* 参数5 每个像素所占字节大小
	*/
	int init_esip_imglayer(
		Esip_ImgLayer* imglayer,
		int width,
		int height,
		int type,
		int bytes_perpixel
	);
	/*构造Esip_LayerInfo结构体
	* 参数1 Esip_LayerInfo结构体指针
	* 参数2 数据获取方式（1：仿真数据，2：真实数据）
	* 参数3 方位向采样间隔
	* 参数4 距离向采样间隔
	* 参数5 方位向起始时间
	* 参数6 距离向起始时间
	*/
	int init_esip_layerinfo(
		Esip_LayerInfo* layerinfo,
		int data_simulation,
		double azimuthSample,
		double rangeSample,
		double azimuthStartTime,
		double rangeStartTime
	);


	/*******************************************************/
	/*                     图像存储工具集                  */
	/*******************************************************/

	/*量化保存SLC功率图
	* 参数1 目标文件名
	* 参数2 功率量化参数（可视范围dB）
	* 参数3 单视复图像
	*/
	int saveSLC(const char* filename, double db, ComplexMat& SLC);
	/*保存干涉相位图
	* 参数1 目标文件名
	* 参数2 颜色映射（jet/hsv/cool/parula等）
	* 参数3 待保存相位
	*/
	int savephase(const char* filename, const char* colormap, Mat phase);
	/*图像重采样
	* 参数1 原图像
	* 参数2 目标图像
	* 参数3 目标图像高度
	* 参数4 目标图像宽度
	*/
	int resampling(const char* Src_file, const char* Dst_file, int dst_height, int dst_width);
	/*量化SAR图像与干涉相位叠加
	* 参数1 量化SAR图像
	* 参数2 干涉相位图
	* 参数3 叠加图像
	* 参数4 SAR图像占比
	*/
	int amplitude_phase_blend(
		const char* amplitude_file,
		const char* phase_file,
		const char* blended_file,
		double SAR_ratio = 0.9
	);



	/*******************************************************/
	/*                Delaunay三角网相关函数库             */
	/*******************************************************/

	/*从.edge文件读取Delaunay三角网的边信息
	* 参数1 .edge文件
	* 参数2 指向边结构体的指针（返回值，内存需要手动释放）
	* 参数3 指向边个数的指针（返回值）
	* 参数4 统计每个节点的邻接边数（返回值，内存需要手动释放）
	* 参数5 节点数
	*/
	int read_edges(const char* filename, tri_edge** edges, long* num_edges, int** neighbours, long num_nodes);
	/** @brief 从.edge文件读取Delaunay三角网的边信息
	
	@param edge_file               .edge文件
	@param num_nodes               节点数
	@param edges                   Delaunay三角网边数组（返回值）
	@param node_neighbours         每个节点的邻接边数（返回值）
	@return  成功返回0， 否则返回-1
	*/
	int read_edges(
		const char* edge_file,
		vector<tri_edge>& edges,
		vector<int>& node_neighbours,
		long num_nodes
	);
	/*初始化Delaunay三角网节点
	* 参数1 节点数组（返回值）
	* 参数2 相位(double型)
	* 参数3 相位mask（int 型）
	* 参数4 edges结构体数组
	* 参数5 edges个数
	* 参数6 每个节点的邻接边信息
	* 参数7 节点数
	*/
	int init_tri_node(
		vector<tri_node>& node_array,
		Mat& phase,
		Mat& mask,
		tri_edge* edges,
		long num_edges,
		int* num_neighbour,
		int num_nodes
	);
	/** @brief 初始化Delaunay三角网节点
	
	@param node_array                 节点数组（返回值）
	@param phase                      相位值
	@param mask                       相位掩膜
	@param edges                      Delaunay三角网络边结构体数组
	@param node_neighbours            每个节点的邻边个数
	@param num_nodes                  节点数
	@return 成功返回0，否则返回-1
	*/
	int init_tri_node(
		vector<tri_node>& node_array,
		const Mat& phase,
		const Mat& mask,
		const vector<tri_edge>& edges,
		const vector<int>& node_neighbours,
		int num_nodes
	);
	/** @brief 初始化Delaunay三角网络边相位差
	
	@param edges                  Delaunay三角网络边数组（已经使用read_edges函数初始化过的）
	@param node_array             Delaunay三角网络节点数组（已经使用init_tri_node函数初始化过的）
	@return 成功返回0，否则返回-1
	*/
	int init_edge_phase_diff(
		vector<tri_edge>& edges,
		const vector<tri_node>& node_array
	);
	/*初始化Delaunay三角网边的相位质量
	* 参数1 相位质量图
	* 参数2 Delaunay三角网边结构体数组指针
	* 参数3 Delaunay三角网边结构体数组大小
	* 参数4 Delaunay三角网节点数组
	*/
	int init_edges_quality(Mat& quality, tri_edge* edges, int num_edges, vector<tri_node>& nodes);
	/** @brief 初始化Delaunay三角网边的相位质量指数
	
	@param quality_index                  相位质量图指数（与相位质量相反）
	@param edges                          Delaunay三角网边结构体数组
	@param nodes                          Delaunay三角网节点数组
	@return 成功返回0， 否则返回-1
	*/
	int init_edges_quality(
		const Mat& quality_index,
		vector<tri_edge>& edges,
		const vector<tri_node>& nodes
	);
	/*从.ele文件和.neigh文件读取Delaunay三角网的三角形信息
	* 参数1 .ele文件
	* 参数2 .neigh文件
	* 参数3 三角形结构体数组指针（返回值, 内存需要手动释放）
	* 参数4 三角形个数（返回值）
	* 参数5 Delaunay三角网节点数组
	* 参数6 Delaunay三角网边数组
	* 参数7 Delaunay三角网边数量
	*/
	int read_triangle(
		const char* ele_file,
		const char* neigh_file,
		triangle** tri,
		int* num_triangle,
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edgs
	);
	/** @brief 从.ele文件和.neigh文件读取Delaunay三角网的三角形信息
	
	@param ele_file                        .ele文件
	@param neigh_file                      .neigh文件
	@param triangle                        三角形结构体数组（返回值）
	@param nodes                           Delaunay三角网节点数组
	@param edges                           Delaunay三角网边数组
	@return 成功返回0，否则返回-1
	*/
	int read_triangle(
		const char* ele_file,
		const char* neigh_file,
		vector<triangle>& triangle,
		vector<tri_node>& nodes,
		vector<tri_edge>& edges
	);
	/*生成Delaunay三角网
	* 参数1 .node文件
	* 参数2 triangle.exe程序路径
	*/
	int gen_delaunay(const char* filename, const char* exe_path);
	/*写.node文件
	* 参数1 .node文件
	* 参数2 节点数组
	*/
	int write_node_file(const char* filename, const Mat& mask);




	/*********************************************************/
    /*                PS-InSAR 常用函数                      */
    /*********************************************************/

	/*振幅离差指数法筛选PS点（D_A）
	* 参数1 SAR幅度矩阵组
	* 参数2 振幅离差阈值
	* 参数3 mask（满足条件的PS点位置mask为1，其他为0）
	*/
	int PS_amp_dispersion(const vector<Mat>& amplitude, double thresh, Mat& mask);
	/*fifth-order butterworth filter（五阶巴特沃斯滤波器）
	* 参数1 grid_size
	* 参数2 n_win
	* 参数3 low_pass_wavelength
	* 参数4 滤波器系数（返回值）
	*/
	int butter_lowpass(int grid_size, int n_win, double low_pass_wavelength, Mat& lowpass);
	/*circle_shift
	*/
	int circshift(Mat& out, const cv::Point& delta);
	/*fftshift2
	*/
	int fftshift2(Mat& out);
	/*ifftshift
	*/
	int ifftshift(Mat& out);
	/*二维傅里叶变换
	* 参数1 输入矩阵
	* 参数2 输出结果
	*/
	int fft2(Mat& Src, Mat& Dst);
	/*复数二维傅里叶变换
	* 参数1 输入矩阵
	* 参数2 输出结果
	*/
	int fft2(ComplexMat& src, ComplexMat& dst);
	/*逆二维傅里叶变换
	* 参数1 输入矩阵
	* 参数2 输出结果
	*/
	int ifft2(ComplexMat& src, ComplexMat& dst);
	/*求标准差
	* 参数1 输入矩阵
	* 参数2 标准差返回值
	*/
	int std(const Mat& input, double* std);
	/*裁剪感兴趣的SAR图像区域（AOI）
	* 参数1 SAR图像序列文件名（img格式）
	* 参数2 裁剪后的SAR图像序列保存路径
	* 参数2 AOI中心纬度/经度/高度(1×3)
	* 参数3 AOI宽（列数）
	* 参数4 AOI高（行数）
	*/
	int PS_cut_AOI(
		vector<string>& SAR_images_files,
		const char* save_path,
		Mat& llh,
		int rows,
		int cols
	);
	/*SAR图像干涉相位序列去平地
	* 参数1 干涉相位序列（原地操作）
	* 参数2 干涉组合
	* 参数3 卫星轨道参数（插值后， n_images×6）
	* 参数4 地面控制点信息（纬度/经度/高度/行/列, n_gcps × 5）
	* 参数5 SAR图像左上角在原SAR图像中的行数(可以直接设置为1)
	* 参数6 SAR图像左上角在原SAR图像中的列数(可以直接设置为1)
	* 参数7 收发方式（1单发单收， 2单发双收）
	* 参数8 波长
	*/
	int PS_deflat(
		vector<Mat>& interf_phase,
		Mat& interf_combination,
		vector<Mat>& pos,
		vector<Mat>& gcps,
		Mat& start_row,
		Mat& start_col,
		int mode,
		double lambda
	);
	/*去平地(线性拟合法)
	* 参数1 干涉相位（原地操作）
	* 参数2 主星轨道参数（插值后， n_images×6）
	* 参数3 辅星轨道参数（插值后， n_images×6）
	* 参数4 地面控制点信息（纬度/经度/高度/行/列, n_gcps × 5）
	* 参数5 SAR图像左上角在原SAR图像中的行数(可以直接设置为1)
	* 参数6 SAR图像左上角在原SAR图像中的列数(可以直接设置为1)
	* 参数7 收发方式（1单发单收， 2单发双收）
	* 参数8 波长
	*/
	int _PS_deflat(
		Mat& phase,
		Mat& pos1,
		Mat& pos2,
		Mat& gcps,
		int start_row,
		int start_col,
		int mode,
		double lambda
	);
	/** @brief 时序SAR图像联合配准(所有slc同时载入内存)
	
	@param SAR_images            时序SAR图像（inplace，原地操作）
	@param offset                配准后左上角偏移量(尺寸：n_images × 2) 
	@param Master_index          主图像序号(序号从1开始)
	@param coh_method            采用实相关还是复相关（0代表实相关， 1代表复相关）
	@param interp_times          插值倍数（2的n次幂）
	@param blocksize             子块大小（2的n次幂）
	*/
	int stack_coregistration(
		vector<ComplexMat>& SAR_images,
		Mat& offset,
		int Master_index,
		int coh_method,
		int interp_times,
		int blocksize
	);
	/** @brief 时序SAR图像联合配准(slc串行载入内存，以节省内存)
	
	@param SAR_images            时序SAR图像文件
	@param SAR_images_out        配准结果文件
	@param offset                配准后左上角偏移量(尺寸：n_images × 2)
	@param Master_index          主图像序号(序号从1开始)
	@param interp_times          插值倍数（2的n次幂）
	@param blocksize             子块大小（2的n次幂）
	*/
	int stack_coregistration(
		vector<string>& SAR_images,
		vector<string>& SAR_images_out,
		Mat& offset,
		int Master_index,
		int interp_times,
		int blocksize
	);
	/*生成干涉图组合
	* 参数1 SAR图像序列
	* 参数2 干涉相位序列（返回值）
	* 参数3 时间基线组合（返回值, 原地操作）
	* 参数4 空间基线自合（返回值, 原地操作）
	* 参数5 时间基线阈值（年）
	* 参数6 空间基线阈值（米）
	* 参数7 主图像序号(如果小于等于0则按照自由组合方式干涉)
	* 参数8 干涉组合
	*/
	int PS_gen_interferogram(
		vector<ComplexMat>& SAR_images,
		vector<Mat>& interf_phase,
		Mat& time_baseline,
		Mat& spatial_baseline,
		double time_thresh,
		double spatial_thresh,
		int Mast_index, 
		Mat& interf_combination
	);
	/*空间基线估计
	* 参数1 卫星轨道参数（插值后位置和速度， 尺寸：n*6）
	* 参数2 地面控制点信息（纬经高，尺寸：n_gcps * 3）
	* 参数3 波长
	* 参数4 基线长度（返回值）
	*/
	int PS_spatialbaseline_est(
		vector<Mat>& pos,
		Mat gcps,
		double lambda,
		Mat& MB_effect
	);
	/*Combined Low-pass Adaptive Phase filtering（StaMPS自适应滤波器）
	* 参数1 差分相位（复数形式）
	* 参数2 滤波后相位（复数形式）
	* 参数3 滤波器参数（一般为1）
	* 参数4 滤波器参数（一般为0.3）
	* 参数5 窗口大小（必须为偶数）
	* 参数6 补零长度（必须为偶数）
	* 参数7 低通滤波器系数(尺寸 = 窗口大小 + 补零长度)
	*/
	int Clap(
		ComplexMat& phase,
		ComplexMat& phase_filter,
		double alpha,
		double beta,
		int n_win,
		int n_pad,
		Mat& lowpass
	);
	/*topofit
	* 参数1 差分相位
	* 参数2 垂直基线长度
	* 参数3 num_trial_wraps
	* 参数4 K0（返回值）
	* 参数5 C0（返回值）
	* 参数6 coh0（返回值）
	* 参数7 残余相位（返回值）
	*/
	int PS_topofit(
		Mat& dif_phase, Mat& bperp,
		double num_trial_wraps,
		double* K0, double* C0,
		double* coh0,
		Mat& phase_residue
	);
	/*hist函数（统计直方图函数）
	* 参数1 待统计数据
	* 参数2 统计标准下限区间中心
	* 参数3 统计标准上限区间中心
	* 参数4 区间半径（n * 区间半径 = （统计标准上限区间中心 - 统计标准下限区间中心））
	* 参数5 统计输出
	*/
	int hist(
		Mat& input,
		double lowercenter,
		double uppercenter,
		double interval,
		Mat& out
	);
	/*经纬高转换到局部坐标系
	* 参数1 经纬高坐标
	* 参数2 局部中心经度
	* 参数3 局部中心纬度
	* 参数4 局部坐标（返回值）
	*/
	int llh2local(const Mat& llh, double lon0, double lat0, Mat& xy);
	/*将候选PS点重采样到规则网格上
	* 参数1 局部坐标系坐标
	* 参数2 网格间距大小（米）
	* 参数3 候选PS点在网格中的坐标（返回值）
	* 参数4 重采样网格点行数
	* 参数5 重采样网格点列数
	*/
	int gen_grid(
		const Mat& xy,
		double grid_size,
		Mat& grid_ij,
		int* rows,
		int* cols
	);
	/*计算随机相位的时间相干系数概率分布函数
	* 参数1 随机数个数（建议值30000）
	* 参数2 差分干涉图幅数
	* 参数3 垂直基线
	* 参数4 搜索区间半径
	* 参数5 概率分布（返回值）
	* 参数6 概率密度不为零的最小横坐标值
	*/
	int gen_CohDistributionOfRandomPhase(
		int nrand,
		int n_ifg,
		const Mat& bperp,
		double num_trial_wraps,
		Mat& Nr,
		int* Nr_max_nz_ix
	);
	/*第一次估计PS点的时间相干系数
	* 参数1/参数2     差分相位（行：PS候选点序号，列：干涉图序号） / 振幅离差指数(按列排)
	* 参数3/参数4     PS点在重采样网格中的坐标 / 垂直基线（与差分相位同尺寸）
	* 参数5/参数6     低通滤波器系数 / 随机相位时间相干系数分布（按行排100行1列）
	* 参数7/参数8     自适应滤波器参数（alpha建议值1.0， beta建议值0.3）
	* 参数9/参数10    优化搜索区间半径 / 迭代最大次数
	* 参数11/参数12   Nr参数 / 迭代收敛判断阈值
	* 参数13/参数14   K_ps / C_ps（返回值）
	* 参数15/参数16   coh_ps（时间相干系数返回值）/ fft窗口大小
	*/
	int PS_est_gamma_quick(
		const Mat& dif_phase,
		const Mat& D_A,
		const Mat& grid_ij,
		const Mat& bperp_mat,
		Mat& lowpass,
		Mat& Nr,
		double alpha,
		double beta,
		double n_trial_wraps,
		int gamma_max_iterations,
		int Nr_max_nz_ix,
		double gamma_change_convergence,
		Mat& K_ps,
		Mat& C_ps,
		Mat& coh_ps,
		int n_win
	);
	/*线性拟合时间相干系数阈值与D_A，并剔除小于阈值的PS候选点
	* 参数1 D_A（振幅离差, 按列排）
	* 参数2 PS候选点时间相干系数（按列排）
	* 参数3 Nr（按行排）
	* 参数4 拟合阈值（返回值）
	* 参数5 low_coh_thresh
	* 参数6 虚警率
	*/
	int coh_thresh_fit(
		const Mat& D_A,
		Mat& coh_ps,
		Mat& Nr,
		Mat& coh_thresh_D_A,
		int low_coh_thresh,
		double false_alarm
	);
	/*二维网格搜索相邻PS点差分相位模型相干系数最大值点，估计出delta_vel和delta_epsilon
	* 搜索区间间隔(1mm/y, 2m)
	* 
	* 参数1 每个PS_pair的相位序列（按行排开）
	* 参数2 与时间基线序列相关的参数（4 * pi * T_i / lambda， 与参数1行列相同）
	* 参数3 与空间基线序列相关的参数（4 * pi * berp_i / lambda / R_i / sin_theta_i， 与参数1行列相同）
	* 参数4 形变速率差搜索区间半径(mm/y)
	* 参数5 高程差搜索区间半径(m)
	* 参数6 MC值（返回值）
	* 参数7 形变速率差（返回值）
	* 参数8 高程差（返回值）
	*/
	int PS_topovel_fit_search(
		ComplexMat& ph,
		Mat& coef_delta_vel,
		Mat& coef_delta_height,
		double radius_delta_vel,
		double radius_delta_height,
		double* MC,
		double* delta_vel,
		double* delta_height
	);
	/*二维网格搜索(精搜索)相邻PS点差分相位模型相干系数最大值点，估计出delta_vel和delta_epsilon
	* 搜索区间间隔(0.1mm/y, 0.2m)
	* 
	* 参数1 每个PS_pair的相位序列（按行排开）
	* 参数2 与时间基线序列相关的参数（4 * pi * T_i / lambda， 与参数1行列相同）
	* 参数3 与空间基线序列相关的参数（4 * pi * berp_i / lambda / R_i / sin_theta_i， 与参数1行列相同）
	* 参数4 粗搜索结果（delta_vel）
	* 参数5 粗搜索结果（delta_height）
	* 参数6 形变速率差搜索区间半径(mm/y)
	* 参数7 高程差搜索区间半径(m)
	* 参数8 MC值（返回值）
	* 参数9 形变速率差（返回值）
	* 参数10 高程差（返回值）
	*/
	int PS_topovel_fit_search_2(
		ComplexMat& ph,
		Mat& coef_delta_vel,
		Mat& coef_delta_height,
		double rough_delta_vel_center,
		double rough_delta_height_center,
		double radius_delta_vel,
		double radius_delta_height,
		double* MC,
		double* delta_vel,
		double* delta_height
	);
	/*通过模型相干系数（Model Coherence）对PS边进行筛选，小于阈值被丢弃
	* 参数1 PS节点数组
	* 参数2 PS边结构体数组
	* 参数3 PS边结构体数
	* 参数4 模型相干系数阈值
	*/
	int Dump_unqualified_pairnodes(
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges,
		double MC_thresh
	);
	/*增量积分集成（积分获取每个PS点的形变速率和高程误差）
	* 参数1 PS点数组
	* 参数2 PS边结构体数组
	* 参数3 PS边数量
	* 参数4 PS积分边长阈值（超过此阈值不通过此边积分集成）
	* 参数5 MC阈值（超过此阈值不通过此边积分集成）
	*/
	int PS_topovel_incre(
		vector<tri_node>& nodes,
		tri_edge* edges,
		int num_edges, 
		double distance_thresh,
		double MC_thresh
	);


	/*
	* 卫星轨道插值
	* 参数1：卫星轨道参数（未插值）
	* 参数2：插值时间间隔（s）
	* 参数3：插值结果
	*/
	int stateVec_interp(Mat& stateVec, double time_interval, Mat& stateVec_interp);
	/*
	* 从h5（SLC）图像文件中裁剪AOI
	* 参数1：h5文件
	* 参数2：AOI左上角经度
	* 参数3：AOI左上角纬度
	* 参数4：AOI右下角经度
	* 参数5：AOI右下角纬度
	* 参数6：输出单视复图像
	* 参数7：AOI左上角行偏移量（从0开始，0代表不偏移）
	* 参数8：AOI左上角列偏移量（从0开始，0代表不偏移）
	*/
	int get_AOI_from_h5SLC(
		const char* h5_file,
		double lon_topleft,
		double lat_topleft,
		double lon_bottomright,
		double lat_bottomright,
		ComplexMat& slc,
		int* offset_row = NULL,
		int* offset_col = NULL
	);
	/** @brief 从h5（SLC）文件中裁剪出AOI区域
	
	@param h5_file h5文件
	@param lon_center AOI中心经度
	@param lat_center AOI中心纬度
	@param width AOI宽度（m）
	@param height AOI高度（m）
	@param slc 裁剪结果slc
	@param offset_row AOI左上角在原图像中行偏移量（从0开始，0代表不偏移）
	@param offset_col AOI左上角在原图像中行偏移量（从0开始，0代表不偏移）
	*/
	int get_AOI_from_h5slc(
		const char* h5_file,
		double lon_center,
		double lat_center,
		double width,
		double height,
		ComplexMat& slc,
		int* offset_row = NULL,
		int* offset_col = NULL
	);
	/** @brief 从h5（SLC）文件中获取裁剪AOI区域的尺寸

	@param h5_file          h5文件
	@param lon_center       AOI中心经度
	@param lat_center       AOI中心纬度
	@param width            AOI宽度（m）
	@param height           AOI高度（m）
	@param AOI_rows         AOI行数
	@param AOI_cols         AOI列数
	@param offset_row       AOI左上角在原图像中行偏移量（从0开始，0代表不偏移）
	@param offset_col       AOI左上角在原图像中行偏移量（从0开始，0代表不偏移）
	*/
	int get_AOI_size(
		const char* h5_file,
		double lon_center,
		double lat_center,
		double width,
		double height,
		int* AOI_rows,
		int* AOI_cols,
		int* offset_row = NULL,
		int* offset_col = NULL
	);
	/** @brief 坐标转换工具函数
	
	@param coefficient       转换系数矩阵
	@param coord_in_1        原坐标矩阵1(1和2的顺序很重要，经度/行坐标在前)
	@param coord_in_2        原坐标矩阵2
	@param coord_out         转换结果矩阵
	*/
	int coord_conversion(
		Mat& coefficient,
		Mat& coord_in_1,
		Mat& coord_in_2,
		Mat& coord_out
	);
	/** @brief 基线估计

	@param stateVec1               主星轨道（未插值）
	@param stateVec2               辅星轨道（未插值）
	@param lon_coef                主星坐标转换系数（图像坐标-->经度）
	@param lat_coef                主星坐标转换系数（图像坐标-->纬度）
	@param offset_row              主图像左上角在原始图像中的行偏移量
	@param offset_col              主图像左上角在原始图像中的列偏移量
	@param scene_height            场景高度(像素行数)
	@param scene_width             场景宽度(像素列数)
	@param interp_interval1        主星轨道插值时间间隔（1/prf）
	@param interp_interval2        辅星轨道插值时间间隔（1/prf）
	@param B_effect                垂直基线长度（返回值）
	@param B_parallel              平行基线长度（返回值）
	@param sigma_B_effect          垂直基线估计标准差（返回值）
	@param sigma_B_parallel        平行基线估计标准差（返回值）
	*/
	int baseline_estimation(
		const Mat& stateVec1,
		const Mat& stateVec2,
		const Mat& lon_coef,
		const Mat& lat_coef,
		int offset_row,
		int offset_col,
		int scene_height,
		int scene_width,
		double interp_interval1,
		double interp_interval2,
		double* B_effect,
		double* B_parallel,
		double* sigma_B_effect = NULL,
		double* sigma_B_parallel = NULL
	);
	/** @brief 统计同质检验
	
	@param pixel1            待检验像元1幅度序列(size: n_images×1)
	@param pixel2            待检验像元2幅度序列(size: n_images×1)
	@param homo_flag         是否为同质像元(返回0则为同质像元，-1则为非同质像元)
	@param alpha             显著性水平（可以设定的值为 0.20,0.15,0.10,0.05,0.025,0.01,0.005,0.001。默认为0.05）
	@param method            检验方法（"KS":Kolmogorov-Smirnov检验，"AD":Anderson-Darling检验, 默认为KS检验）
	@return                  正常运行返回0，报错返回-1
	*/
	int homogeneous_test(
		const Mat& pixel1,
		const Mat& pixel2,
		int* homo_flag,
		double alpha = 0.05,
		const char* method = "KS"
	);
	/** @brief Hermitian矩阵特征值分解
	
	@param input               输入复矩阵（n×n, double型）
	@param eigenvalue          特征值（n×1实矩阵,从大到小排列）
	@param eigenvector         特征向量（n×n复矩阵， 列向量为特征向量）
	@return                    成功返回0，否则返回-1              
	*/
	int HermitianEVD(
		const ComplexMat& input,
		Mat& eigenvalue,
		ComplexMat& eigenvector
	);
	/** @brief 时序SAR图像复相关矩阵估计
	
	@param slc_series               slc数据堆栈
	@param coherence_matrix         相关矩阵（复数, 返回值）
	@param est_window_width         估计窗口宽度（奇数）
	@param est_window_height        估计窗口高度（奇数）
	@param ref_row                  （若进行统计同质检验）参考点行坐标，不进行同质检验则不需要此参数
	@param ref_col                  （若进行统计同质检验）参考点列坐标，不进行同质检验则不需要此参数
	@param b_homogeneous_test       是否进行统计同质检验（同质检验参考像素默认为中间点像素）
	@param b_normalize              估计相关矩阵时slc序列是否归一化处理
	@return                         成功返回0，否则返回-1 
	*/
	int coherence_matrix_estimation(
		const vector<ComplexMat>& slc_series,
		ComplexMat& coherence_matrix,
		int est_window_width,
		int est_window_height,
		int ref_row,
		int ref_col,
		bool b_homogeneous_test = true,
		bool b_normalize = true
	);
	/** @brief 多基线时间序列相位估计（分块读取、计算、储存）
	
	@param coregis_slc_files              配准后SAR图像数据堆栈（文件）
	@param phase_files                    时间序列干涉相位（文件，与coregis_slc_files数量相同，主图像相位为0）
	@param coherence_files                各辅图像与主图像之间的相关系数文件（是否估计相关系数取决于输入参数b_coh_est）
	@param master_indx                    主图像序号（从1开始）
	@param blocksize_row                  子块尺寸（行，必须大于同质检验搜索窗口半径）
	@param blocksize_col                  子块尺寸（列，必须大于同质检验搜索窗口半径）
	@param out_mask                       掩膜输出（标记经过EVD法估计的像素点，与参数thresh_c1_to_c2有关）
	@param b_coh_est                      是否估计相关系数（默认是）
	@param homogeneous_test_wnd           同质检验搜索窗口大小（奇数，homogeneous_test_wnd×homogeneous_test_wnd， 默认为21×21）
	@param thresh_c1_to_c2                协方差矩阵第2特征值与第1特征值比值阈值（0-1之间，默认为0.7）
	@param b_flat                         是否去平地相位（默认是）
	@param b_normalize                    协方差矩阵是否归一化（默认是）
	*/
	int MB_phase_estimation(
		vector<string> coregis_slc_files,
		vector<string> phase_files,
		vector<string> coherence_files,
		int master_indx,
		int blocksize_row,
		int blocksize_col,
		Mat& out_mask,
		bool b_coh_est = true,
		int homogeneous_test_wnd = 21,
		double thresh_c1_to_c2 = 0.7,
		bool b_flat = true,
		bool b_normalize = true
	);
	/** @brief 区域生长法解缠（delaunay三角网）
	
	@param nodes                       Delaunay三角网络节点数组
	@param edges                       Delaunay三角网络边结构体数组
	@param start_edge                  积分起始边序号（从1开始）
	@param distance_thresh             边长阈值，超过此阈值不通过此边积分
	@param quality_thresh              质量阈值，低于此阈值不通过此边积分
	@return 成功返回0，否则返回-1
	*/
	int unwrap_region_growing(
		vector<tri_node>& nodes,
		const vector<tri_edge>& edges,
		size_t start_edge,
		double distance_thresh,
		double quality_thresh
	);
	/** @brief 3D相位解缠，1D时间-->2D空间
	
	@param mask                          需要解缠的像素点掩膜(int型)
	@param quality_map                   质量图
	@param wrapped_phase_series          待解缠相位时间序列
	@param unwrapped_phase_series        解缠相位时间序列
	@param delaunay_exe_path             生成Delaunay三角网的可执行程序（delaunay.exe）所在路径
	@param tmp_file_path                 产生的临时文件储存路径
	@param distance_thresh               积分解缠距离阈值，超过此阈值不通过该边积分解缠，不能小于1
	@param quality_thresh                质量阈值，低于此阈值不通过此边积分
	@return 成功返回0，否则返回-1
	*/
	int unwrap_3D(
		const Mat& mask,
		const vector<Mat>& quality_map,
		vector<Mat>& wrapped_phase_series,
		vector<Mat>& unwrapped_phase_series,
		const char* delaunay_exe_path,
		const char* tmp_file_path,
		double distance_thresh,
		double quality_thresh
	);
	/** @brief 3D相位解缠（Delaunay三角网最小费用流法）
	
	@param mask                          需要解缠的像素点掩膜(int型)
	@param quality_map                   质量图
	@param wrapped_phase_series          待解缠相位时间序列
	@param unwrapped_phase_series        解缠相位时间序列
	@param delaunay_exe_path             生成Delaunay三角网的可执行程序（delaunay.exe）所在路径
	@param mcf_exe_path                  最小费用流求解器可执行程序（mcf.exe）所在路径
	@param tmp_file_path                 产生的临时文件储存路径
	@param distance_thresh               积分解缠距离阈值，超过此阈值不通过该边积分解缠，不能小于1
	@return 成功返回0，否则返回-1
	*/
	int unwrap_3D_mcf(
		const Mat& mask,
		const vector<Mat>& quality_map,
		vector<Mat>& wrapped_phase_series,
		vector<Mat>& unwrapped_phase_series,
		const char* delaunay_exe_path,
		const char* mcf_exe_path,
		const char* tmp_file_path,
		double distance_thresh
	);
	/** @brief 自适应分块3D相位解缠
	
	@param mask                          需要解缠的像素点掩膜(int型)
	@param quality_map                   质量图
	@param wrapped_phase_series          待解缠相位时间序列
	@param unwrapped_phase_series        解缠相位时间序列
	@param delaunay_exe_path             生成Delaunay三角网的可执行程序（delaunay.exe）所在路径
	@param mcf_exe_path                  最小费用流求解器可执行程序（mcf.exe）所在路径
	@param tmp_file_path                 产生的临时文件储存路径
	@param distance_thresh               积分解缠距离阈值，超过此阈值不通过该边积分解缠，不能小于1
	@param quality_thresh                质量阈值，低于此阈值不通过此边积分
	@return 成功返回0，否则返回-1
	*/
	int unwrap_3D_adaptive_tiling(
		const Mat& mask,
		const vector<Mat>& quality_map,
		vector<Mat>& wrapped_phase_series,
		vector<Mat>& unwrapped_phase_series,
		const char* delaunay_exe_path,
		const char* mcf_exe_path,
		const char* tmp_file_path,
		double distance_thresh,
		double quality_thresh
	);
private:
	char error_head[256];
	char parallel_error_head[256];

};


#endif