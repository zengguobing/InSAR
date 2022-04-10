#include"SBAS.h"
#include"..\include\Utils.h"
#include"..\include\FormatConversion.h"
#include"..\include\ComplexMat.h"

#ifdef _DEBUG
#pragma comment(lib, "ComplexMat_d.lib")
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "FormatConversion_d.lib")
#else
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "FormatConversion.lib")
#endif // _DEBUG


inline bool return_check(int ret, const char* detail_info, const char* error_head)
{
	if (ret < 0)
	{
		fprintf(stderr, "%s %s\n\n", error_head, detail_info);
		return true;
	}
	else
	{
		return false;
	}
}

SBAS_node::SBAS_node()
{
	this->B_spatial = 0.0;
	this->B_temporal = 0.0;
	this->b_unwrapped = false;
	this->x = 0;
	this->y = 0;
	this->deformation_vel = 0.0;
	this->epsilon_height = 0.0;
	this->neigh_edges = NULL;
	this->num_neigh_edges = 0;
	this->phase = 0.0;
}

SBAS_node::SBAS_node(const SBAS_node& cp)
{
	this->B_spatial = cp.B_spatial;
	this->B_temporal = cp.B_temporal;
	this->b_unwrapped = cp.b_unwrapped;
	this->x = cp.x;
	this->y = cp.y;
	this->deformation_vel = cp.deformation_vel;
	this->epsilon_height = cp.epsilon_height;
	int num_node = cp.num_neigh_edges <= 0 ? 1 : cp.num_neigh_edges;
	if (cp.neigh_edges == NULL)
	{
		this->neigh_edges = NULL;
	}
	else if (cp.neigh_edges == this->neigh_edges)
	{

	}
	else
	{
		this->neigh_edges = (int*)malloc(sizeof(int) * num_node);
		if (this->neigh_edges != NULL)
		{
			std::memcpy(this->neigh_edges, cp.neigh_edges, sizeof(int) * num_node);
		}
	}
	this->num_neigh_edges = cp.num_neigh_edges;
	this->phase = cp.phase;
}

SBAS_node::SBAS_node(int num_neigh_edge)
{
	this->B_spatial = 0.0;
	this->B_temporal = 0.0;
	this->b_unwrapped = false;
	this->y = 0;
	this->x = 0;
	this->deformation_vel = 0.0;
	this->epsilon_height = 0.0;
	this->num_neigh_edges = num_neigh_edge;
	this->phase = 0.0;
	if (num_neigh_edge > 0)
	{
		this->neigh_edges = (int*)malloc(sizeof(int) * num_neigh_edge);
	}
	else
	{
		this->neigh_edges = NULL;
	}
	if (this->neigh_edges != NULL)
	{
		for (int i = 0; i < num_neigh_edge; i++)
		{
			*(this->neigh_edges + i) = -1;//初始化邻接边序号都为-1
		}
	}

}

SBAS_node::~SBAS_node()
{
	if (this->neigh_edges != NULL)
	{
		free(this->neigh_edges);
		this->neigh_edges = NULL;
	}
}
SBAS_node SBAS_node::operator=(const SBAS_node& src)
{
	if (src.neigh_edges == this->neigh_edges && this->neigh_edges != NULL)//两者相等
	{
		return *this;
	}
	else
	{
		if (this->neigh_edges)
		{
			free(this->neigh_edges);
			this->neigh_edges = NULL;
		}
		if (src.num_neigh_edges > 0)
		{
			this->neigh_edges = (int*)malloc(src.num_neigh_edges * sizeof(int));
			if (this->neigh_edges != NULL && src.neigh_edges != NULL)
			{
				memcpy(this->neigh_edges, src.neigh_edges, src.num_neigh_edges * sizeof(int));
			}
		}
		this->B_spatial = src.B_spatial;
		this->B_temporal = src.B_temporal;
		this->b_unwrapped = src.b_unwrapped;
		this->y = src.y;
		this->x = src.x;
		this->num_neigh_edges = src.num_neigh_edges;
		this->phase = src.phase;
		this->epsilon_height = src.epsilon_height;
		this->deformation_vel = src.deformation_vel;
		return *this;
	}
}

SBAS::SBAS()
{
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "SBAS_DLL_ERROR: error happens when using ");
}

SBAS::~SBAS()
{
}

int SBAS::write_spatialTemporal_node(
	const char* nodeFile, 
	Mat& B_temporal,
	Mat& B_effect
)
{
	if (!nodeFile ||
		B_temporal.rows != 1 ||
		B_temporal.cols < 2 ||
		B_effect.rows != 1 ||
		B_effect.cols != B_temporal.cols ||
		B_effect.type() != CV_64F ||
		B_temporal.type() != CV_64F
		)
	{
		fprintf(stderr, "write_spatialTemporal_node(): input check failed!\n");
		return -1;
	}
	double B_span, T_span, B_min, B_max, T_min, T_max;
	cv::minMaxLoc(B_temporal, &T_min, &T_max);
	cv::minMaxLoc(B_effect, &B_min, &B_max);
	B_span = B_max - B_min; 
	T_span = T_max - T_min;
	if (B_span < 1e-10 || T_span < 1e-10)
	{
		fprintf(stderr, "write_spatialTemporal_node(): B_span < 1e-10 || T_span < 1e-10!\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(nodeFile, "wt");
	if (fp == NULL)
	{
		fprintf(stderr, "write_spatialTemporal_node(): can't open %s\n", nodeFile);
		return -1;
	}
	int nonzero = B_temporal.cols;
	int dim, attr1, attr2, count;
	dim = 2;
	attr1 = 0;
	attr2 = 0;
	fprintf(fp, "%d %d %d %d\n", nonzero, dim, attr1, attr2);
	int cols = B_temporal.cols;
	count = 1;
	double x, y;
	for (int i = 0; i < cols; i++)
	{
		x = (B_temporal.at<double>(0, i) - T_min) / T_span;
		y = (B_effect.at<double>(0, i) - B_min) / B_span;
		fprintf(fp, "%d %lf %lf\n", count, x, y);
		count++;
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int SBAS::set_spatialTemporalBaseline(
	vector<SBAS_node>& nodes,
	Mat& B_temporal,
	Mat& B_effect
)
{
	int num_nodes = nodes.size();
	if (num_nodes != B_effect.cols ||
		B_temporal.rows != 1 ||
		B_temporal.cols < 2 ||
		B_effect.rows != 1 ||
		B_effect.cols != B_temporal.cols ||
		B_effect.type() != CV_64F ||
		B_temporal.type() != CV_64F
		)
	{
		fprintf(stderr, "set_spatialTemporalBaseline(): input check failed!\n");
		return -1;
	}
	double B_span, T_span, B_min, B_max, T_min, T_max;
	cv::minMaxLoc(B_temporal, &T_min, &T_max);
	cv::minMaxLoc(B_effect, &B_min, &B_max);
	B_span = B_max - B_min;
	T_span = T_max - T_min;
	if (B_span < 1e-10 || T_span < 1e-10)
	{
		fprintf(stderr, "set_spatialTemporalBaseline(): B_span < 1e-10 || T_span < 1e-10!\n");
		return -1;
	}
	for (int i = 0; i < num_nodes; i++)
	{
		nodes[i].B_spatial = B_effect.at<double>(0, i);
		nodes[i].B_temporal = B_temporal.at<double>(0, i);
		nodes[i].x = (B_temporal.at<double>(0, i) - T_min) / T_span;
		nodes[i].y = (B_effect.at<double>(0, i) - B_min) / B_span;
	}
	return 0;
}

int SBAS::read_edges(
	const char* edge_file,
	int num_nodes,
	vector<SBAS_edge>& edges,
	vector<int>& node_neighbours
)
{
	if (edge_file == NULL ||
		num_nodes < 3)
	{
		fprintf(stderr, "read_edges(): input check failed!\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(edge_file, "rt");
	if (fp == NULL)
	{
		fprintf(stderr, "read_edges(): can't open %s\n", edge_file);
		return -1;
	}
	char str[1024];
	char* ptr;
	fgets(str, 1024, fp);
	long long num_edges = 0;
	num_edges = strtol(str, &ptr, 0);
	if (num_edges <= 0)
	{
		fprintf(stderr, "read_edges(): %s is unknown format!\n", edge_file);
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	edges.clear(); node_neighbours.clear();
	edges.resize(num_edges);
	node_neighbours.resize(num_nodes);
	long end1, end2, edges_number, boundry_marker;
	for (int i = 0; i < num_edges; i++)
	{
		fgets(str, 1024, fp);
		edges_number = strtol(str, &ptr, 0);
		end1 = strtol(ptr, &ptr, 0);
		end2 = strtol(ptr, &ptr, 0);
		boundry_marker = strtol(ptr, &ptr, 0);
		edges[i].end1 = end1;
		edges[i].end2 = end2;
		edges[i].num = i + 1;
		edges[i].isBoundry = (boundry_marker == 1);
		if (end1 < 1 ||
			end1 > num_nodes ||
			end2 < 1 ||
			end2 > num_nodes)
		{
			fprintf(stderr, "read_edges(): endpoints exceed 1~num_nodes!\n");
			if (fp)
			{
				fclose(fp);
				fp = NULL;
			}
			return -1;
		}
		node_neighbours[end1 - 1] += 1;//统计每个节点有多少邻接边
		node_neighbours[end2 - 1] += 1;//统计每个节点有多少邻接边
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int SBAS::init_SBAS_node(
	vector<SBAS_node>& node_array,
	vector<SBAS_edge>& edges,
	vector<int>& node_neighbours
)
{
	if (edges.size() < 3 ||
		node_neighbours.size() < 3)
	{
		fprintf(stderr, "init_SBAS_node(): input check failed!\n");
		return -1;
	}
	int num_nodes = node_neighbours.size();
	int num_edges = edges.size();
	node_array.clear();
	node_array.resize(num_nodes);
	SBAS_node* ptr = NULL;
	for (int i = 0; i < num_nodes; i++)
	{
		ptr = new SBAS_node(node_neighbours[i]);
		node_array[i] = *ptr;
		delete ptr;
		ptr = NULL;
	}

	int* neighbour_ptr = NULL;
	int dummy, ret;
	SBAS_edge tmp;
	for (int i = 0; i < num_edges; i++)
	{
		tmp = edges[i];
		neighbour_ptr = node_array[tmp.end1 - 1].neigh_edges;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;

		neighbour_ptr = node_array[tmp.end2 - 1].neigh_edges;
		while (neighbour_ptr != NULL && *neighbour_ptr != -1)
		{
			neighbour_ptr = neighbour_ptr + 1;
		}
		*neighbour_ptr = i + 1;
	}
	return 0;
}

int SBAS::init_SBAS_triangle(
	const char* ele_file,
	const char* neigh_file,
	vector<SBAS_triangle>& triangle,
	vector<SBAS_edge>& edges,
	vector<SBAS_node>& nodes
)
{
	if (!ele_file ||
		!neigh_file ||
		edges.size() < 3 ||
		nodes.size() < 3)
	{
		fprintf(stderr, "init_SBAS_triangle(): input check failed!\n");
		return -1;
	}
	FILE* fp_ele, * fp_neigh;
	fp_ele = NULL;
	fp_neigh = NULL;
	fp_ele = fopen(ele_file, "rt");
	if (fp_ele == NULL)
	{
		fprintf(stderr, "init_SBAS_triangle(): can't open %s\n", ele_file);
		return -1;
	}
	fp_neigh = fopen(neigh_file, "rt");
	if (fp_neigh == NULL)
	{
		fprintf(stderr, "init_SBAS_triangle(): can't open %s\n", neigh_file);
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		return -1;
	}

	char str[1024];
	char* ptr;
	fgets(str, 1024, fp_ele);
	int num_triangle = strtol(str, &ptr, 0);
	if (num_triangle < 1)
	{
		fprintf(stderr, "init_SBAS_triangle(): number of triangles exceed legal range!\n");
		if (fp_ele)
		{
			fclose(fp_ele);
			fp_ele = NULL;
		}
		if (fp_neigh)
		{
			fclose(fp_neigh);
			fp_neigh = NULL;
		}
		return -1;
	}
	fgets(str, 1024, fp_neigh);
	triangle.clear();
	triangle.resize(num_triangle);

	int p1, p2, p3, neigh1, neigh2, neigh3, num1, num2;
	for (int i = 0; i < num_triangle; i++)
	{
		fgets(str, 1024, fp_ele);
		num1 = strtol(str, &ptr, 0);
		p1 = strtol(ptr, &ptr, 0);
		p2 = strtol(ptr, &ptr, 0);
		p3 = strtol(ptr, &ptr, 0);

		fgets(str, 1024, fp_neigh);
		num2 = strtol(str, &ptr, 0);
		neigh1 = strtol(ptr, &ptr, 0);
		neigh2 = strtol(ptr, &ptr, 0);
		neigh3 = strtol(ptr, &ptr, 0);
		triangle[i].p1 = p1;
		triangle[i].p2 = p2;
		triangle[i].p3 = p3;
		triangle[i].neigh1 = neigh1;
		triangle[i].neigh2 = neigh2;
		triangle[i].neigh3 = neigh3;
		triangle[i].num = num1;
	}
	if (fp_ele)
	{
		fclose(fp_ele);
		fp_ele = NULL;
	}
	if (fp_neigh)
	{
		fclose(fp_neigh);
		fp_neigh = NULL;
	}
	//获取三角形的边序号
	int* ptr_neigh = NULL;
	int num_neigh, count;
	int edge[3];
	memset(edge, 0, sizeof(int) * 3);
	for (int j = 0; j < num_triangle; j++)
	{
		count = 0;
		ptr_neigh = nodes[triangle[j].p1 - 1].neigh_edges;
		num_neigh = nodes[triangle[j].p1 - 1].num_neigh_edges;
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p2) ||
				(edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p3) ||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p2) ||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p3)
				)
			{
				edge[count] = *(ptr_neigh + i);
				count++;
			}
		}
		ptr_neigh = nodes[triangle[j].p2 - 1].neigh_edges;
		num_neigh = nodes[triangle[j].p2 - 1].num_neigh_edges;
		for (int i = 0; i < num_neigh; i++)
		{
			if ((edges[*(ptr_neigh + i) - 1].end1 == triangle[j].p3) ||
				(edges[*(ptr_neigh + i) - 1].end2 == triangle[j].p3))
			{
				edge[count] = *(ptr_neigh + i);
			}
		}
		triangle[j].edge1 = edge[0];
		triangle[j].edge2 = edge[1];
		triangle[j].edge3 = edge[2];
	}

	return 0;
}

int SBAS::compute_spatialTemporal_residue(
	vector<SBAS_node>& nodes,
	vector<SBAS_edge>& edges,
	vector<SBAS_triangle>& triangles
)
{
	if (nodes.size() < 3 ||
		edges.size() < 3 ||
		triangles.size() < 1)
	{
		fprintf(stderr, "compute_spatialTemporal_residue(): input check failed!\n");
		return -1;
	}
	int num_triangle = triangles.size();
	int num_nodes = nodes.size();
	int end1, end2, end3, tmp;
	double x21, y21, x32, y32, direction, delta, residue,x1, x2, x3, y1, y2, y3;
	delta = 0.0;
	for (int i = 0; i < num_triangle; i++)
	{
		end1 = triangles[i].p1;
		end2 = triangles[i].p2;
		end3 = triangles[i].p3;
		x1 = nodes[end1 - 1].x;
		y1 = nodes[end1 - 1].y;
		x2 = nodes[end2 - 1].x;
		y2 = nodes[end2 - 1].y;
		x3 = nodes[end3 - 1].x;
		y3 = nodes[end3 - 1].y;

		x21 = x2 - x1;
		y21 = y2 - y1;
		x32 = x3 - x2;
		y32 = y3 - y2;

		direction = x1 * y2 - x2 * y1;

		//edge1处于end1和end2之间
		if ((edges[triangles[i].edge1 - 1].end1 == end1 && edges[triangles[i].edge1 - 1].end2 == end2) ||
			(edges[triangles[i].edge1 - 1].end1 == end2 && edges[triangles[i].edge1 - 1].end2 == end1)
			)
		{
			if (nodes[end2 - 1].B_temporal > nodes[end1 - 1].B_temporal)
			{
				delta += edges[triangles[i].edge1 - 1].phase_gradient;
			}
			else
			{
				delta -= edges[triangles[i].edge1 - 1].phase_gradient;
			}

		    //edge2处于end1和end3之间
			if (edges[triangles[i].edge2 - 1].end1 == end1 || edges[triangles[i].edge2 - 1].end2 == end1)
			{
				if (nodes[end1 - 1].B_temporal > nodes[end3 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}


				//edge3处于end2和end3之间
				if (nodes[end3 - 1].B_temporal > nodes[end2 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}
			else//edge2处于end2和end3之间
			{
				if (nodes[end3 - 1].B_temporal > nodes[end2 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}
				//edge3处于end1和end3之间
				if (nodes[end1 - 1].B_temporal > nodes[end3 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}
		}
		//edge1处于end1和end3之间
		else if ((edges[triangles[i].edge1 - 1].end1 == end1 && edges[triangles[i].edge1 - 1].end2 == end3) ||
			(edges[triangles[i].edge1 - 1].end1 == end3 && edges[triangles[i].edge1 - 1].end2 == end1)
			)
		{
			if (nodes[end1 - 1].B_temporal > nodes[end3 - 1].B_temporal)
			{
				delta += edges[triangles[i].edge1 - 1].phase_gradient;
			}
			else
			{
				delta -= edges[triangles[i].edge1 - 1].phase_gradient;
			}
			//edge2处于end1和end2之间
			if (edges[triangles[i].edge2 - 1].end1 == end1 || edges[triangles[i].edge2 - 1].end2 == end1)
			{
				if (nodes[end2 - 1].B_temporal > nodes[end1 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}
				//edge3处于end2和end3之间
				if (nodes[end3 - 1].B_temporal > nodes[end2 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}

			//edge2处于end2和end3之间
			else
			{
				if (nodes[end3 - 1].B_temporal > nodes[end2 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}

				//edge3处于end1和end2之间
				if (nodes[end2 - 1].B_temporal > nodes[end1 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}
		}
		//edge1处于end2和end3之间
		else
		{
			if (nodes[end3 - 1].B_temporal > nodes[end2 - 1].B_temporal)
			{
				delta += edges[triangles[i].edge1 - 1].phase_gradient;
			}
			else
			{
				delta -= edges[triangles[i].edge1 - 1].phase_gradient;
			}
			//edge2处于end1和end3之间
			if (edges[triangles[i].edge2 - 1].end1 == end3 || edges[triangles[i].edge2 - 1].end2 == end3)
			{
				if (nodes[end1 - 1].B_temporal > nodes[end3 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}

				//edge3处于end1和end2之间
				if (nodes[end2 - 1].B_temporal > nodes[end1 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}
			//edge2处于end1和end2之间
			else
			{
				if (nodes[end2 - 1].B_temporal > nodes[end1 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge2 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge2 - 1].phase_gradient;
				}

				//edge3处于end1和end3之间
				if (nodes[end1 - 1].B_temporal > nodes[end3 - 1].B_temporal)
				{
					delta += edges[triangles[i].edge3 - 1].phase_gradient;
				}
				else
				{
					delta -= edges[triangles[i].edge3 - 1].phase_gradient;
				}
			}
		}
		

		double res = round(delta / 2.0 / PI);
		if (direction < 0.0)//在目标三角形中逆残差方向(残差积分方向定义为逆时针方向)
		{
			triangles[i].residue = -res;
		}
		else
		{
			triangles[i].residue = res;
		}
	}
	return 0;
}

int SBAS::generate_interferograms(
	vector<string>& SLCH5Files,
	vector<SBAS_edge>& edges, 
	vector<SBAS_node>& nodes,
	int multilook_az, 
	int multilook_rg, 
	const char* ifgSavePath
)
{
	if (edges.size() < 3 ||
		nodes.size() < 3 ||
		multilook_rg < 1 ||
		multilook_az < 1 ||
		SLCH5Files.size() != nodes.size() ||
		!ifgSavePath
		)
	{
		fprintf(stderr, "generate_interferograms(): input check failed!\n");
		return -1;
	}
	int ret;
	Utils util; FormatConversion conversion;
	string path(ifgSavePath), h5file;
	std::replace(path.begin(), path.end(), '/', '\\');
	ComplexMat master, slave;
	Mat phase;
	char str[256];
	int master_ix, slave_ix, offset_row, offset_col;
	double B_temporal, B_spatial;
	for (int i = 0; i < edges.size(); i++)
	{
		if (nodes[edges[i].end1 - 1].B_temporal > nodes[edges[i].end2 - 1].B_temporal)
		{
			master_ix = edges[i].end1;slave_ix = edges[i].end2;
			B_temporal = nodes[edges[i].end1 - 1].B_temporal - nodes[edges[i].end2 - 1].B_temporal;
			B_spatial = nodes[edges[i].end1 - 1].B_spatial > nodes[edges[i].end2 - 1].B_spatial;
		}
		else
		{
			master_ix = edges[i].end2; slave_ix = edges[i].end1;
			B_temporal = nodes[edges[i].end2 - 1].B_temporal - nodes[edges[i].end1 - 1].B_temporal;
			B_spatial = nodes[edges[i].end2 - 1].B_spatial > nodes[edges[i].end1 - 1].B_spatial;
		}
		ret = conversion.read_slc_from_h5(SLCH5Files[master_ix - 1].c_str(), master);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		ret = conversion.read_int_from_h5(SLCH5Files[master_ix - 1].c_str(), "offset_row", &offset_row);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		ret = conversion.read_int_from_h5(SLCH5Files[master_ix - 1].c_str(), "offset_col", &offset_col);
		if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
		ret = conversion.read_slc_from_h5(SLCH5Files[slave_ix - 1].c_str(), slave);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (master.type() != CV_64F) master.convertTo(master, CV_64F);
		if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
		ret = util.Multilook(master, slave, multilook_rg, multilook_az, phase);
		if (return_check(ret, "Multilook()", error_head)) return -1;
		sprintf(str, "\\%d.h5", i + 1);
		h5file = path + str;
		ret = conversion.creat_new_h5(h5file.c_str());
		ret = conversion.write_array_to_h5(h5file.c_str(), "phase", phase);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "offset_row", offset_row);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "offset_col", offset_col);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "range_len", phase.cols);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "azimuth_len", phase.rows);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "multilook_rg", multilook_rg);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_int_to_h5(h5file.c_str(), "multilook_az", multilook_az);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_double_to_h5(h5file.c_str(), "B_temporal", B_temporal);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.write_double_to_h5(h5file.c_str(), "B_spatial", B_spatial);
		if (return_check(ret, "write_int_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(SLCH5Files[master_ix - 1].c_str(), h5file.c_str());
		if (return_check(ret, "Copy_para_from_h5_2_h5()", error_head)) return -1;

	}
	return 0;
}

int SBAS::generate_high_coherence_mask(
	vector<string>& phaseFiles, 
	int wndsize_rg,
	int wndsize_az, 
	double coherence_thresh,
	Mat& mask
)
{
	if (phaseFiles.size() < 1 ||
		wndsize_az % 2 != 1 ||
		wndsize_rg % 2 != 1
		)
	{
		fprintf(stderr, "generate_high_coherence_mask(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Utils util;
	int ret;
	Mat phase, coherence;
	coherence_thresh = coherence_thresh < 0.3 ? 0.3 : coherence_thresh;
	coherence_thresh = coherence_thresh > 0.95 ? 0.95 : coherence_thresh;
	int num_images = phaseFiles.size();
	for (int i = 0; i < num_images; i++)
	{
		ret = conversion.read_array_from_h5(phaseFiles[i].c_str(), "phase", phase);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		int rows = phase.rows; int cols = phase.cols;
		ret = util.phase_coherence(phase, wndsize_rg, wndsize_az, coherence);
		if (return_check(ret, "phase_coherence()", error_head)) return -1;
		if (i == 0)
		{
			Mat tmp = Mat::ones(rows, cols, CV_32S);
			tmp.copyTo(mask);
		}
#pragma omp parallel for schedule(guided)
		for (int j = 0; j < rows; j++)
		{
			for (int k = 0; k < cols; k++)
			{
				if (coherence.at<double>(j, k) < coherence_thresh) mask.at<int>(j, k) = 0.0;
			}
		}
	}
	return 0;
}

int SBAS::write_high_coherence_node(Mat& mask, const char* filename)
{
	if (filename == NULL ||
		mask.rows < 2 ||
		mask.cols < 2 ||
		mask.channels() != 1 ||
		mask.type() != CV_32S
		)
	{
		fprintf(stderr, "write_high_coherence_node(): input check failed!\n\n");
		return -1;
	}
	FILE* fp = NULL;
	fp = fopen(filename, "wt");
	if (fp == NULL)
	{
		fprintf(stderr, "write_high_coherence_node(): can't open %s\n", filename);
		return -1;
	}
	int nonzero = cv::countNonZero(mask);
	if (nonzero <= 0)
	{
		fprintf(stderr, "write_high_coherence_node(): no node exist!\n");
		if (fp)
		{
			fclose(fp);
			fp = NULL;
		}
		return -1;
	}
	int dim, attr1, attr2, count;
	dim = 2;
	attr1 = 0;
	attr2 = 0;
	fprintf(fp, "%d %d %d %d\n", nonzero, dim, attr1, attr2);
	int rows = mask.rows;
	int cols = mask.cols;
	count = 1;
	double x, y;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				x = double(i + 1);
				y = double(j + 1);
				fprintf(fp, "%d %lf %lf\n", count, x, y);
				count++;
			}
		}
	}
	if (fp)
	{
		fclose(fp);
		fp = NULL;
	}
	return 0;
}

int SBAS::set_high_coherence_node_coordinate(
	Mat& mask,
	vector<SBAS_node>& nodes
)
{
	if (mask.empty() ||
		mask.type() != CV_32S ||
		nodes.size() < 3
		)
	{
		fprintf(stderr, "set_high_coherence_node_coordinate(): input check failed!\n");
		return -1;
	}
	int rows = mask.rows;
	int cols = mask.cols;
	int nonzero = cv::countNonZero(mask);
	int num_nodes = nodes.size();
	if (num_nodes != nonzero)
	{
		fprintf(stderr, "set_high_coherence_node_coordinate(): nodes and mask mismatch!\n");
		return -1;
	}
	int count = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				nodes[count].x = j + 1;
				nodes[count].y = i + 1;
				count++;
			}
		}
	}
	return 0;
}

int SBAS::set_high_coherence_node_phase(
	Mat& mask, 
	vector<SBAS_node>& nodes, 
	Mat& phase
)
{
	if (mask.empty() ||
		mask.type() != CV_32S ||
		nodes.size() < 3 || 
		phase.rows != mask.rows ||
		phase.cols != mask.cols ||
		phase.type() != CV_64F
		)
	{
		fprintf(stderr, "set_high_coherence_node_phase(): input check failed!\n");
		return -1;
	}
	int rows = mask.rows;
	int cols = mask.cols;
	int nonzero = cv::countNonZero(mask);
	int num_nodes = nodes.size();
	if (num_nodes != nonzero)
	{
		fprintf(stderr, "set_high_coherence_node_phase(): nodes and mask mismatch!\n");
		return -1;
	}
	int count = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mask.at<int>(i, j) > 0)
			{
				nodes[count].phase = phase.at<double>(i, j);
				count++;
			}
		}
	}
	return 0;
}
