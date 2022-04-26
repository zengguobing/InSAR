//#include<io.h>
//#include<Windows.h>
#include"pch.h"
#include"gdal_priv.h"
#include"..\include\FormatConversion.h"
#include"..\include\Registration.h"
#include"Utils.h"
//#include<atlconv.h>
//#include<tchar.h>
#include<urlmon.h>
#pragma comment(lib,"URlmon")

#ifdef _DEBUG
#pragma comment(lib,"ComplexMat_d.lib")
#pragma comment(lib,"Registration_d.lib")
#pragma comment(lib,"Utils_d.lib")
#else
#pragma comment(lib,"ComplexMat.lib")
#pragma comment(lib,"Registration.lib")
#pragma comment(lib,"Utils.lib")

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

inline float ReverseFloat(const float inFloat)
{
	float retVal;
	unsigned char* floatToConvert = (unsigned char*)&inFloat;
	unsigned char* returnFloat = (unsigned char*)&retVal;

	// swap the bytes into a temporary buffer
	returnFloat[0] = floatToConvert[3];
	returnFloat[1] = floatToConvert[2];
	returnFloat[2] = floatToConvert[1];
	returnFloat[3] = floatToConvert[0];

	return retVal;
}

int UTC2GPS(const char* utc_time, double* gps_time)
{
	if (utc_time == NULL || gps_time == NULL)
	{
		fprintf(stderr, "UTC2GPS(): input check failed!\n");
		return -1;
	}
	int ret, year, month, day, hour, minute, second, s;
	double sec;
	ret = sscanf(utc_time, "%d-%d-%dT%d:%d:%lf\n", &year, &month, &day, &hour, &minute, &sec);
	if (ret != 6)
	{
		fprintf(stderr, "UTC2GPS(): %s: unknown format!\n", utc_time);
		return -1;
	}
	second = int(floor(sec));
	sec = sec - (double)second;
	tm TM;
	TM.tm_year = year - 1900;
	TM.tm_mon = month - 1;
	TM.tm_mday = day;
	TM.tm_hour = hour;
	TM.tm_min = minute;
	TM.tm_sec = second;
	*gps_time = double(mktime(&TM) - 315964809) + sec;
	return 0;
}

FormatConversion::FormatConversion()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using parallel computing in function: ");
}

FormatConversion::~FormatConversion()
{

}

int FormatConversion::utc2gps(const char* utc_time, double* gps_time)
{
	if (utc_time == NULL || gps_time == NULL)
	{
		fprintf(stderr, "utc2gps(): input check failed!\n");
		return -1;
	}
	int ret, year, month, day, hour, minute, second, s;
	double sec;
	ret = sscanf(utc_time, "%d-%d-%dT%d:%d:%lf\n", &year, &month, &day, &hour, &minute, &sec);
	if (ret != 6)
	{
		fprintf(stderr, "utc2gps(): %s: unknown format!\n", utc_time);
		return -1;
	}
	second = int(floor(sec));
	sec = sec - (double)second;
	tm TM;
	TM.tm_year = year - 1900;
	TM.tm_mon = month - 1;
	TM.tm_mday = day;
	TM.tm_hour = hour;
	TM.tm_min = minute;
	TM.tm_sec = second;
	*gps_time = double(mktime(&TM) - 315964809) + sec;
	return 0;
}

int FormatConversion::creat_new_h5(const char* filename)
{
	if (filename == NULL)
	{
		fprintf(stderr, "creat_new_h5(): invalid filename!\n");
		return -1;
	}
	hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "creat_new_h5(): failed to create %s!\n", filename);
		return -1;
	}
	herr_t err = H5Fclose(file_id);
	if (err < 0)
	{
		fprintf(stderr, "creat_new_h5(): unable to close %s!\n", filename);
		return -1;
	}
	return 0;
}

int FormatConversion::write_array_to_h5(const char* filename, const char* dataset_name, const Mat& input_array)
{
	if (filename == NULL ||
		dataset_name == NULL ||
		input_array.empty() ||
		input_array.channels() != 1 ||
		(input_array.type() != CV_64F && input_array.type() != CV_16S && input_array.type() != CV_32S && input_array.type() != CV_32F)
		)
	{
		fprintf(stderr, "write_array_to_h5(): input check  failed!\n");
		return -1;
	}
	hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "write_array_to_h5(): can't open %s\n", filename);
		return -1;
	}
	string s = "/"; 
	s.append(dataset_name);	
	hid_t dataset_id;
	if ((H5Lexists(file_id, dataset_name, H5P_DEFAULT)) == 0)
	{
		hsize_t dims[2];
		dims[0] = input_array.rows;
		dims[1] = input_array.cols;
		hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
		if (input_array.type() == CV_16S)
		{
			dataset_id = H5Dcreate(file_id, s.c_str(), H5T_NATIVE_INT16, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		else if(input_array.type() == CV_64F)
		{
			dataset_id = H5Dcreate(file_id, s.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		else if(input_array.type() == CV_32S)
		{
			dataset_id = H5Dcreate(file_id, s.c_str(), H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		else
		{
			dataset_id = H5Dcreate(file_id, s.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		
		if (dataset_id < 0)
		{
			fprintf(stderr, "write_array_to_h5(): failed to create dataset %s !\n", dataset_name);
			H5Sclose(dataspace_id);
			H5Fclose(file_id);
			return -1;
		}
		herr_t status;
		if (input_array.type() == CV_16S)
		{
			status = H5Dwrite(dataset_id, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)input_array.data);
		}
		else if(input_array.type() == CV_64F)
		{
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)input_array.data);
		}
		else if(input_array.type() == CV_32S)
		{
			status = H5Dwrite(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)input_array.data);
		}
		else
		{
			status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)input_array.data);
		}
		if (status < 0)
		{
			fprintf(stderr, "write_array_to_h5(): failed to write to dataset %s !\n", dataset_name);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			H5Fclose(file_id);
			return -1;
		}
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Fclose(file_id);
		
	}
	else
	{
		fprintf(stderr, "write_array_to_h5(): dataset %s already exists!\n", dataset_name);
		H5Fclose(file_id);
		return -1;
	}
	return 0;
}

int FormatConversion::write_double_to_h5(const char* h5File, const char* datasetName, double data)
{
	if (!h5File || !datasetName)
	{
		fprintf(stderr, "write_double_to_h5(): input check failed!\n");
		return -1;
	}
	Mat tmp(1, 1, CV_64F);
	tmp.at<double>(0, 0) = data;
	int ret = write_array_to_h5(h5File, datasetName, tmp);
	if (return_check(ret, "write_array_to_h5", error_head)) return -1;
	return 0;
}

int FormatConversion::write_int_to_h5(const char* h5File, const char* datasetName, int data)
{
	if (!h5File || !datasetName)
	{
		fprintf(stderr, "write_double_to_h5(): input check failed!\n");
		return -1;
	}
	Mat tmp(1, 1, CV_32S);
	tmp.at<int>(0, 0) = data;
	int ret = write_array_to_h5(h5File, datasetName, tmp);
	if (return_check(ret, "write_array_to_h5", error_head)) return -1;
	return 0;
}

int FormatConversion::read_array_from_h5(const char* filename, const char* dataset_name, Mat& out_array)
{
	if (filename == NULL ||
		dataset_name == NULL
		)
	{
		fprintf(stderr, "read_array_from_h5(): input check  failed!\n");
		return -1;
	}
	hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "read_array_from_h5(): failed to open %s!\n", filename);
		return -1;
	}
	string s = "/";
	s.append(dataset_name);
	hid_t dataset_id = H5Dopen(file_id, s.c_str(), H5P_DEFAULT);
	if (dataset_id < 0)
	{
		fprintf(stderr, "read_array_from_h5(): failed to open dataset %s!\n", dataset_name);
		H5Fclose(file_id);
		return -1;
	}
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id < 0)
	{
		fprintf(stderr, "read_array_from_h5(): failed to open dataspace of %s!\n", dataset_name);
		H5Dclose(dataset_id);
		H5Fclose(file_id);
		return -1;
	}
	hsize_t dims[2];
	int ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
	hid_t type = H5Dget_type(dataset_id);
	herr_t status;
	if (H5Tequal(type, H5T_NATIVE_INT16) > 0)
	{
		out_array.create(dims[0], dims[1], CV_16S);
		status = H5Dread(dataset_id, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)out_array.data);
	}
	else if(H5Tequal(type, H5T_NATIVE_DOUBLE) > 0)
	{
		out_array.create(dims[0], dims[1], CV_64F);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)out_array.data);
	}
	else if (H5Tequal(type, H5T_NATIVE_FLOAT) > 0)
	{
		out_array.create(dims[0], dims[1], CV_32F);
		status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)out_array.data);
	}
	else if (H5Tequal(type, H5T_NATIVE_INT) > 0)
	{
		out_array.create(dims[0], dims[1], CV_32S);
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)out_array.data);
	}
	if (status < 0)
	{
		fprintf(stderr, "read_array_from_h5(): failed to read from %s!\n", dataset_name);
		H5Dclose(dataset_id);
		H5Sclose(space_id);
		H5Fclose(file_id);
		H5Tclose(type);
		return -1;
	}
	H5Dclose(dataset_id);
	H5Sclose(space_id);
	H5Fclose(file_id);
	H5Tclose(type);
	return 0;
}

int FormatConversion::read_double_from_h5(const char* h5File, const char* datasetName, double* data)
{
	if (!h5File || !datasetName)
	{
		fprintf(stderr, "read_double_from_h5(): input check failed!\n");
		return -1;
	}
	Mat tmp;
	int ret = read_array_from_h5(h5File, datasetName, tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (tmp.type() != CV_64F) tmp.convertTo(tmp, CV_64F);
	*data = tmp.at<double>(0, 0);
	return 0;
}

int FormatConversion::read_int_from_h5(const char* h5File, const char* datasetName, int* data)
{
	if (!h5File || !datasetName)
	{
		fprintf(stderr, "read_int_from_h5(): input check failed!\n");
		return -1;
	}
	Mat tmp;
	int ret = read_array_from_h5(h5File, datasetName, tmp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	if (tmp.type() != CV_32S) tmp.convertTo(tmp, CV_32S);
	*data = tmp.at<int>(0, 0);
	return 0;
}

int FormatConversion::read_subarray_from_h5(const char* filename, const char* dataset_name, int offset_row, int offset_col, int rows_subarray, int cols_subarray, Mat& out_array)
{
	if (filename == NULL ||
		dataset_name == NULL ||
		offset_row < 0 ||
		offset_col < 0 ||
		rows_subarray < 1 ||
		cols_subarray < 1
		)
	{
		fprintf(stderr, "read_subarray_from_h5(): input check failed!\n");
		return -1;
	}

	hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "read_subarray_from_h5(): failed to open %s!\n", filename);
		return -1;
	}
	string s = "/";
	s.append(dataset_name);
	if (0 == H5Lexists(file_id, dataset_name, H5P_DEFAULT))
	{
		fprintf(stderr, "read_subarray_from_h5(): dataset %s doesn't exist!\n", dataset_name);
		H5Fclose(file_id);
		return -1;
	}
	hid_t dataset_id = H5Dopen(file_id, s.c_str(), H5P_DEFAULT);
	// 读取文件中dataset的dataspace空间
	hid_t dataspace_id = H5Dget_space(dataset_id);
	hsize_t dim[2];
	int ndims = H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
	if ((int)dim[0] < (offset_row - 1) ||
		(int)dim[1] < (offset_col - 1) ||
		(offset_row + rows_subarray) > (int)dim[0] ||
		(offset_col + cols_subarray) > (int)dim[1])
	{
		fprintf(stderr, "read_subarray_from_h5(): invalide subarray index!\n");
		H5Fclose(file_id);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		return -1;
	}
	hid_t type = H5Dget_type(dataset_id);
	if (0 < H5Tequal(type, H5T_NATIVE_INT16))
	{
		out_array.create(rows_subarray, cols_subarray, CV_16S);
	}
	else if(0 < H5Tequal(type, H5T_NATIVE_DOUBLE))
	{
		out_array.create(rows_subarray, cols_subarray, CV_64F);
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_FLOAT))
	{
		out_array.create(rows_subarray, cols_subarray, CV_32F);
	}
	else
	{
		fprintf(stderr, "read_subarray_from_h5(): datatype not support yet!\n");
		H5Fclose(file_id);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Tclose(type);
		return -1;
	}
	// 定义子集四大件，补偿，个数，间隔和块大小
	hsize_t count[2];   // 块的大小
	hsize_t offset[2];  // 补偿，即开始位置
	hsize_t stride[2];  // 间隔
	hsize_t block[2];   // 块的个数

	offset[0] = offset_row;
	offset[1] = offset_col;

	count[0] = rows_subarray;
	count[1] = cols_subarray;

	stride[0] = 1;
	stride[1] = 1;

	block[0] = 1;
	block[1] = 1;

	// 创建内存中的dataspce空间
	hsize_t dimsm[2];
	dimsm[0] = rows_subarray;
	dimsm[1] = cols_subarray;
	hid_t memspace_id;
	memspace_id = H5Screate_simple(2, dimsm, NULL);
	H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
	if (0 < H5Tequal(type, H5T_NATIVE_INT16))
	{
		H5Dread(dataset_id, H5T_NATIVE_INT16, memspace_id, dataspace_id, H5P_DEFAULT, out_array.data);
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_DOUBLE))
	{
		H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, out_array.data);
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_FLOAT))
	{
		H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, out_array.data);
	}
	else
	{

	}
	H5Fclose(file_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(type);
	return 0;
}

int FormatConversion::write_subarray_to_h5(const char* h5_filename, const char* dataset_name, Mat& subarray, int offset_row, int offset_col, int rows_subarray, int cols_subarray)
{
	if (h5_filename == NULL ||
		dataset_name == NULL ||
		offset_row < 0 ||
		offset_col < 0 ||
		rows_subarray < 1 ||
		cols_subarray < 1
		)
	{
		fprintf(stderr, "write_subarray_to_h5(): input check failed!\n");
		return -1;
	}

	hid_t file_id = H5Fopen(h5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "write_subarray_to_h5(): failed to open %s!\n", h5_filename);
		return -1;
	}
	string s = "/";
	s.append(dataset_name);
	if (0 == H5Lexists(file_id, dataset_name, H5P_DEFAULT))
	{
		fprintf(stderr, "write_subarray_to_h5(): dataset %s doesn't exist!\n", dataset_name);
		H5Fclose(file_id);
		return -1;
	}
	hid_t dataset_id = H5Dopen(file_id, s.c_str(), H5P_DEFAULT);
	// 读取文件中dataset的dataspace空间
	hid_t dataspace_id = H5Dget_space(dataset_id);
	hsize_t dim[2];
	int ndims = H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
	if ((int)dim[0] < (offset_row - 1) ||
		(int)dim[1] < (offset_col - 1) ||
		(offset_row + rows_subarray) > (int)dim[0] ||
		(offset_col + cols_subarray) > (int)dim[1])
	{
		fprintf(stderr, "write_subarray_to_h5(): invalide subarray index!\n");
		H5Fclose(file_id);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		return -1;
	}
	hid_t type = H5Dget_type(dataset_id);
	if (0 < H5Tequal(type, H5T_NATIVE_INT16))
	{
		if (subarray.type() != CV_16S)
		{
			fprintf(stderr, "write_subarray_to_h5(): datatype mismatch!\n");
			H5Fclose(file_id);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			H5Tclose(type);
			return -1;
		}
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_DOUBLE))
	{
		if (subarray.type() != CV_64F)
		{
			fprintf(stderr, "write_subarray_to_h5(): datatype mismatch!\n");
			H5Fclose(file_id);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			H5Tclose(type);
			return -1;
		}
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_FLOAT))
	{
		if (subarray.type() != CV_32F)
		{
			fprintf(stderr, "write_subarray_to_h5(): datatype mismatch!\n");
			H5Fclose(file_id);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
			H5Tclose(type);
			return -1;
		}
	}
	else
	{
		fprintf(stderr, "write_subarray_to_h5(): datatype not support yet!\n");
		H5Fclose(file_id);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Tclose(type);
		return -1;
	}
	// 定义子集四大件，补偿，个数，间隔和块大小
	hsize_t count[2];   // 块的大小
	hsize_t offset[2];  // 补偿，即开始位置
	hsize_t stride[2];  // 间隔
	hsize_t block[2];   // 块的个数

	offset[0] = offset_row;
	offset[1] = offset_col;

	count[0] = rows_subarray;
	count[1] = cols_subarray;

	stride[0] = 1;
	stride[1] = 1;

	block[0] = 1;
	block[1] = 1;

	// 创建内存中的dataspce空间
	hsize_t dimsm[2];
	dimsm[0] = rows_subarray;
	dimsm[1] = cols_subarray;
	hid_t memspace_id;
	memspace_id = H5Screate_simple(2, dimsm, NULL);
	H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
	if (0 < H5Tequal(type, H5T_NATIVE_INT16))
	{
		H5Dwrite(dataset_id, H5T_NATIVE_INT16, memspace_id, dataspace_id, H5P_DEFAULT, subarray.data);
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_DOUBLE))
	{
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, subarray.data);
	}
	else if (0 < H5Tequal(type, H5T_NATIVE_FLOAT))
	{
		H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, subarray.data);
	}
	else
	{

	}
	H5Fclose(file_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(type);
	return 0;
}

int FormatConversion::write_str_to_h5(const char* filename, const char* dataset_name, const char* Str)
{
	if (filename == NULL ||
		dataset_name == NULL||
		Str == NULL
		)
	{
		fprintf(stderr, "write_str_to_h5(): input check  failed!\n");
		return -1;
	}
	hid_t file_id, dataset_id, space_id, filetype, memtype;
	herr_t status;
	hsize_t dims[1] = { 1 };
	size_t sdim;
	string s("/");
	string str(Str);
	s.append(dataset_name);
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "write_str_to_h5(): failed to open %s !\n", filename);
		return -1;
	}
	filetype = H5Tcopy(H5T_FORTRAN_S1);
	status = H5Tset_size(filetype, str.length());
	memtype = H5Tcopy(H5T_C_S1);
	status = H5Tset_size(memtype, str.length());
	space_id = H5Screate_simple(1, dims, NULL);
	if ((H5Lexists(file_id, dataset_name, H5P_DEFAULT)) == 0)
	{
		dataset_id = H5Dcreate(file_id, s.c_str(), filetype, space_id, H5P_DEFAULT, H5P_DEFAULT,
			H5P_DEFAULT);
		if (dataset_id < 0)
		{
			fprintf(stderr, "write_str_to_h5(): failed to create dataset %s !\n", dataset_name);
			status = H5Sclose(space_id);
			status = H5Tclose(filetype);
			status = H5Tclose(memtype);
			status = H5Fclose(file_id);
			return -1;
		}
		status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str.c_str());
		if (status < 0)
		{
			fprintf(stderr, "write_str_to_h5(): failed to write to dataset %s !\n", dataset_name);
			status = H5Dclose(dataset_id);
			status = H5Sclose(space_id);
			status = H5Tclose(filetype);
			status = H5Tclose(memtype);
			status = H5Fclose(file_id);
			return -1;
		}
	}
	else
	{
		fprintf(stderr, "write_str_to_h5(): dataset %s already exists!\n", dataset_name);
		status = H5Sclose(space_id);
		status = H5Tclose(filetype);
		status = H5Tclose(memtype);
		status = H5Fclose(file_id);
		return -1;
	}
	status = H5Dclose(dataset_id);
	status = H5Sclose(space_id);
	status = H5Tclose(filetype);
	status = H5Tclose(memtype);
	status = H5Fclose(file_id);
	return 0;
}

int FormatConversion::read_str_from_h5(const char* filename, const char* dataset_name, string& Str)
{
	if (filename == NULL ||
		dataset_name == NULL
		)
	{
		fprintf(stderr, "read_str_from_h5(): input check failed!\n");
		return -1;
	}
	hid_t file_id, dataset_id, space_id, filetype, memtype;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		fprintf(stderr, "read_str_from_h5(): failed to open %s !\n", filename);
		return -1;
	}
	string s("/");
	s.append(dataset_name);
	dataset_id = H5Dopen(file_id, s.c_str(), H5P_DEFAULT);
	if (dataset_id < 0)
	{
		fprintf(stderr, "read_str_from_h5(): failed to open dataset %s !\n", dataset_name);
		H5Fclose(file_id);
		return -1;
	}
	filetype = H5Dget_type(dataset_id);
	size_t sdim = H5Tget_size(filetype) + 1;
	space_id = H5Dget_space(dataset_id);
	hsize_t dims[1] = { 1 };
	int ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
	memtype = H5Tcopy(H5T_C_S1);
	status = H5Tset_size(memtype, sdim);
	char* rdata = (char*)malloc(dims[0] * sdim * sizeof(char));
	status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
	if (status < 0)
	{
		fprintf(stderr, "read_str_from_h5(): failed to read from dataset %s!\n", dataset_name);
		H5Sclose(space_id);
		H5Dclose(dataset_id);
		H5Tclose(filetype);
		H5Tclose(memtype);
		H5Fclose(file_id);
		free(rdata);
	}
	string tmp(rdata);
	free(rdata);
	Str = tmp;
	H5Sclose(space_id);
	H5Dclose(dataset_id);
	H5Tclose(filetype);
	H5Tclose(memtype);
	H5Fclose(file_id);
	return 0;
}

int FormatConversion::write_slc_to_h5(const char* filename, const ComplexMat& slc)
{
	if (filename == NULL ||
		slc.isempty()/*||
		slc.type() != CV_64F*/
		)
	{
		fprintf(stderr, "write_slc_to_h5(): input check failed!\n");
		return -1;
	}
	int ret;
	ret = write_array_to_h5(filename, "s_re", slc.re);
	if (return_check(ret, "write_slc_to_h5()", this->error_head)) return -1;
	ret = write_array_to_h5(filename, "s_im", slc.im);
	if (return_check(ret, "write_slc_to_h5()", this->error_head)) return -1;
	return 0;
}

int FormatConversion::read_slc_from_h5(const char* filename, ComplexMat& slc)
{
	if (filename == NULL)
	{
		fprintf(stderr, "read_slc_from_h5(): input check failed!\n");
		return -1;
	}
	int ret;
	ret = read_array_from_h5(filename, "s_re", slc.re);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	ret = read_array_from_h5(filename, "s_im", slc.im);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	return 0;
}

int FormatConversion::read_slc_from_TSXcos(const char* filename, ComplexMat& slc)
{
	if (filename == NULL)
	{
		fprintf(stderr, "read_slc_from_TSXcos(): input check failed!\n");
		return -1;
	}
	GDALAllRegister();	//注册已知驱动
	GDALDataset* poDataset = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);	//打开cos文件
	if (poDataset == NULL)
	{
		fprintf(stderr, "read_slc_from_TSXcos(): failed to open %s!\n", filename);
		GDALDestroyDriverManager();
		return -1;
	}
	int nBand = poDataset->GetRasterCount();	//获取波段数（cos应为1）
	int xsize = 0;
	int ysize = 0;
	if (nBand == 1)
	{
		GDALRasterBand* poBand = poDataset->GetRasterBand(1);	//获取指向波段1的指针
		xsize = poBand->GetXSize();		//cols
		ysize = poBand->GetYSize();		//rows
		if (xsize < 0 || ysize < 0)
		{
			fprintf(stderr, "read_slc_from_TSXcos(): band rows and cols error!\n");
			GDALClose(poDataset);
			GDALDestroyDriverManager();
			return -1;
		}
		GDALDataType dataType = poBand->GetRasterDataType();	//数据存储类型，cos应为GDT_CInt16
		int* pbuf = NULL;
		pbuf = (int*)malloc(sizeof(int) * xsize * ysize);		//分配数据指针空间
		if (!pbuf)
		{
			fprintf(stderr, "read_slc_from_TSXcos(): out of memory!\n");
			GDALClose(poDataset);
			GDALDestroyDriverManager();
			return -1;
		}
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize, pbuf, xsize, ysize, dataType, 0, 0);		//读取复图像数据到pbuf中
		int i, j;
		//ComplexMat CMat;
		//CMat.re = Mat::zeros(ysize, xsize, CV_64F);
		//CMat.im = Mat::zeros(ysize, xsize, CV_64F);
		slc.re.create(ysize, xsize, CV_16S);
		slc.im.create(ysize, xsize, CV_16S);
		for (i = 0; i < ysize; i++)
			for (j = 0; j < xsize; j++)
			{
				/*将数据按照实部、虚部读入Mat中
				由于cos按照大端存储，RasterIO会自动转化为小端，但会导致实部虚部位置颠倒
				因此高16位为虚部，低16位为实部
				*/
				slc.re.ptr<short>(i)[j] = (pbuf[j + i * xsize] << 16) >> 16;
				slc.im.ptr<short>(i)[j] = (pbuf[j + i * xsize] >> 16);
			}
		if (pbuf)
		{
			free(pbuf);
			pbuf = NULL;
		}
		GDALClose(poDataset);
		GDALDestroyDriverManager();
	}
	else
	{
		fprintf(stderr, "read_slc_from_TSXcos(): number of Bands != 1\n");
		GDALClose(poDataset);
		GDALDestroyDriverManager();
		return -1;
	}
	
	return 0;
}

int FormatConversion::TSX2h5(const char* cosar_filename, const char* xml_filename, const char* GEOREF_filename, const char* dst_h5_filename)
{
	if (cosar_filename == NULL ||
		xml_filename == NULL ||
		GEOREF_filename == NULL ||
		dst_h5_filename == NULL
		)
	{
		fprintf(stderr, "TSX2h5(): input check failed!\n");
		return -1;
	}
	/*
	* 检查h5文件是否已经存在
	*/

	int ret;
	ret = creat_new_h5(dst_h5_filename);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;

	/*
	* 写入slc数据
	*/

	ComplexMat slc;
	int rows, cols;
	ret = read_slc_from_TSXcos(cosar_filename, slc);
	if (return_check(ret, "read_slc_from_TSXcos()", error_head)) return -1;
	rows = slc.GetRows(); cols = slc.GetCols();
	ret = write_array_to_h5(dst_h5_filename, "s_re", slc.re);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "s_im", slc.im);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	slc.re.release();
	slc.im.release();

	/*
	* 写入控制点数据
	*/

	Mat gcps;
	XMLFile xmldoc;
	ret = xmldoc.XMLFile_load(GEOREF_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	ret = xmldoc.get_gcps_from_TSX(gcps);
	if (return_check(ret, "get_gcps_from_TSX", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "gcps", gcps);
	if (return_check(ret, "write_array_to_h5", error_head)) return -1;

	/*
	* 根据控制点数据拟合经纬度、下视角与像素坐标（行、列）之间的多项式关系
	*/

	double mean_lon, mean_lat, mean_inc, max_lon, max_lat, max_inc, min_lon, min_lat, min_inc;
	Mat lon, lat, inc, row, col;
	gcps(cv::Range(0, gcps.rows), cv::Range(0, 1)).copyTo(lon);
	gcps(cv::Range(0, gcps.rows), cv::Range(1, 2)).copyTo(lat);
	gcps(cv::Range(0, gcps.rows), cv::Range(3, 4)).copyTo(row);
	gcps(cv::Range(0, gcps.rows), cv::Range(4, 5)).copyTo(col);
	gcps(cv::Range(0, gcps.rows), cv::Range(5, 6)).copyTo(inc);
	mean_lon = cv::mean(lon)[0];
	mean_lat = cv::mean(lat)[0];
	mean_inc = cv::mean(inc)[0];
	cv::minMaxLoc(lon, &min_lon, &max_lon);
	cv::minMaxLoc(lat, &min_lat, &max_lat);
	cv::minMaxLoc(inc, &min_inc, &max_inc);
	lon = (lon - mean_lon) / (max_lon - min_lon + 1e-10);
	lat = (lat - mean_lat) / (max_lat - min_lat + 1e-10);
	inc = (inc - mean_inc) / (max_inc - min_inc + 1e-10);
	row = (row - double(rows) * 0.5) / (double(rows) + 1e-10);
	col = (col - double(cols) * 0.5) / (double(cols) + 1e-10);

	//拟合经度

	Mat A, B, b, temp, coefficient, error, eye, b_t, a, a_t;
	double rms;
	//eye = Mat::zeros(lon.rows, lon.rows, CV_64F);
	//for (int i = 0; i < lon.rows; i++)
	//{
	//	eye.at<double>(i, i) = 1.0;
	//}
	lon.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		/*error = b_t * (eye - a * error * a_t) * b;*/
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lon;
		temp.at<double>(0, 1) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 2) = double(rows) * 0.5;
		temp.at<double>(0, 3) = double(rows) + 1e-10;
		temp.at<double>(0, 4) = double(cols) * 0.5;
		temp.at<double>(0, 5) = double(cols) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "lon_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合纬度

	lat.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lat;
		temp.at<double>(0, 1) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 2) = double(rows) * 0.5;
		temp.at<double>(0, 3) = double(rows) + 1e-10;
		temp.at<double>(0, 4) = double(cols) * 0.5;
		temp.at<double>(0, 5) = double(cols) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "lat_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合下视角

	inc.copyTo(b);
	A = Mat::ones(inc.rows, 6, CV_64F);
	col.copyTo(A(cv::Range(0, inc.rows), cv::Range(1, 2)));
	temp = col.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(2, 3)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(3, 4)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(4, 5)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(5, 6)));
	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 11, CV_64F);
		temp.at<double>(0, 0) = mean_inc;
		temp.at<double>(0, 1) = max_inc - min_inc + 1e-10;
		temp.at<double>(0, 2) = double(cols) * 0.5;
		temp.at<double>(0, 3) = double(cols) + 1e-10;
		temp.at<double>(0, 10) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(4, 10)));
		ret = write_array_to_h5(dst_h5_filename, "inc_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合行坐标

	row.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(rows) * 0.5; 
		temp.at<double>(0, 1) = double(rows) + 1e-10; 
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "row_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合列坐标

	col.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(cols) * 0.5;
		temp.at<double>(0, 1) = double(cols) + 1e-10;
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "col_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	/*
	* 写入轨道数据
	*/

	Mat stateVec;
	ret = xmldoc.XMLFile_load(xml_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	ret = xmldoc.get_stateVec_from_TSX(stateVec);
	if (return_check(ret, "get_stateVec_from_TSX()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "state_vec", stateVec);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	/*
	* 写入多普勒中心频率参数
	*/

	Mat Dc;
	ret = xmldoc.get_dopplerCentroid_from_TSX(Dc);
	if (return_check(ret, "get_dopplerCentroid_from_TSX()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "doppler_centroid", Dc);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	/*
	* 写入其他辅助参数
	*/
	string file_type, sensor, polarization, imaging_mode,
		lookside, orbit_dir, acquisition_start_time, acquisition_stop_time, process_state;
	
	//数据类型
	ret = write_str_to_h5(dst_h5_filename, "file_type", "SLC");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//卫星名称
	ret = xmldoc.get_str_para("mission", sensor);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "sensor", sensor.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//极化
	ret = xmldoc.get_str_para("polLayer", polarization);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "polarization", polarization.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//工作模式
	ret = xmldoc.get_str_para("imagingMode", imaging_mode);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "imaging_mode", imaging_mode.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//视向
	ret = xmldoc.get_str_para("lookDirection", lookside);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "lookside", lookside.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//轨道方向
	ret = xmldoc.get_str_para("orbitDirection", orbit_dir);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "orbit_dir", orbit_dir.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;

	//处理等级
	ret = write_str_to_h5(dst_h5_filename, "process_state", "InSAR_0");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//处理描述
	ret = write_str_to_h5(dst_h5_filename, "comment", "import from TerraSAR-X Single Look Complex, unprocessed.");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;

	double orbit_altitude, carrier_frequency, incidence_center,slant_range_first_pixel,
		slant_range_last_pixel, heading, prf, scene_center_lon, scene_center_lat, scene_topleft_lon,
		scene_topleft_lat, scene_bottomleft_lon, scene_bottomleft_lat, scene_topright_lon, scene_topright_lat,
		scene_bottomright_lon, scene_bottomright_lat, azimuth_resolution,
		range_resolution, azimuth_spacing, range_spacing;
	Mat tmp = Mat::zeros(1, 1, CV_64F);
	//轨道高度,TerraSAR没提供，设置为-1
	tmp.at<double>(0, 0) = -1;
	ret = write_array_to_h5(dst_h5_filename, "orbit_altitude", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	
	TiXmlElement* pnode, * pchild;

	//拍摄起始时间
	ret = xmldoc.find_node("start", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = xmldoc._find_node(pnode, "timeUTC", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_start_time", pchild->GetText());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//拍摄结束时间
	ret = xmldoc.find_node("stop", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = xmldoc._find_node(pnode, "timeUTC", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_stop_time", pchild->GetText());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//载频
	ret = xmldoc.find_node("instrument", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = xmldoc._find_node(pnode, "centerFrequency", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = sscanf(pchild->GetText(), "%lf", &carrier_frequency);
	if (ret != 1)
	{
		fprintf(stderr, "TSX2h5(): carrier frequency not found in %s!\n", xml_filename);
		return -1;
	}
	tmp.at<double>(0, 0) = carrier_frequency;
	ret = write_array_to_h5(dst_h5_filename, "carrier_frequency", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//中心下视角
	ret = xmldoc.get_double_para("incidenceAngle", &incidence_center);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = incidence_center;
	ret = write_array_to_h5(dst_h5_filename, "incidence_center", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//最近斜距
	ret = xmldoc.get_double_para("firstPixel", &slant_range_first_pixel);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = slant_range_first_pixel * 299792458.0 / 2;
	ret = write_array_to_h5(dst_h5_filename, "slant_range_first_pixel", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//最远斜距
	ret = xmldoc.get_double_para("lastPixel", &slant_range_last_pixel);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = slant_range_last_pixel * 299792458.0 / 2;
	ret = write_array_to_h5(dst_h5_filename, "slant_range_last_pixel", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//headingAngle
	ret = xmldoc.get_double_para("headingAngle", &heading);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = heading;
	ret = write_array_to_h5(dst_h5_filename, "heading", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//脉冲重复频率
	ret = xmldoc.get_double_para("commonPRF", &prf);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = prf;
	ret = write_array_to_h5(dst_h5_filename, "prf", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//方位向分辨率
	ret = xmldoc.get_double_para("azimuthResolution", &azimuth_resolution);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = azimuth_resolution;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_resolution", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向分辨率
	ret = xmldoc.get_double_para("slantRangeResolution", &range_resolution);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = range_resolution;
	ret = write_array_to_h5(dst_h5_filename, "range_resolution", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//方位向采样间隔
	ret = xmldoc.find_node("productSpecific", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = xmldoc._find_node(pnode, "projectedSpacingAzimuth", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = sscanf(pchild->GetText(), "%lf", &azimuth_spacing);
	if (ret != 1)
	{
		fprintf(stderr, "TSX2h5(): projectedSpacingAzimuth not found in %s!\n", xml_filename);
		return -1;
	}
	tmp.at<double>(0, 0) = azimuth_spacing;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_spacing", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向采样间隔
	ret = xmldoc._find_node(pnode, "commonRSF", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = sscanf(pchild->GetText(), "%lf", &range_spacing);
	range_spacing = VEL_C / range_spacing / 2.0;
	if (ret != 1)
	{
		fprintf(stderr, "TSX2h5(): groundNear not found in %s!\n", xml_filename);
		return -1;
	}
	tmp.at<double>(0, 0) = range_spacing;
	ret = write_array_to_h5(dst_h5_filename, "range_spacing", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;


	int azimuth_len, range_len;
	Mat tmp_int = Mat::zeros(1, 1, CV_32S);
	//方位向像素点数
	ret = xmldoc.get_int_para("numberOfRows", &azimuth_len);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	tmp_int.at<int>(0, 0) = azimuth_len;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_len", tmp_int);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向像素点数
	ret = xmldoc.get_int_para("numberOfColumns", &range_len);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	tmp_int.at<int>(0, 0) = range_len;
	ret = write_array_to_h5(dst_h5_filename, "range_len", tmp_int);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	return 0;
}

int FormatConversion::TSX2h5(const char* xml_filename, const char* dst_h5_filename)
{
	if (xml_filename == NULL ||
		dst_h5_filename == NULL)
	{
		fprintf(stderr, "TSX2h5(): input check failed!\n");
		return -1;
	}
	string main_xml(xml_filename);
	std::replace(main_xml.begin(), main_xml.end(), '/', '\\');
	string folder;
	if (main_xml.length() > main_xml.rfind("\\") && main_xml.rfind("\\") >= 0)
	{
		folder = main_xml.substr(0, main_xml.rfind("\\"));
	}
	else if(main_xml.length() > main_xml.rfind("/") && main_xml.rfind("/") >= 0)
	{
		folder = main_xml.substr(0, main_xml.rfind("/"));
	}
	else
	{
		fprintf(stderr, "TSX2h5(): invalide file %s!\n", main_xml.c_str());
		return -1;
	}
	string GEOREF = folder + "\\ANNOTATION\\GEOREF.xml";
	string COSAR = folder + "\\IMAGEDATA\\";
	XMLFile xmldoc;
	int ret = xmldoc.XMLFile_load(xml_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	TiXmlElement* pRoot = NULL, * pnode = NULL;
	ret = xmldoc.find_node("imageData", pRoot);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = xmldoc._find_node(pRoot, "filename", pnode);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	COSAR = COSAR + pnode->GetText();
	ret = TSX2h5(COSAR.c_str(), xml_filename, GEOREF.c_str(), dst_h5_filename);
	if (return_check(ret, "TSX2h5()", error_head)) return -1;
	return 0;
}

int FormatConversion::read_POD(const char* POD_filename, double start_time, double stop_time, const char* dst_h5_filename)
{
	if (POD_filename == NULL ||
		dst_h5_filename == NULL ||
		start_time >= stop_time ||
		start_time < 0.0
		)
	{
		fprintf(stderr, "read_POD():input check failed!\n");
		return -1;
	}

	/*
	* 读取精密轨道数据
	*/

	int ret, numOfstateVec;
	TiXmlElement* pnode, * pchild;
	XMLFile xmldoc;
	ret = xmldoc.XMLFile_load(POD_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	ret = xmldoc.find_node("List_of_OSVs", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &numOfstateVec);
	if (ret != 1)
	{
		fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
		return -1;
	}

	Mat tmp = Mat::zeros(numOfstateVec, 7, CV_64F);
	ret = xmldoc.find_node("OSV", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	string str;
	double gps_time, x, y, z, vx, vy, vz;
	bool start = false; bool stop = false;
	int count = 0;
	for (int i = 0; i < numOfstateVec; i++)
	{
		if (!pnode || stop) break;
		ret = xmldoc._find_node(pnode, "UTC", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "read_POD(): node UTC not found!\n");
			return -1;
		}
		str = pchild->GetText();
		str = str.substr(4);
		ret = utc2gps(str.c_str(), &gps_time);
		if (return_check(ret, "utc2gps()", error_head)) return -1;
		if (gps_time <= start_time && fabs(gps_time - start_time) <= 100.0) start = true;
		if (gps_time >= stop_time && fabs(gps_time - stop_time) >= 100.0) stop = true;

		if (start && !stop)//开始记录
		{
			ret = xmldoc._find_node(pnode, "X", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node X not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &x);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}
			ret = xmldoc._find_node(pnode, "Y", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node Y not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &y);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}
			ret = xmldoc._find_node(pnode, "Z", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node Z not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &z);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}
			ret = xmldoc._find_node(pnode, "VX", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node VX not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &vx);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}
			ret = xmldoc._find_node(pnode, "VY", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node VY not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &vy);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}
			ret = xmldoc._find_node(pnode, "VZ", pchild);
			if (ret < 0)
			{
				fprintf(stderr, "read_POD(): node VZ not found!\n");
				return -1;
			}
			ret = sscanf(pchild->GetText(), "%lf", &vz);
			if (ret != 1)
			{
				fprintf(stderr, "read_POD(): %s: unknown data format!\n", POD_filename);
				return -1;
			}

			tmp.at<double>(count, 0) = gps_time;
			tmp.at<double>(count, 1) = x;
			tmp.at<double>(count, 2) = y;
			tmp.at<double>(count, 3) = z;
			tmp.at<double>(count, 4) = vx;
			tmp.at<double>(count, 5) = vy;
			tmp.at<double>(count, 6) = vz;
			count++;
		}

		pnode = pnode->NextSiblingElement();
	}
	Mat stateVec;
	if (count < 1)
	{
		fprintf(stderr, "read_POD(): orbit mismatch! please check if POD file!\n");
		return -1;
	}
	tmp(cv::Range(0, count), cv::Range(0, 7)).copyTo(stateVec);

	/*
	* 写入精密轨道数据
	*/

	ret = write_array_to_h5(dst_h5_filename, "fine_state_vec", stateVec);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	return 0;
}

int FormatConversion::read_slc_from_Sentinel(
	const char* filename, 
	const char* xml_filename,
	ComplexMat& slc,
	Mat& gcps_line
)
{
	if (filename == NULL ||
		xml_filename == NULL
		)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): input check failed!\n");
		return -1;
	}
	XMLFile xmldoc;
	int ret;
	ret = xmldoc.XMLFile_load(xml_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	TiXmlElement* pnode = NULL;
	int linesPerBurst, samplesPerBurst, burst_count;
	/*
	* 读取xml文件里的burst参数
	*/
	ret = xmldoc.get_int_para("linesPerBurst", &linesPerBurst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	ret = xmldoc.get_int_para("samplesPerBurst", &samplesPerBurst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;

	ret = xmldoc.find_node("burstList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &burst_count);
	if (ret != 1)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): %s: unknown data format!\n", xml_filename);
		return -1;
	}



	//逐个burst读取内容
	TiXmlElement* pchild = NULL;
	ret = xmldoc._find_node(pnode, "burst", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): node 'burst' not found!\n");
		return -1;
	}
	int64 bytesoffset;
	FILE* fp = NULL;
	fopen_s(&fp, filename, "rb");
	if (!fp)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): failed to open %s!\n", filename);
		return -1;
	}
	Mat gcps_merged_line_num = Mat::zeros(1, burst_count + 1, CV_32S);
	ComplexMat last_burst, this_burst;
	ret = get_a_burst(pchild, xmldoc, fp, linesPerBurst, samplesPerBurst, slc);
	if (return_check(ret, "get_a_burst()", error_head)) return -1;
	gcps_merged_line_num.at<int>(0, 1) = slc.GetRows();
	pchild = pchild->NextSiblingElement();
	int overlapSize;
	for (int i = 1; i < burst_count; i++)
	{
		if (!pchild) break;
		ret = get_a_burst(pchild, xmldoc, fp, linesPerBurst, samplesPerBurst, this_burst);
		if (return_check(ret, "get_a_burst()", error_head)) return -1;
		ret = deburst_overlapSize(slc, this_burst, &overlapSize);
		if (return_check(ret, "deburst_overlapSize()", error_head)) return -1;
		ret = burst_stitch(this_burst, slc, overlapSize);
		if (return_check(ret, "burst_stitch()", error_head)) return -1;
		gcps_merged_line_num.at<int>(0, i + 1) = slc.GetRows();
		pchild = pchild->NextSiblingElement();
	}
	if (fp)fclose(fp);
	gcps_merged_line_num.copyTo(gcps_line);
	return 0;
}

int FormatConversion::sentinel_deburst(const char* xml_filename, ComplexMat& slc, Mat& Sentinel)
{
	if (xml_filename == NULL ||
		slc.isempty() ||
		slc.type() != CV_16S
		)
	{
		fprintf(stderr, "sentinel_deburst(): input check failed!\n");
		return -1;
	}
	XMLFile xmldoc;
	int ret;
	ret = xmldoc.XMLFile_load(xml_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	TiXmlElement* pnode = NULL;
	int linesPerBurst, burst_count, invalideLines;
	ret = xmldoc.get_int_para("linesPerBurst", &linesPerBurst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	ret = xmldoc.find_node("burstList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &burst_count);
	if (ret != 1)
	{
		fprintf(stderr, "sentinel_deburst(): %s: unknown data format!\n", xml_filename);
		return -1;
	}
	if (slc.GetRows() != (burst_count * linesPerBurst))
	{
		fprintf(stderr, "sentinel_deburst(): %s: input slc size mismatch!\n", xml_filename);
		return -1;
	}
	TiXmlElement* pchild = NULL;
	ret = xmldoc._find_node(pnode, "burst", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "sentinel_deburst(): node 'burst' not found!\n");
		return -1;
	}

	/*
	* 找到所有无效行
	*/
	long firstValidSample;
	char* ptr;
	const char* p;
	invalideLines = 0;
	int count = 0;
	Mat sentinel = Mat::zeros(1, slc.GetRows(), CV_64F);
	for (int i = 0; i < burst_count; i++)
	{
		if (!pchild) break;
		ret = xmldoc._find_node(pchild, "firstValidSample", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "sentinel_deburst(): node 'firstValidSample' not found!\n");
			return -1;
		}
		p = pnode->GetText();
		firstValidSample = strtol(p, &ptr, 0);
		if (firstValidSample < 0)
		{
			invalideLines++;
			sentinel.at<double>(0, count) = -1;
		}
		count++;
		for (int j = 0; j < linesPerBurst - 1; j++)
		{
			firstValidSample = strtol(ptr, &ptr, 0);
			if (firstValidSample < 0)
			{
				invalideLines++;
				sentinel.at<double>(0, count) = -1;
			}
			count++;
		}
		pchild = pchild->NextSiblingElement();
	}
	count = 0;
	ComplexMat tmp; Mat sentinel_accu = Mat::zeros(slc.GetRows() - invalideLines, 1, CV_32S);
	tmp.re.create(slc.GetRows() - invalideLines, slc.GetCols(), CV_16S);
	tmp.im.create(slc.GetRows() - invalideLines, slc.GetCols(), CV_16S);
	int samplesPerLine = slc.GetCols();
	ComplexMat c;
	int invalideLine_accu = 0;
	for (int i = 0; i < slc.GetRows(); i++)
	{
		if (sentinel.at<double>(0, i) > -0.5)
		{
			sentinel_accu.at<int>(count, 0) = invalideLine_accu;
			c = slc(cv::Range(i, i + 1), cv::Range(0, samplesPerLine));
			tmp.SetValue(cv::Range(count, count + 1), cv::Range(0, samplesPerLine), c);
			count++;
		}
		else
		{
			invalideLine_accu++;
		}
	}
	slc = tmp;
	Mat t = Mat::zeros(burst_count, 1, CV_32S);
	count = 0;
	for (int i = 0; i < sentinel_accu.rows - 1; i++)
	{
		if (sentinel_accu.at<int>(i, 0) != sentinel_accu.at<int>(i + 1, 0))
		{
			t.at<int>(count, 0) = sentinel_accu.at<int>(i + 1, 0);
			count++;
		}
	}
	if (invalideLine_accu != sentinel_accu.at<int>(sentinel_accu.rows - 1, 0))
	{
		t.at<int>(count, 0) = invalideLine_accu;
	}
	else
	{
		t.at<int>(count, 0) = t.at<int>(count - 1, 0);
	}
	t.copyTo(Sentinel);
	return 0;
}

int FormatConversion::sentinel2h5(const char* tiff_filename, const char* xml_filename, const char* dst_h5_filename, const char* POD_file)
{
	if (tiff_filename == NULL ||
		xml_filename == NULL ||
		dst_h5_filename == NULL
		)
	{
		fprintf(stderr, "sentinel2h5(): input check failed!\n");
		return -1;
	}

	/*
	* 检查h5文件是否已经存在
	*/

	int ret;
	ret = creat_new_h5(dst_h5_filename);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;

	/*
	* 写入slc数据
	*/

	ComplexMat slc;Mat gcps_line_index;
	int rows, cols;
	ret = read_slc_from_Sentinel(tiff_filename, xml_filename, slc, gcps_line_index);//需要deburst
	if (return_check(ret, "read_slc_from_Sentinel()", error_head)) return -1;
	//ret = sentinel_deburst(xml_filename, slc, sentinel);
	//if (return_check(ret, "sentinel_deburst()", error_head)) return -1;
	rows = slc.GetRows(); cols = slc.GetCols();
	ret = write_array_to_h5(dst_h5_filename, "s_re", slc.re);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "s_im", slc.im);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	slc.re.release();
	slc.im.release();

	/*
	* 写入控制点数据
	*/

	Mat gcps;
	XMLFile xmldoc;
	int linesPerburst;
	ret = xmldoc.XMLFile_load(xml_filename);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	ret = xmldoc.get_gcps_from_sentinel(gcps);
	if (return_check(ret, "get_gcps_from_sentinel()", error_head)) return -1;
	ret = xmldoc.get_int_para("linesPerBurst", &linesPerburst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	for (int i = 0; i < gcps.rows; i++)
	{
		int temp_r, xx;
		temp_r = (int)gcps.at<double>(i, 3);
		xx = temp_r / linesPerburst;
		if (xx > 0 && xx < gcps_line_index.cols)
		{
			gcps.at<double>(i, 3) = (double)gcps_line_index.at<int>(0, xx);
		}
	}
	ret = write_array_to_h5(dst_h5_filename, "gcps", gcps);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	/*
	* 根据控制点数据拟合经纬度、下视角与像素坐标（行、列）之间的多项式关系
	*/

	double mean_lon, mean_lat, mean_inc, max_lon, max_lat, max_inc, min_lon, min_lat, min_inc;
	Mat lon, lat, inc, row, col;
	gcps(cv::Range(0, gcps.rows), cv::Range(0, 1)).copyTo(lon);
	gcps(cv::Range(0, gcps.rows), cv::Range(1, 2)).copyTo(lat);
	gcps(cv::Range(0, gcps.rows), cv::Range(3, 4)).copyTo(row);
	gcps(cv::Range(0, gcps.rows), cv::Range(4, 5)).copyTo(col);
	gcps(cv::Range(0, gcps.rows), cv::Range(5, 6)).copyTo(inc);
	mean_lon = cv::mean(lon)[0];
	mean_lat = cv::mean(lat)[0];
	mean_inc = cv::mean(inc)[0];
	cv::minMaxLoc(lon, &min_lon, &max_lon);
	cv::minMaxLoc(lat, &min_lat, &max_lat);
	cv::minMaxLoc(inc, &min_inc, &max_inc);
	lon = (lon - mean_lon) / (max_lon - min_lon + 1e-10);
	lat = (lat - mean_lat) / (max_lat - min_lat + 1e-10);
	inc = (inc - mean_inc) / (max_inc - min_inc + 1e-10);
	row = (row + 1 - double(rows) * 0.5) / (double(rows) + 1e-10);//sentinel行列起点为0，+1统一为1.
	col = (col + 1 - double(cols) * 0.5) / (double(cols) + 1e-10);
	
	//拟合经度

	Mat A, B, b, temp, coefficient, error, eye, b_t, a, a_t;
	double rms;
	/*eye = Mat::zeros(lon.rows, lon.rows, CV_64F);
	for (int i = 0; i < lon.rows; i++)
	{
		eye.at<double>(i, i) = 1.0;
	}*/
	lon.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lon;
		temp.at<double>(0, 1) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 2) = double(rows) * 0.5;
		temp.at<double>(0, 3) = double(rows) + 1e-10;
		temp.at<double>(0, 4) = double(cols) * 0.5;
		temp.at<double>(0, 5) = double(cols) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "lon_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合纬度

	lat.copyTo(b);

	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lat;
		temp.at<double>(0, 1) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 2) = double(rows) * 0.5;
		temp.at<double>(0, 3) = double(rows) + 1e-10;
		temp.at<double>(0, 4) = double(cols) * 0.5;
		temp.at<double>(0, 5) = double(cols) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "lat_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合下视角

	inc.copyTo(b);
	A = Mat::ones(inc.rows, 6, CV_64F);
	col.copyTo(A(cv::Range(0, inc.rows), cv::Range(1, 2)));
	temp = col.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(2, 3)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(3, 4)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(4, 5)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(5, 6)));
	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 11, CV_64F);
		temp.at<double>(0, 0) = mean_inc;
		temp.at<double>(0, 1) = max_inc - min_inc + 1e-10;
		temp.at<double>(0, 2) = double(cols) * 0.5;
		temp.at<double>(0, 3) = double(cols) + 1e-10;
		temp.at<double>(0, 10) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(4, 10)));
		ret = write_array_to_h5(dst_h5_filename, "inc_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合行坐标

	row.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(rows) * 0.5; 
		temp.at<double>(0, 1) = double(rows) + 1e-10; 
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "row_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	//拟合列坐标

	col.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(cols) * 0.5;
		temp.at<double>(0, 1) = double(cols) + 1e-10;
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		ret = write_array_to_h5(dst_h5_filename, "col_coefficient", temp);
		if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	}

	/*
	* 写入轨道数据
	*/

	Mat stateVec;
	ret = xmldoc.get_stateVec_from_sentinel(stateVec);
	if (return_check(ret, "get_stateVec_from_sentinel()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "state_vec", stateVec);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	

	/*
	* 写入多普勒中心频率数据
	*/

	Mat Dc;
	ret = xmldoc.get_dopplerCentroid_from_sentinel(Dc);
	if (return_check(ret, "get_dopplerCentroid_from_sentinel()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5_filename, "doppler_centroid", Dc);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	/*
	* 其他辅助数据
	*/

	string file_type, sensor, polarization, imaging_mode,
		lookside, orbit_dir, acquisition_start_time, acquisition_stop_time, process_state;

	//数据类型
	ret = xmldoc.get_str_para("productType", file_type);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "file_type", file_type.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//卫星名称
	ret = xmldoc.get_str_para("missionId", sensor);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "sensor", sensor.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//极化
	ret = xmldoc.get_str_para("polarisation", polarization);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "polarization", polarization.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//工作模式
	ret = xmldoc.get_str_para("mode", imaging_mode);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "imaging_mode", imaging_mode.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//视向
	ret = write_str_to_h5(dst_h5_filename, "lookside", "RIGHT");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//轨道方向
	ret = xmldoc.get_str_para("pass", orbit_dir);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "orbit_dir", orbit_dir.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//拍摄起始时间
	ret = xmldoc.get_str_para("startTime", acquisition_start_time);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_start_time", acquisition_start_time.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//拍摄结束时间
	ret = xmldoc.get_str_para("stopTime", acquisition_stop_time);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_stop_time", acquisition_stop_time.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;

	//精密轨道数据
	if (POD_file)
	{
		double start_t, end_t;
		utc2gps(acquisition_start_time.c_str(), &start_t);
		utc2gps(acquisition_stop_time.c_str(), &end_t);
		read_POD(POD_file, start_t, end_t, dst_h5_filename);
	}

	//处理等级
	ret = write_str_to_h5(dst_h5_filename, "process_state", "InSAR_0");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//处理描述
	ret = write_str_to_h5(dst_h5_filename, "comment", "import from sentinel1 Single Look Complex, unprocessed.");
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;

	double orbit_altitude, carrier_frequency, incidence_center, slant_range_first_pixel,
		slant_range_last_pixel, heading, prf, scene_center_lon, scene_center_lat, scene_topleft_lon,
		scene_topleft_lat, scene_bottomleft_lon, scene_bottomleft_lat, scene_topright_lon, scene_topright_lat,
		scene_bottomright_lon, scene_bottomright_lat, azimuth_resolution,
		range_resolution, azimuth_spacing, range_spacing;
	Mat tmp = Mat::zeros(1, 1, CV_64F);

	//轨道高度,sentinel没提供，设置为-1
	tmp.at<double>(0, 0) = -1;
	ret = write_array_to_h5(dst_h5_filename, "orbit_altitude", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//载频
	ret = xmldoc.get_double_para("radarFrequency", &carrier_frequency);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = carrier_frequency;
	ret = write_array_to_h5(dst_h5_filename, "carrier_frequency", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//中心下视角
	ret = xmldoc.get_double_para("incidenceAngleMidSwath", &incidence_center);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = incidence_center;
	ret = write_array_to_h5(dst_h5_filename, "incidence_center", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//最近斜距
	ret = xmldoc.get_double_para("slantRangeTime", &slant_range_first_pixel);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = slant_range_first_pixel;
	ret = write_array_to_h5(dst_h5_filename, "slant_range_first_pixel", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//最远斜距，sentinel未提供，设置为-1
	tmp.at<double>(0, 0) = -1;
	ret = write_array_to_h5(dst_h5_filename, "slant_range_last_pixel", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//heading
	ret = xmldoc.get_double_para("platformHeading", &heading);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = heading;
	ret = write_array_to_h5(dst_h5_filename, "heading", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//脉冲重复频率
	ret = xmldoc.get_double_para("prf", &prf);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = prf;
	ret = write_array_to_h5(dst_h5_filename, "prf", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//方位向分辨率
	tmp.at<double>(0, 0) = 20;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_resolution", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向分辨率
	tmp.at<double>(0, 0) = 5;
	ret = write_array_to_h5(dst_h5_filename, "range_resolution", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//方位向采样间隔
	ret = xmldoc.get_double_para("azimuthPixelSpacing", &azimuth_spacing);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = azimuth_spacing;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_spacing", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向采样间隔
	ret = xmldoc.get_double_para("rangePixelSpacing", &range_spacing);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	tmp.at<double>(0, 0) = range_spacing;
	ret = write_array_to_h5(dst_h5_filename, "range_spacing", tmp);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;


	int azimuth_len, range_len;
	Mat tmp_int = Mat::zeros(1, 1, CV_32S);
	//方位向像素点数
	tmp_int.at<int>(0, 0) = rows;
	ret = write_array_to_h5(dst_h5_filename, "azimuth_len", tmp_int);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	//距离向像素点数
	tmp_int.at<int>(0, 0) = cols;
	ret = write_array_to_h5(dst_h5_filename, "range_len", tmp_int);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	return 0;
}

int FormatConversion::import_sentinel(
	const char* manifest,
	const char* subswath_name,
	const char* polarization,
	const char* dest_h5_file,
	const char* PODFile
)
{
	if (manifest == NULL ||
		subswath_name == NULL ||
		polarization == NULL ||
		dest_h5_file == NULL)
	{
		fprintf(stderr, "import_sentinel(): input check failed!\n");
		return -1;
	}
	int ret;
	string xmlhead, tiffhead, subswath, polar;
	if (0 == strcmp("iw2", subswath_name))subswath = "iw2";
	else if (0 == strcmp("iw3", subswath_name))subswath = "iw3";
	else subswath = "iw1";

	if (0 == strcmp("vh", polarization)) polar = "vh";
	else polar = "vv";

	XMLFile xmldoc;
	ret = xmldoc.XMLFile_load(manifest);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	TiXmlElement* root = NULL, * pnode = NULL;
	ret = xmldoc.find_node("dataObjectSection", root);
	if (return_check(ret, "find_node()", error_head)) return -1;

	ret = xmldoc._find_node(root, "dataObject", pnode);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	string tmp(pnode->FirstAttribute()->Value());
	if (0 == strcmp(tmp.substr(0, 10).c_str(), "products1a"))//sentinel-1A
	{
		xmlhead = "products1a" + subswath + "slc" + polar;
		tiffhead = "s1a" + subswath + "slc" + polar;
	}
	else//sentinel-1B
	{
		xmlhead = "products1b" + subswath + "slc" + polar;
		tiffhead = "s1b" + subswath + "slc" + polar;
	}

	string tiff_filename, xml_filename;
	string main_xml(manifest);
	string folder;
	if (main_xml.length() > main_xml.rfind("\\") && main_xml.rfind("\\") >= 0)
	{
		folder = main_xml.substr(0, main_xml.rfind("\\"));
	}
	else if (main_xml.length() > main_xml.rfind("/") && main_xml.rfind("/") >= 0)
	{
		folder = main_xml.substr(0, main_xml.rfind("/"));
	}
	else
	{
		fprintf(stderr, "import_sentinel(): invalid manifest file %s !\n", manifest);
		return -1;
	}
	TiXmlElement* pchild = NULL;

	while (pnode)
	{
		tmp = pnode->FirstAttribute()->Value();
		if (tmp.length() > 18)
		{
			if (0 == strcmp(tmp.substr(0, 18).c_str(), xmlhead.c_str()))
			{
				ret = xmldoc._find_node(pnode, "fileLocation", pchild);
				if (return_check(ret, "_find_node()", error_head)) return -1;
				xml_filename = pchild->Attribute("href");
				xml_filename = folder + xml_filename.substr(1, xml_filename.length() - 1);
			}
		}
		if (tmp.length() > 11)
		{
			if (0 == strcmp(tmp.substr(0, 11).c_str(), tiffhead.c_str()))
			{
				ret = xmldoc._find_node(pnode, "fileLocation", pchild);
				if (return_check(ret, "_find_node()", error_head)) return -1;
				tiff_filename = pchild->Attribute("href");
				tiff_filename = folder + tiff_filename.substr(1, tiff_filename.length() - 1);
			}
		}
		pnode = pnode->NextSiblingElement();
	}
	Sentinel1Reader reader(xml_filename.c_str(), tiff_filename.c_str(), PODFile);
	ret = reader.writeToh5(dest_h5_file);
	if (return_check(ret, "writeToh5()", error_head)) return -1;
	return 0;
}

int FormatConversion::get_a_burst(
	TiXmlElement* pnode,
	XMLFile& xmldoc,
	FILE*& fp,
	int linesPerBurst,
	int samplesPerBurst, 
	ComplexMat& burst
)
{
	if (pnode == NULL ||
		linesPerBurst < 1 ||
		samplesPerBurst < 1 ||
		fp == NULL
		)
	{
		fprintf(stderr, "get_a_burst(): input check failed!\n");
		if (fp) {
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	int ret;
	TiXmlElement* pchild = NULL;
	size_t bytesoffset;
	INT16* buf = NULL;
	buf = (INT16*)malloc(linesPerBurst * samplesPerBurst * 2 * sizeof(INT16));
	ret = xmldoc._find_node(pnode, "byteOffset", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "get_a_burst(): node 'byteOffset' not found!\n");
		if (buf) free(buf);
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	ret = sscanf(pchild->GetText(), "%lld", &bytesoffset);
	if (ret != 1)
	{
		fprintf(stderr, "get_a_burst(): unknown data format!\n");
		if (buf) free(buf);
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	fseek(fp, bytesoffset, SEEK_SET);
	fread(buf, sizeof(INT16), linesPerBurst * samplesPerBurst * 2, fp);
	size_t offset = 0;
	burst.re.create(linesPerBurst, samplesPerBurst, CV_16S);
	burst.im.create(linesPerBurst, samplesPerBurst, CV_16S);
	for (int j = 0; j < linesPerBurst; j++)
	{
		for (int k = 0; k < samplesPerBurst; k++)
		{
			burst.re.at<short>(j, k) = buf[offset];
			offset++;
			burst.im.at<short>(j, k) = buf[offset];
			offset++;
		}
	}
	if (buf) free(buf);

	//剔除无效数据
	ret = xmldoc._find_node(pnode, "firstValidSample", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "get_a_burst(): node 'firstValidSample' not found!\n");
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	long firstValidSample;
	char* ptr;
	const char* p;
	int invalideLines = 0;
	int count = 0;
	Mat sentinel = Mat::ones(1, linesPerBurst, CV_64F);
	p = pchild->GetText();
	firstValidSample = strtol(p, &ptr, 0);
	if (firstValidSample < 0)
	{
		invalideLines++;
		sentinel.at<double>(0, count) = -1;
	}
	count++;
	for (int j = 0; j < linesPerBurst - 1; j++)
	{
		firstValidSample = strtol(ptr, &ptr, 0);
		if (firstValidSample < 0)
		{
			invalideLines++;
			sentinel.at<double>(0, count) = -1;
		}
		count++;
	}
	int start, end;
	if (sentinel.at<double>(0, 0) > 0.0)
	{
		start = 0;
	}
	else
	{
		for (int i = 1; i < linesPerBurst; i++)
		{
			if (sentinel.at<double>(0, i - 1) * sentinel.at<double>(0, i) < 0.0)
			{
				start = i; break;
			}
		}
	}
	end = linesPerBurst - (invalideLines - start);
	burst = burst(cv::Range(start, end), cv::Range(0, samplesPerBurst));
	return 0;
}

int FormatConversion::get_burst_sentinel(
	int burst_num,
	const char* xml_file,
	const char* tiff_file,
	ComplexMat& burst,
	int* overlapSize
)
{
	if (burst_num < 1 ||
		xml_file == NULL ||
		tiff_file == NULL ||
		overlapSize == NULL)
	{
		fprintf(stderr, "get_burst_sentinel():  input check failed!\n");
		return -1;
	}


	XMLFile xmldoc;
	int ret;
	ret = xmldoc.XMLFile_load(xml_file);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	TiXmlElement* pnode = NULL;
	int linesPerBurst, samplesPerBurst, burst_count;
	double azimuthTimeInterval;
	/*
	* 读取xml文件里的burst参数
	*/
	ret = xmldoc.get_int_para("linesPerBurst", &linesPerBurst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	ret = xmldoc.get_int_para("samplesPerBurst", &samplesPerBurst);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	ret = xmldoc.get_double_para("azimuthTimeInterval", &azimuthTimeInterval);
	if (return_check(ret, "get_double_para()", error_head)) return -1;

	ret = xmldoc.find_node("burstList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &burst_count);
	if (ret != 1)
	{
		fprintf(stderr, "get_burst_sentinel(): %s: unknown data format!\n", xml_file);
		return -1;
	}
	if (burst_num > burst_count)
	{
		fprintf(stderr, "get_burst_sentinel(): burst_num out of range!\n");
		return -1;
	}


	//读取burst内容
	TiXmlElement* pchild = NULL, * pchild2 = NULL;
	int count = 1;
	ret = xmldoc._find_node(pnode, "burst", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "get_burst_sentinel(): node 'burst' not found!\n");
		return -1;
	}
	while (pchild && count != burst_num)
	{
		pchild = pchild->NextSiblingElement();
		count++;
	}
	

	size_t bytesoffset;
	double azimuthAnxTime, azimuthAnxTime2;
	ret = xmldoc._find_node(pchild, "azimuthAnxTime", pchild2);
	if (ret < 0)
	{
		fprintf(stderr, "get_burst_sentinel(): node 'azimuthAnxTime' not found!\n");
		return -1;
	}
	ret = sscanf(pchild2->GetText(), "%lf", &azimuthAnxTime);
	if (ret != 1)
	{
		fprintf(stderr, "get_burst_sentinel(): unknown data format!\n");
		return -1;
	}
	FILE* fp = NULL;
	fopen_s(&fp, tiff_file, "rb");
	if (!fp)
	{
		fprintf(stderr, "get_burst_sentinel(): failed to open %s!\n", tiff_file);
		return -1;
	}

	short* buf = NULL;
	buf = (short*)malloc(linesPerBurst * samplesPerBurst * 2 * sizeof(short));
	if (!buf)
	{
		fprintf(stderr, "get_burst_sentinel(): out of memory!\n");
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	ret = xmldoc._find_node(pchild, "byteOffset", pchild2);
	if (ret < 0)
	{
		fprintf(stderr, "get_burst_sentinel(): node 'byteOffset' not found!\n");
		if (buf) free(buf);
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	ret = sscanf(pchild2->GetText(), "%lld", &bytesoffset);
	if (ret != 1)
	{
		fprintf(stderr, "get_burst_sentinel(): unknown data format!\n");
		if (buf) free(buf);
		if (fp)
		{
			fclose(fp); fp = NULL;
		}
		return -1;
	}
	fseek(fp, bytesoffset, SEEK_SET);
	fread(buf, sizeof(short), linesPerBurst * samplesPerBurst * 2, fp);
	if (fp)
	{
		fclose(fp); fp = NULL;
	}
	size_t offset = 0;
	burst.re.create(linesPerBurst, samplesPerBurst, CV_16S);
	burst.im.create(linesPerBurst, samplesPerBurst, CV_16S);
	for (int j = 0; j < linesPerBurst; j++)
	{
		for (int k = 0; k < samplesPerBurst; k++)
		{
			burst.re.at<short>(j, k) = buf[offset];
			offset++;
			burst.im.at<short>(j, k) = buf[offset];
			offset++;
		}
	}
	if (buf)
	{
		free(buf); buf = NULL;
	}

	//剔除无效数据
	ret = xmldoc._find_node(pchild, "firstValidSample", pchild2);
	if (ret < 0)
	{
		fprintf(stderr, "get_burst_sentinel(): node 'firstValidSample' not found!\n");
		return -1;
	}
	long firstValidSample;
	char* ptr;
	const char* p;
	int invalideLines = 0;
	count = 0;
	Mat sentinel = Mat::ones(1, linesPerBurst, CV_64F);
	p = pchild2->GetText();
	firstValidSample = strtol(p, &ptr, 0);
	if (firstValidSample < 0)
	{
		invalideLines++;
		sentinel.at<double>(0, count) = -1;
	}
	count++;
	for (int j = 0; j < linesPerBurst - 1; j++)
	{
		firstValidSample = strtol(ptr, &ptr, 0);
		if (firstValidSample < 0)
		{
			invalideLines++;
			sentinel.at<double>(0, count) = -1;
		}
		count++;
	}
	int start, end;
	if (sentinel.at<double>(0, 0) > 0.0)
	{
		start = 0;
	}
	else
	{
		for (int i = 1; i < linesPerBurst; i++)
		{
			if (sentinel.at<double>(0, i - 1) * sentinel.at<double>(0, i) < 0.0)
			{
				start = i; break;
			}
		}
	}
	end = linesPerBurst - (invalideLines - start);
	burst = burst(cv::Range(start, end), cv::Range(0, samplesPerBurst));

	//求取overlapSize
	if (!pchild->NextSiblingElement()) *overlapSize = -1;//-1表示最后一个burst
	else
	{
		pchild = pchild->NextSiblingElement();
		ret = xmldoc._find_node(pchild, "azimuthAnxTime", pchild2);
		if (ret < 0)
		{
			fprintf(stderr, "get_burst_sentinel(): node 'azimuthAnxTime' not found!\n");
			return -1;
		}
		ret = sscanf(pchild2->GetText(), "%lf", &azimuthAnxTime2);
		if (ret != 1)
		{
			fprintf(stderr, "get_burst_sentinel(): unknown data format!\n");
			return -1;
		}
		int temp = std::round((azimuthAnxTime2 - azimuthAnxTime) / azimuthTimeInterval);
		*overlapSize = linesPerBurst - invalideLines - temp;
	}
	return 0;
}

int FormatConversion::deburst_overlapSize(ComplexMat& last_burst, ComplexMat& this_burst, int* overlapSize)
{
	if (last_burst.isempty() ||
		this_burst.isempty() ||
		last_burst.GetCols() != this_burst.GetCols() ||
		overlapSize == NULL
		)
	{
		fprintf(stderr, "deburst_overlapSize(): input check failed!\n");
		return -1;
	}

	int nr, nc, nr2, nc2, match_wnd_rows, match_wnd_cols, col_start, col_end, offset_row, offset_col, ret;
	ComplexMat match_wnd, match_wnd2;
	nr = last_burst.GetRows(); nc = last_burst.GetCols(); nr2 = this_burst.GetRows(); nc2 = this_burst.GetCols();
	match_wnd_cols = 1000;
	match_wnd_rows = (nr < nr2 ? nr : nr2) / 3;
	col_start = nc / 2 - match_wnd_cols; 
	col_start = col_start < 0 ? 0 : col_start;
	col_end = col_start + match_wnd_cols;
	col_end = col_end > nc ? nc : col_end;
	match_wnd_cols = col_end - col_start;
	//match_wnd_rows = 600;
	match_wnd = last_burst(cv::Range(nr - match_wnd_rows, nr), cv::Range(0, nc));
	match_wnd2 = this_burst(cv::Range(0, match_wnd_rows), cv::Range(0, nc));
	match_wnd.convertTo(match_wnd, CV_64F);
	match_wnd2.convertTo(match_wnd2, CV_64F);
	

	//求取偏移量
	Registration regis;
	//ComplexMat t1, t2;
	//regis.interp_paddingzero(match_wnd, t1, 8);
	//regis.interp_paddingzero(match_wnd2, t2, 8);
	//Utils util;
	//util.saveSLC("E:/working_dir/projects/software/InSAR/bin/match_wnd.jpg", 65, last_burst);
	//util.saveSLC("E:/working_dir/projects/software/InSAR/bin/match_wnd2.jpg", 65, this_burst);
	ret = regis.real_coherent(match_wnd, match_wnd2, &offset_row, &offset_col);
	if (return_check(ret, "real_coherent()", error_head)) return -1;

	//水平搬移this_burst
	
	//if (offset_col > 0)
	//{
	//	ComplexMat tmp_burst, tmp_burst2;
	//	tmp_burst.re = Mat::zeros(this_burst.GetRows(), this_burst.GetCols(), CV_16S);
	//	tmp_burst.im = Mat::zeros(this_burst.GetRows(), this_burst.GetCols(), CV_16S);
	//	tmp_burst2 = this_burst(cv::Range(0, nr2), cv::Range(0, nc - offset_col));
	//	tmp_burst.SetValue(cv::Range(0, nr2), cv::Range(offset_col, nc), tmp_burst2);
	//	this_burst = tmp_burst;
	//}
	//if (offset_col < 0)
	//{
	//	ComplexMat tmp_burst, tmp_burst2;
	//	tmp_burst.re = Mat::zeros(this_burst.GetRows(), this_burst.GetCols(), CV_16S);
	//	tmp_burst.im = Mat::zeros(this_burst.GetRows(), this_burst.GetCols(), CV_16S);
	//	tmp_burst2 = this_burst(cv::Range(0, nr2), cv::Range(-offset_col, nc));
	//	tmp_burst.SetValue(cv::Range(0, nr2), cv::Range(0, nc + offset_col), tmp_burst2);
	//	this_burst = tmp_burst;
	//}

	*overlapSize = (offset_row > 0 ? offset_row : -offset_row);

	return 0;
}

int FormatConversion::burst_stitch(
	ComplexMat& src_burst,
	ComplexMat& dst_burst,
	int overlapSize,
	const char* stitch_type
)
{
	if (src_burst.isempty() ||
		dst_burst.isempty() ||
		src_burst.GetCols() != dst_burst.GetCols() ||
		overlapSize < 0||
		overlapSize > dst_burst.GetRows()||
		overlapSize > src_burst.GetRows() ||
		src_burst.type() != CV_16S||
		dst_burst.type() != CV_16S
		)
	{
		fprintf(stderr, "burst_stitch(): input check failed!\n");
		return -1;
	}

	//stitch
	
	int nr, nc, nr2, nc2;
	nr = dst_burst.GetRows(); nc = dst_burst.GetCols(); nr2 = src_burst.GetRows(); nc2 = nc;
	Mat src_lowerpart_real, src_lowerpart_imag;
	if (strcmp(stitch_type, "low") == 0)
	{
		if (overlapSize < nr2)
		{
			src_burst.re(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_real);
			src_burst.im(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_imag);
			cv::vconcat(dst_burst.re, src_lowerpart_real, dst_burst.re);
			cv::vconcat(dst_burst.im, src_lowerpart_imag, dst_burst.im);
		}
	}
	else if (strcmp(stitch_type, "mid") == 0)
	{
		int mid = overlapSize / 2;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < mid; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				dst_burst.re.at<short>(nr - mid + i, j) = src_burst.re.at<short>(i + overlapSize - mid, j);
				dst_burst.im.at<short>(nr - mid + i, j) = src_burst.im.at<short>(i + overlapSize - mid, j);
			}
		}
		if (mid < nr2)
		{
			src_burst.re(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_real);
			src_burst.im(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_imag);
			cv::vconcat(dst_burst.re, src_lowerpart_real, dst_burst.re);
			cv::vconcat(dst_burst.im, src_lowerpart_imag, dst_burst.im);
		}

	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < overlapSize; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				dst_burst.re.at<short>(nr - overlapSize + i, j) = src_burst.re.at<short>(i, j);
				dst_burst.im.at<short>(nr - overlapSize + i, j) = src_burst.im.at<short>(i, j);
			}
		}
		if (overlapSize < nr2)
		{
			src_burst.re(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_real);
			src_burst.im(cv::Range(overlapSize, nr2), cv::Range(0, nc)).copyTo(src_lowerpart_imag);
			cv::vconcat(dst_burst.re, src_lowerpart_real, dst_burst.re);
			cv::vconcat(dst_burst.im, src_lowerpart_imag, dst_burst.im);
		}

	}
	return 0;
}

int FormatConversion::read_slc_from_ALOS(const char* img_file, ComplexMat& slc)
{
	if (img_file == NULL)
	{
		fprintf(stderr, "read_slc_from_ALOS(): input check failed!\n");
		return -1;
	}
	FILE* fp = NULL;
	fopen_s(&fp, img_file, "rb");
	char buf[2048];
	memset(buf, 0, 2048);
	if (!fp)
	{
		fprintf(stderr, "read_slc_from_ALOS(): failed to open %s!\n", img_file);
		return -1;
	}
	int rows, cols, sarfd_record_length, record_length, sardata_offset;
	//检查数据格式
	fseek(fp, 9 - 1, SEEK_SET);
	if (fread(&sarfd_record_length, 4, 1, fp) != 1)
	{
		fprintf(stderr, "read_slc_from_ALOS(): %s: unknown format!\n", img_file);
		if (fp) fclose(fp);
		return -1;
	}
	sarfd_record_length = Big2Little32(sarfd_record_length);
	if (sarfd_record_length != 720)
	{
		fprintf(stderr, "read_slc_from_ALOS(): %s: unknown format!\n", img_file);
		if (fp) fclose(fp);
		return -1;
	}
	char* ptr = NULL;
	//确定slc数据尺寸
	fseek(fp, 237 - 1, SEEK_SET);
	fread(buf, 1, 8, fp);
	rows = strtol(buf, &ptr, 0);
	fseek(fp, sarfd_record_length + 9 - 1, SEEK_SET);
	fread(&record_length, 4, 1, fp);
	record_length = Big2Little32(record_length);
	fseek(fp, sarfd_record_length + 25 - 1, SEEK_SET);
	fread(&cols, 4, 1, fp);
	cols = Big2Little32(cols);
	//读取数据偏移量
	memset(buf, 0, 2048);
	fseek(fp, 277 - 1, SEEK_SET);
	fread(buf, 1, 4, fp);
	sardata_offset = strtol(buf, &ptr, 0);
	//读取数据
	slc.re.create(rows, cols, CV_32F);
	slc.im.create(rows, cols, CV_32F);
	float* p = (float*)malloc(record_length);
	if (!p)
	{
		fprintf(stderr, "read_slc_from_ALOS(): out of memory!\n");
		if (fp) fclose(fp);
		return -1;
	}
	fseek(fp, sarfd_record_length + sardata_offset, SEEK_SET);
	for (int i = 0; i < rows; i++)
	{
		fread(p, 1, record_length, fp);
		for (int j = 0; j < cols; j++)
		{
			slc.re.at<float>(i, j) = ReverseFloat(p[2 * j]);
			slc.im.at<float>(i, j) = ReverseFloat(p[2 * j + 1]);
		}
	}
	if (p) free(p);
	if (fp) fclose(fp);
	return 0;
}

int FormatConversion::read_stateVec_from_ALOS(const char* LED_file, Mat& stateVec)
{
	if (LED_file == NULL)
	{
		fprintf(stderr, "read_stateVec_from_ALOS(): input check failed!\n");
		return -1;
	}

	//检查文件格式
	int File_descriptor_len = 720;
	int Data_set_summary_len = 4096;
	int num_stateVec = 28;
	int tmp;
	FILE* fp = NULL;
	fopen_s(&fp, LED_file, "rb");
	if (!fp)
	{
		fprintf(stderr, "read_stateVec_from_ALOS(): failed to open %s!\n", LED_file);
		return -1;
	}
	fseek(fp, File_descriptor_len + Data_set_summary_len , SEEK_SET);
	fread(&tmp, 4, 1, fp);
	if (Big2Little32(tmp) != 3)
	{
		fprintf(stderr, "read_stateVec_from_ALOS():  %s : unknown format!\n", LED_file);
		if (fp) fclose(fp);
		return -1;
	}
	
	stateVec.create(num_stateVec, 7, CV_64F);

	//计算GPS时间
	
	char str[10240];
	memset(str, 0, 10240);
	string year, month, day, hour, minute, second;
	hour = "0"; minute = "0"; second = "0.0";
	fseek(fp, File_descriptor_len + Data_set_summary_len + 145 - 1, SEEK_SET);
	fread(str, 1, 4, fp);
	year = str;
	memset(str, 0, 10240);
	fread(str, 1, 4, fp);
	month = str + 2;
	memset(str, 0, 10240);
	fread(str, 1, 4, fp);
	day = str + 2;
	string time = year + "-" + month + "-" + day + "T" + hour + ":" + minute + ":" + second;
	double gps_time_start0, gps_time_start, time_interval; char* ptr;
	if (return_check(utc2gps(time.c_str(), &gps_time_start0), "utc2gps()", error_head))
	{
		if (fp) fclose(fp);
		return -1;
	}
	fseek(fp, File_descriptor_len + Data_set_summary_len + 161 - 1, SEEK_SET);
	memset(str, 0, 10240);
	fread(str, 1, 22, fp);
	gps_time_start = strtod(str, &ptr);
	gps_time_start = gps_time_start + gps_time_start0;
	memset(str, 0, 10240);
	fread(str, 1, 22, fp);
	time_interval = strtod(str, &ptr);
	fseek(fp, File_descriptor_len + Data_set_summary_len + 387 - 1, SEEK_SET);
	memset(str, 0, 10240);
	fread(str, 1, 22 * 6 * num_stateVec, fp);
	if (fp)fclose(fp);
	ptr = str;
	for (int i = 0; i < num_stateVec; i++)
	{
		stateVec.at<double>(i, 0) = gps_time_start + i * time_interval;
		stateVec.at<double>(i, 1) = strtod(ptr, &ptr);
		stateVec.at<double>(i, 2) = strtod(ptr, &ptr);
		stateVec.at<double>(i, 3) = strtod(ptr, &ptr);
		stateVec.at<double>(i, 4) = strtod(ptr, &ptr);
		stateVec.at<double>(i, 5) = strtod(ptr, &ptr);
		stateVec.at<double>(i, 6) = strtod(ptr, &ptr);
	}
	return 0;
}

int FormatConversion::read_conversion_coefficient_from_ALOS(const char* LED_file, Mat& lon_coefficient, Mat& lat_coefficient, Mat& row_coefficient, Mat& col_coefficient)
{
	if (LED_file == NULL)
	{
		fprintf(stderr, "read_conversion_coefficient_from_ALOS(): input check failed!\n");
		return -1;
	}
	//检查文件

	int ret, record_len, data_set_summary_len, platform_pos_len,
		attitue_data_len, radiometric_data_len, facility_related_record_len, offset;
	Mat tmp = Mat::zeros(1, 32, CV_64F);
	tmp.copyTo(lon_coefficient);
	tmp.copyTo(lat_coefficient);
	tmp.copyTo(row_coefficient);
	tmp.copyTo(col_coefficient);
	char* ptr; char str[20480];
	memset(str, 0, 20480);
	FILE* fp = NULL;
	fopen_s(&fp, LED_file, "rb");
	if (!fp)
	{
		fprintf(stderr, "read_conversion_coefficient_from_ALOS(): failed to open %s!\n", LED_file);
		return -1;
	}

	fseek(fp, 9 - 1, SEEK_SET);
	fread(&record_len, 4, 1, fp);
	record_len = Big2Little32(record_len);
	if (record_len != 720)
	{
		fprintf(stderr, "read_conversion_coefficient_from_ALOS(): %s: unknown format!\n", LED_file);
		if (fp) fclose(fp);
		return -1;
	}

	//计算偏移量

	fseek(fp, 187 - 1, SEEK_SET);
	fread(str, 1, 6, fp);
	data_set_summary_len = strtol(str, &ptr, 0);
	memset(str, 0, 20480);
	fseek(fp, 211 - 1, SEEK_SET);
	fread(str, 1, 6, fp);
	platform_pos_len = strtol(str, &ptr, 0);
	memset(str, 0, 20480);
	fseek(fp, 223 - 1, SEEK_SET);
	fread(str, 1, 6, fp);
	attitue_data_len = strtol(str, &ptr, 0);
	memset(str, 0, 20480);
	fseek(fp, 235 - 1, SEEK_SET);
	fread(str, 1, 6, fp);
	radiometric_data_len = strtol(str, &ptr, 0);
	facility_related_record_len = 0;
	offset = record_len + data_set_summary_len + platform_pos_len + attitue_data_len + radiometric_data_len;

	//读取数据

	while (facility_related_record_len != 5000)
	{
		offset += facility_related_record_len;
		fseek(fp, offset + 9 - 1, SEEK_SET);
		fread(&facility_related_record_len, 4, 1, fp);
		facility_related_record_len = Big2Little32(facility_related_record_len);
	}
	fseek(fp, offset + 1025 - 1, SEEK_SET);
	memset(str, 0, 20480);
	fread(str, 1, 1040 * 2, fp);
	if (fp)fclose(fp);
	ptr = str;
	for (int i = 0; i < 25; i++)
	{
		lat_coefficient.at<double>(0, 30 - i) = strtod(ptr, &ptr);
	}
	for (int i = 0; i < 25; i++)
	{
		lon_coefficient.at<double>(0, 30 - i) = strtod(ptr, &ptr);
	}
	lat_coefficient.at<double>(0, 0) = 0;
	lat_coefficient.at<double>(0, 1) = 1.0;
	lat_coefficient.at<double>(0, 2) = strtod(ptr, &ptr);
	lat_coefficient.at<double>(0, 3) = 1.0;
	lat_coefficient.at<double>(0, 4) = strtod(ptr, &ptr);
	lat_coefficient.at<double>(0, 5) = 1.0;

	lon_coefficient.at<double>(0, 0) = 0;
	lon_coefficient.at<double>(0, 1) = 1.0;
	lon_coefficient.at<double>(0, 2) = lat_coefficient.at<double>(0, 2);
	lon_coefficient.at<double>(0, 3) = 1.0;
	lon_coefficient.at<double>(0, 4) = lat_coefficient.at<double>(0, 4);
	lon_coefficient.at<double>(0, 5) = 1.0;

	for (int i = 0; i < 25; i++)
	{
		col_coefficient.at<double>(0, 30 - i) = strtod(ptr, &ptr);
	}
	for (int i = 0; i < 25; i++)
	{
		row_coefficient.at<double>(0, 30 - i) = strtod(ptr, &ptr);
	}
	col_coefficient.at<double>(0, 0) = 0;
	col_coefficient.at<double>(0, 1) = 1.0;
	col_coefficient.at<double>(0, 4) = strtod(ptr, &ptr);
	col_coefficient.at<double>(0, 3) = 1.0;
	col_coefficient.at<double>(0, 2) = strtod(ptr, &ptr);
	col_coefficient.at<double>(0, 5) = 1.0;

	row_coefficient.at<double>(0, 0) = 0;
	row_coefficient.at<double>(0, 1) = 1.0;
	row_coefficient.at<double>(0, 2) = col_coefficient.at<double>(0, 2);
	row_coefficient.at<double>(0, 3) = 1.0;
	row_coefficient.at<double>(0, 4) = col_coefficient.at<double>(0, 4);
	row_coefficient.at<double>(0, 5) = 1.0;

	return 0;
}

int FormatConversion::ALOS2h5(const char* IMG_file, const char* LED_file, const char* dst_h5)
{
	if (IMG_file == NULL ||
		LED_file == NULL ||
		dst_h5 == NULL)
	{
		fprintf(stderr, "ALOS2h5(): input check failed!\n");
		return -1;
	}

	///////////////////////创建h5文件////////////

	int ret;
	if (return_check(creat_new_h5(dst_h5), "creat_new_h5()", error_head)) return -1;

	//////////////读取slc数据并写入到目标文件中/////////////

	ComplexMat slc;
	ret = read_slc_from_ALOS(IMG_file, slc);
	if (return_check(ret, "read_slc_from_ALOS()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "s_re", slc.re);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "s_im", slc.im);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	int rows = slc.GetRows(); int cols = slc.GetCols();
	slc.re.release();
	slc.im.release();

	//////////////读取并写入轨道数据//////////////////////

	Mat stateVec;
	ret = read_stateVec_from_ALOS(LED_file, stateVec);
	if (return_check(ret, "read_stateVec_from_ALOS()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "state_vec", stateVec);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	////////////读取并写入图像坐标与经纬坐标转换关系////////////////////

	Mat lon_coef, lat_coef, row_coef, col_coef;
	ret = read_conversion_coefficient_from_ALOS(LED_file, lon_coef, lat_coef, row_coef, col_coef);
	if (return_check(ret, "read_conversion_coefficient_from_ALOS()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "lon_coefficient", lon_coef);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "lat_coefficient", lat_coef);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "row_coefficient", row_coef);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;
	ret = write_array_to_h5(dst_h5, "col_coefficient", col_coef);
	if (return_check(ret, "write_array_to_h5()", error_head)) return -1;


	//////////打开LED文件/////////////

	FILE* fp = NULL;
	fopen_s(&fp, LED_file, "rb");
	if (!fp)
	{
		fprintf(stderr, "ALOS2h5(): failed to open %s!\n", LED_file);
		return -1;
	}

	///////////读取并写入多普勒中心频率（ALOS只给了多普勒中心频率沿距离向拟合系数）/////////

	char* ptr;
	Mat tmp = Mat::zeros(1, 1, CV_64F);
	char str[2048]; memset(str, 0, 2048);
	fseek(fp, 720 + 1479 - 1, SEEK_SET);
	fread(str, 1, 32, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr);
	write_array_to_h5(dst_h5, "doppler_coefficient_a", tmp);
	tmp.at<double>(0, 0) = strtod(ptr, &ptr);
	write_array_to_h5(dst_h5, "doppler_coefficient_b", tmp);

	///////////下视角拟合系数//////////////

	Mat inc_coefficient_r = Mat::zeros(1, 11, CV_64F);
	fseek(fp, 720 + 1887 - 1, SEEK_SET); memset(str, 0, 2048);
	fread(str, 1, 120, fp);
	double factor = 180.0 / 3.1415926535;
	inc_coefficient_r.at<double>(0, 0) = 0.0;
	inc_coefficient_r.at<double>(0, 1) = 1.0;
	inc_coefficient_r.at<double>(0, 2) = 0.0;
	inc_coefficient_r.at<double>(0, 3) = 1.0;
	inc_coefficient_r.at<double>(0, 4) = strtod(str, &ptr) * factor;
	inc_coefficient_r.at<double>(0, 5) = strtod(ptr, &ptr) * factor;
	inc_coefficient_r.at<double>(0, 6) = strtod(ptr, &ptr) * factor;
	inc_coefficient_r.at<double>(0, 7) = strtod(ptr, &ptr) * factor;
	inc_coefficient_r.at<double>(0, 8) = strtod(ptr, &ptr) * factor;
	inc_coefficient_r.at<double>(0, 9) = strtod(ptr, &ptr) * factor;
	write_array_to_h5(dst_h5, "inc_coefficient_r", inc_coefficient_r);

	///////////写入其他辅助参数//////////////

	string file_type, sensor, polarization, imaging_mode,
		lookside, orbit_dir, acquisition_start_time, acquisition_stop_time, process_state;
	int record_len, data_set_summary_len, platform_pos_len,
		attitue_data_len, radiometric_data_len;
	record_len = 720; data_set_summary_len = 4096; platform_pos_len = 4680; radiometric_data_len = 9860;
	//数据类型
	ret = write_str_to_h5(dst_h5, "file_type", "SLC");
	if (return_check(ret, "write_array_to_h5()", error_head))
	{
		if (fp) fclose(fp);
		return -1;
	}
	//卫星名称
	fseek(fp, 49 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 3, fp);
	sensor = str;
	write_str_to_h5(dst_h5, "sensor", sensor.c_str());
	//拍摄起始时间
	string year, month, day, hour, minute, second;
	int h, m; double s;
	fseek(fp, record_len + data_set_summary_len + 145 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 4, fp);
	year = str;
	memset(str, 0, 2048);
	fread(str, 1, 4, fp);
	if (str[2] == ' ') str[2] = '0';
	month = str + 2;
	memset(str, 0, 2048);
	fread(str, 1, 4, fp);
	if (str[2] == ' ') str[2] = '0';
	day = str + 2;
	memset(str, 0, 2048);
	fread(str, 1, 26, fp);
	s = strtod(str + 4, &ptr);
	h = (int)floor(s / 3600.0);
	m = (int)floor((s - h * 3600) / 60.0);
	s = s - h * 3600.0 - m * 60.0;
	memset(str, 0, 2048);
	if (h < 10)sprintf(str, "0%d", h);
	else sprintf(str, "%d", h);
	hour = str;
	memset(str, 0, 2048);
	if(m < 10) sprintf(str, "0%d", m);
	else sprintf(str, "%d", m);
	minute = str;
	memset(str, 0, 2048);
	if (s < 10.0)sprintf(str, "%.5lf", s);
	else sprintf(str, "0%.5lf", s);
	second = str;
	string time = year + "-" + month + "-" + day + "T" + hour + ":" + minute + ":" + second;
	write_str_to_h5(dst_h5, "acquisition_start_time", time.c_str());
	//极化方式
	FILE* fp1 = NULL;
	fopen_s(&fp1, IMG_file, "rb");
	if (fp1)
	{
		short trans, recv;
		fseek(fp1, 720 + 53 - 1, SEEK_SET);
		fread(&trans, 2, 1, fp1);
		trans = Big2Little16(trans);
		fread(&recv, 2, 1, fp1);
		recv = Big2Little16(recv);
		string t, r;
		t = (trans == 0 ? "H" : "V");
		r = (recv == 0 ? "H" : "V");
		polarization = r + t;
		write_str_to_h5(dst_h5, "polarization", polarization.c_str());

		//左上角经纬度（第一个像素）
		int lon_topleft, lat_topleft, lon_bottomright, lat_bottomright, center_lat, center_lon, slant_range_first_pixel;
		Mat temp = Mat::zeros(1, 1, CV_64F);
		fseek(fp1, 720 + 193 - 1, SEEK_SET);
		fread(&lat_topleft, 4, 1, fp1);
		lat_topleft = Big2Little32(lat_topleft);
		temp.at<double>(0, 0) = (double)lat_topleft / 1e6;
		write_array_to_h5(dst_h5, "scene_topleft_lat", temp);

		fseek(fp1, 720 + 205 - 1, SEEK_SET);
		fread(&lon_topleft, 4, 1, fp1);
		lon_topleft = Big2Little32(lon_topleft);
		temp.at<double>(0, 0) = (double)lon_topleft / 1e6;
		write_array_to_h5(dst_h5, "scene_topleft_lon", temp);

		//右上角经纬度
		fseek(fp1, 720 + 201 - 1, SEEK_SET);
		fread(&lat_bottomright, 4, 1, fp1);
		lat_bottomright = Big2Little32(lat_bottomright);
		temp.at<double>(0, 0) = (double)lat_bottomright / 1e6;
		write_array_to_h5(dst_h5, "scene_topright_lat", temp);

		fseek(fp1, 720 + 213 - 1, SEEK_SET);
		fread(&lon_bottomright, 4, 1, fp1);
		lon_bottomright = Big2Little32(lon_bottomright);
		temp.at<double>(0, 0) = (double)lon_bottomright / 1e6;
		write_array_to_h5(dst_h5, "scene_topright_lon", temp);

		//最近斜距
		fseek(fp1, 720 + 117 - 1, SEEK_SET);
		fread(&slant_range_first_pixel, 4, 1, fp1);
		slant_range_first_pixel = Big2Little32(slant_range_first_pixel);
		write_int_to_h5(dst_h5, "slant_range_first_pixel", slant_range_first_pixel);

		if (fp1) fclose(fp1);
	}
	//拍摄模式
	fseek(fp, 720 + 413 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 32, fp);
	if (str[4] == '2')//ALOS2
	{
		attitue_data_len = 16384;
		//视向
		string filename = IMG_file;
		filename = filename.substr(filename.find_last_of('\\') + 1);
		//拍摄模式
		imaging_mode = filename.substr(29, 3);
		write_str_to_h5(dst_h5, "imaging_mode", imaging_mode.c_str());
		if (filename[32] == 'L')//左视
		{
			lookside = "LEFT";
		}
		else
		{
			lookside = "RIGHT";
		}
		
		//轨道方向
		if (filename[38] == 'D')
		{
			orbit_dir = "Descending";
		}
		else
		{
			orbit_dir = "Ascending";
		}
	}
	else//ALOS1
	{
		attitue_data_len = 8192;
		lookside = "RIGHT";
		imaging_mode = str[10];
		write_str_to_h5(dst_h5, "imaging_mode", imaging_mode.c_str());
		string filename = IMG_file;
		if (filename[filename.length() - 1] == 'A')
		{
			orbit_dir = "Ascending";
		}
		else
		{
			orbit_dir = "Descending";
		}
	}
	write_str_to_h5(dst_h5, "lookside", lookside.c_str());
	write_str_to_h5(dst_h5, "orbit_dir", orbit_dir.c_str());
	//载频
	fseek(fp, 720 + 501 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 16, fp);
	tmp.at<double>(0, 0) = 3e8 / strtod(str, &ptr);
	write_array_to_h5(dst_h5, "carrier_frequency", tmp);
	//中心下视角
	fseek(fp, 720 + 485 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 8, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr);
	write_array_to_h5(dst_h5, "incidence_center", tmp);
	//PRF
	fseek(fp, 720 + 935 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 16, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr) / 1000;
	write_array_to_h5(dst_h5, "prf", tmp);
	//方位向像素点数
	Mat tmp_int = Mat::zeros(1, 1, CV_32S);
	tmp_int.at<int>(0, 0) = rows;
	write_array_to_h5(dst_h5, "azimuth_len", tmp_int);
	//距离向像素点数
	tmp_int.at<int>(0, 0) = cols;
	write_array_to_h5(dst_h5, "range_len", tmp_int);
	//距离向分辨率
	int offset = record_len + data_set_summary_len + platform_pos_len + attitue_data_len + radiometric_data_len;
	fseek(fp, offset + 127 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 16, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr);
	write_array_to_h5(dst_h5, "range_resolution", tmp);
	//方位向分辨率
	memset(str, 0, 2048);
	fread(str, 1, 16, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr);
	write_array_to_h5(dst_h5, "azimuth_resolution", tmp);
	//方位向采样间隔
	fseek(fp, 720 + 1687 - 1, SEEK_SET);
	memset(str, 0, 2048);
	fread(str, 1, 32, fp);
	tmp.at<double>(0, 0) = strtod(str, &ptr);
	write_array_to_h5(dst_h5, "azimuth_spacing", tmp);
	//距离向采样间隔
	tmp.at<double>(0, 0) = strtod(ptr, &ptr);
	write_array_to_h5(dst_h5, "range_spacing", tmp);

	//处理等级
	ret = write_str_to_h5(dst_h5, "process_state", "InSAR_0");
	//处理描述
	int t = attitue_data_len > 10000 ? 2 : 1;
	char strings[1024]; memset(strings, 0, 1024);
	sprintf(strings, "import from ALOS%d Single Look Complex, unprocessed.", t);
	string sss(strings);
	write_str_to_h5(dst_h5, "comment", sss.c_str());

	if (fp)fclose(fp);
	return 0;
}


















XMLFile::XMLFile()
{
	memset(m_xmlFileName, 0, 2048);
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "XMLFILE_DLL_ERROR: error happens when using ");
}

XMLFile::~XMLFile()
{
}

int XMLFile::XMLFile_creat_new_project(const char* project_path, const char* project_name, const char* project_version)
{
	if (project_path == NULL ||
		project_name == NULL ||
		project_version == NULL
		)
	{
		fprintf(stderr, "XMLFile_creat_new_project(): input check failed!\n");
		return -1;
	}
	TiXmlDeclaration* declaration = new TiXmlDeclaration("1.0", "UTF-8", "yes");
	doc.LinkEndChild(declaration);
	TiXmlElement* Root = new TiXmlElement("Root");
	doc.LinkEndChild(Root);
	TiXmlElement* prj_info_node = new TiXmlElement("project_info");
	Root->LinkEndChild(prj_info_node);
	prj_info_node->SetAttribute("version", "1.0");
	TiXmlElement* prj_name_node = new TiXmlElement("project_name");
	TiXmlText* content = new TiXmlText(project_name);
	prj_name_node->LinkEndChild(content);
	prj_info_node->LinkEndChild(prj_name_node);
	content = new TiXmlText(project_path);
	TiXmlElement* prj_path_node = new TiXmlElement("project_path");
	prj_path_node->LinkEndChild(content);
	prj_info_node->LinkEndChild(prj_path_node);

	string path(project_path); string name(project_name);
	string filename = path + "\\" + name;
	std::replace(filename.begin(), filename.end(), '/', '\\');
	this->doc.SaveFile(filename.c_str());
	return 0;
}

int XMLFile::XMLFile_add_origin(
	const char* datanode_node,
	const char* node_name,
	const char* node_path,
	const char* sensor
)
{
	if (node_name == NULL ||
		node_path == NULL ||
		sensor == NULL ||
		datanode_node == NULL
		)
	{
		fprintf(stderr, "XMLFile_add_origin(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	TiXmlElement* p = doc.RootElement();
	int ret = find_node_with_attribute(p, "DataNode", "name", datanode_node, DataNode);
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		doc.RootElement()->LinkEndChild(DataNode);
		DataNode->SetAttribute("name", datanode_node);
		DataNode->SetAttribute("index", "1");
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "import");
		DataNode->SetAttribute("rank", "complex-0.0");
		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-0.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		Row_Offset->LinkEndChild(new TiXmlText("0"));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		Col_Offset->LinkEndChild(new TiXmlText("0"));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Sensor = new TiXmlElement("Sensor");
		Sensor->LinkEndChild(new TiXmlText(sensor));
		Data_Processing_Parameters->LinkEndChild(Sensor);
	}
	else
	{
		data_node_count = atoi(DataNode->Attribute("data_count")) + 1;
		string tmp;
		tmp = int2str(data_node_count);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");


		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-0.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		Row_Offset->LinkEndChild(new TiXmlText("0"));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		Col_Offset->LinkEndChild(new TiXmlText("0"));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}











int XMLFile::XMLFile_add_cut(
	const char* datanode_name,
	const char* node_name,
	const char* node_path,
	int Row_offset,
	int Col_offset, 
	double lon, 
	double lat, 
	double width,
	double height,
	const char* data_rank
)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL
		)
	{
		fprintf(stderr, "XMLFile_add_cut(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		tmp = root->Value();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			int i = strcmp(root->Attribute("rank"), "import") == 0;
			int j = strcmp(root->Attribute("rank"), "cut") == 0;
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "cut");
		DataNode->SetAttribute("rank", data_rank);
		
		
		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText(data_rank));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Center_lon = new TiXmlElement("Center_lon");
		Center_lon->SetAttribute("unit", "degree");
		char tmp2[100];
		sprintf_s(tmp2, "%.2f", lon);
		Center_lon->LinkEndChild(new TiXmlText(tmp2));
		Data_Processing_Parameters->LinkEndChild(Center_lon);
		TiXmlElement* Center_lat = new TiXmlElement("Center_lat");
		Center_lat->SetAttribute("unit", "degree");
		sprintf_s(tmp2, "%.2f", lat);
		Center_lat->LinkEndChild(new TiXmlText(tmp2));
		Data_Processing_Parameters->LinkEndChild(Center_lat);
		TiXmlElement* Width = new TiXmlElement("Width");
		Width->SetAttribute("unit", "m");
		sprintf_s(tmp2, "%.2f", width);
		Width->LinkEndChild(new TiXmlText(tmp2));
		Data_Processing_Parameters->LinkEndChild(Width);
		TiXmlElement* Height = new TiXmlElement("Height");
		Height->SetAttribute("unit", "m");
		sprintf_s(tmp2, "%.2f", height);
		Height->LinkEndChild(new TiXmlText(tmp2));
		Data_Processing_Parameters->LinkEndChild(Height);
		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText(data_rank));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_regis(const char* datanode_name, const char* node_name, const char* node_path, int Row_offset, int Col_offset, int master_index, int interp_times, int block_size, const char* temporal_baseline, const char* B_effect, const char* B_parallel)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL ||
		B_effect == NULL ||
		B_parallel == NULL
		)
	{
		fprintf(stderr, "XMLFile_add_regis(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "coregistration");
		DataNode->SetAttribute("rank", "complex-2.0");
		
		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Master_index = new TiXmlElement("master_image");
		tmp = int2str(master_index);
		Master_index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Master_index);
		TiXmlElement* Blocksize = new TiXmlElement("blocksize");
		tmp = int2str(block_size);
		Blocksize->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Blocksize);
		TiXmlElement* Interp_times = new TiXmlElement("interp_times");
		tmp = int2str(interp_times);
		Interp_times->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Interp_times);
		TiXmlElement* Temporal_baseline = new TiXmlElement("temporal_baseline_distribution");
		Temporal_baseline->SetAttribute("unit", "day");
		Temporal_baseline->LinkEndChild(new TiXmlText(temporal_baseline));
		Data_Processing_Parameters->LinkEndChild(Temporal_baseline);
		TiXmlElement* V_baseline = new TiXmlElement("effect_baseline_distribution");
		V_baseline->SetAttribute("unit", "m");
		V_baseline->LinkEndChild(new TiXmlText(B_effect));
		Data_Processing_Parameters->LinkEndChild(V_baseline);
		TiXmlElement* H_baseline = new TiXmlElement("parallel_baseline_distribution");
		H_baseline->SetAttribute("unit", "m");
		H_baseline->LinkEndChild(new TiXmlText(B_parallel));
		Data_Processing_Parameters->LinkEndChild(H_baseline);
		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_backgeocoding(const char* dataNode, const char* dataName, const char* dataPath, int masterIndex)
{
	if (!dataNode ||
		!dataName ||
		!dataPath
		)
	{
		fprintf(stderr, "XMLFile_add_backgeocoding(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", dataNode, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", dataNode);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "coregistration");
		DataNode->SetAttribute("rank", "complex-2.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Master_index = new TiXmlElement("master_image");
		tmp = int2str(masterIndex);
		Master_index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Master_index);
		
		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_SLC_deramp(const char* dataNode, const char* dataName, const char* dataPath, int masterIndex)
{
	if (!dataNode ||
		!dataName ||
		!dataPath
		)
	{
		fprintf(stderr, "XMLFile_add_SLC_deramp(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", dataNode, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", dataNode);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-3.0") == 0
				)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "SLC_deramp");
		DataNode->SetAttribute("rank", "complex-3.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-3.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Master_index = new TiXmlElement("master_image");
		tmp = int2str(masterIndex);
		Master_index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Master_index);

		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-3.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_SBAS(const char* dataNode, const char* dataName, const char* dataPath)
{
	if (!dataNode ||
		!dataName ||
		!dataPath
		)
	{
		fprintf(stderr, "XMLFile_add_SBAS(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	string tmp;
	DataNode = new TiXmlElement("DataNode");
	DataNode->SetAttribute("name", dataNode);
	int index = 1;
	TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
	for (; root != NULL; root = root->NextSiblingElement(), index++)
	{
		if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
			strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
			strcmp(root->Attribute("rank"), "complex-3.0") == 0 ||
			strcmp(root->Attribute("rank"), "phase-1.0") == 0 ||
			strcmp(root->Attribute("rank"), "phase-2.0") == 0 ||
			strcmp(root->Attribute("rank"), "phase-3.0") == 0 ||
			strcmp(root->Attribute("rank"), "dem-1.0") == 0 ||
			strcmp(root->Attribute("rank"), "SBAS-1.0") == 0
			)
			continue;
		else
			break;
	}
	string index_str = int2str(index);
	DataNode->SetAttribute("index", index_str.c_str());
	DataNode->SetAttribute("data_count", "1");
	DataNode->SetAttribute("data_processing", "SBAS");
	DataNode->SetAttribute("rank", "SBAS-1.0");

	TiXmlElement* Data = new TiXmlElement("Data");
	DataNode->LinkEndChild(Data);
	TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
	Data_Name->LinkEndChild(new TiXmlText(dataName));
	Data->LinkEndChild(Data_Name);
	TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
	Data_Rank->LinkEndChild(new TiXmlText("SBAS-1.0"));
	Data->LinkEndChild(Data_Rank);
	TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
	Data_Index->LinkEndChild(new TiXmlText("1"));
	Data->LinkEndChild(Data_Index);
	TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
	Data_Path->LinkEndChild(new TiXmlText(dataPath));
	Data->LinkEndChild(Data_Path);
	TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
	tmp = int2str(0);
	Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
	Data->LinkEndChild(Row_Offset);
	TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
	tmp = int2str(0);
	Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
	Data->LinkEndChild(Col_Offset);

	TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
	DataNode->LinkEndChild(Data_Processing_Parameters);
	TiXmlElement* nill = new TiXmlElement("nill");
	tmp = int2str(0);
	nill->LinkEndChild(new TiXmlText(tmp.c_str()));
	Data_Processing_Parameters->LinkEndChild(nill);

	if (!root)
	{
		doc.RootElement()->LinkEndChild(DataNode);
	}
	else
	{
		doc.RootElement()->InsertBeforeChild(root, *DataNode);
		while (root)
		{
			index_str = int2str(++index);
			root->SetAttribute("index", index_str.c_str());
			root = root->NextSiblingElement();
		}

	}
	return 0;
}

int XMLFile::XMLFile_add_S1_Deburst(const char* dataNode, const char* dataName, const char* dataPath)
{
	if (!dataNode ||
		!dataName ||
		!dataPath
		)
	{
		fprintf(stderr, "XMLFile_add_S1_Deburst(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", dataNode, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", dataNode);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "deburst");
		DataNode->SetAttribute("rank", "complex-1.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-1.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* nill = new TiXmlElement("nill");
		tmp = int2str(0);
		nill->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(nill);

		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(dataName));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-1.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(dataPath));
		Data->LinkEndChild(Data_Path);
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(0);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(0);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_interferometric_phase(const char* datanode_name, const char* node_name, const char* node_path,
	const char* master_name, const char* rank, int offset_row, int offset_col, int isdeflat, int istopo_removal, int iscoherence,
	int win_w, int win_h, int multilook_rg, int multilook_az)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL ||
		master_name == NULL ||
		rank == NULL)
	{
		fprintf(stderr, "XMLFile_add_interferometric_phase(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		tmp = root->Value();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-1.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "interferometric_formation");
		DataNode->SetAttribute("rank", "phase-1.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText(rank));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(offset_row);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(offset_col);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Master_image = new TiXmlElement("master_image");
		Master_image->LinkEndChild(new TiXmlText(master_name));
		Data_Processing_Parameters->LinkEndChild(Master_image);
		TiXmlElement* Isdeflat = new TiXmlElement("IsDeflat");
		tmp = int2str(isdeflat);
		Isdeflat->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Isdeflat);
		TiXmlElement* Istopo_removal = new TiXmlElement("IsTopo_Removal");
		tmp = int2str(istopo_removal);
		Istopo_removal->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Istopo_removal);
		TiXmlElement* Iscoherence = new TiXmlElement("IsCoherence");
		tmp = int2str(iscoherence);
		Iscoherence->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Iscoherence);
		TiXmlElement* Win_w = new TiXmlElement("coh_est_width");
		tmp = int2str(win_w);
		Win_w->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Win_w);
		TiXmlElement* Win_h = new TiXmlElement("coh_est_height");
		tmp = int2str(win_h);
		Win_h->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Win_h);
		TiXmlElement* Multilook_rg = new TiXmlElement("multilook_rg");
		tmp = int2str(multilook_rg);
		Multilook_rg->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data_Processing_Parameters->LinkEndChild(Multilook_rg);
		TiXmlElement* Multilook_az = new TiXmlElement("multilook_az");
		tmp = int2str(multilook_az);
		Multilook_az->LinkEndChild(new TiXmlText(tmp.c_str()));

		Data_Processing_Parameters->LinkEndChild(Multilook_az);

		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText(rank));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(offset_row);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(offset_col);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}


int XMLFile::XMLFile_add_denoise(const char* datanode_name, const char* node_name, const char* node_path, int Row_offset, int Col_offset, const char* method, int Slop_win, int Pre_win, int Goldstein_win, int Goldstein_filled_win, double alpha, const char* filter_dl_path, const char* dl_model_file, const char* tmp_path)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL ||
		method == NULL)
	{
		fprintf(stderr, "XMLFile_add_denoise(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		tmp = root->Value();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-2.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "denoise");
		DataNode->SetAttribute("rank", "phase-2.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("phase-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Method = new TiXmlElement("method");
		Method->LinkEndChild(new TiXmlText(method));
		Data_Processing_Parameters->LinkEndChild(Method);

		if (strcmp(method, "Slope") == 0)
		{
			TiXmlElement* Slop_wnd_size = new TiXmlElement("filter_wnd_size");
			tmp = int2str(Slop_win);
			Slop_wnd_size->LinkEndChild(new TiXmlText(tmp.c_str()));
			Data_Processing_Parameters->LinkEndChild(Slop_wnd_size);

			TiXmlElement* Pre_wnd_size = new TiXmlElement("filter_wnd_size");
			tmp = int2str(Pre_win);
			Pre_wnd_size->LinkEndChild(new TiXmlText(tmp.c_str()));
			Data_Processing_Parameters->LinkEndChild(Pre_wnd_size);
		}
		else if (strcmp(method, "Goldstein") == 0)
		{
			char tmp2[100];
			TiXmlElement* Goldstein_wnd_size = new TiXmlElement("FFT_size");
			tmp = int2str(Goldstein_win);
			Goldstein_wnd_size->LinkEndChild(new TiXmlText(tmp.c_str()));
			Data_Processing_Parameters->LinkEndChild(Goldstein_wnd_size);

			TiXmlElement* Goldstein_filled_size = new TiXmlElement("n_pad");
			tmp = int2str(Goldstein_filled_win);
			Goldstein_filled_size->LinkEndChild(new TiXmlText(tmp.c_str()));
			Data_Processing_Parameters->LinkEndChild(Goldstein_filled_size);

			TiXmlElement* Alpha = new TiXmlElement("alpha");
			sprintf_s(tmp2, "%.2f", alpha);
			Alpha->LinkEndChild(new TiXmlText(tmp2));
			Data_Processing_Parameters->LinkEndChild(Alpha);
		}
		else if (strcmp(method, "DL") == 0)
		{

			TiXmlElement* DL_path = new TiXmlElement("filter_dl_path");
			DL_path->LinkEndChild(new TiXmlText(filter_dl_path));
			Data_Processing_Parameters->LinkEndChild(DL_path);

			TiXmlElement* Model_path = new TiXmlElement("dl_model_file");
			Model_path->LinkEndChild(new TiXmlText(dl_model_file));
			Data_Processing_Parameters->LinkEndChild(Model_path);

			TiXmlElement* Tmp_path = new TiXmlElement("tmp_path");
			Tmp_path->LinkEndChild(new TiXmlText(tmp_path));
			Data_Processing_Parameters->LinkEndChild(Tmp_path);
		}

		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("phase-2.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_unwrap(const char* datanode_name, const char* node_name, const char* node_path, int Row_offset, int Col_offset, const char* method, double threshold)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL ||
		method == NULL)
	{
		fprintf(stderr, "XMLFile_add_unwrap(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		tmp = root->Value();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-3.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "unwrap");
		DataNode->SetAttribute("rank", "phase-3.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("phase-3.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Method = new TiXmlElement("method");
		Method->LinkEndChild(new TiXmlText(method));
		Data_Processing_Parameters->LinkEndChild(Method);

		if (strcmp(method, "Combined") == 0)
		{
			char tmp2[200];
			TiXmlElement* Coherence_threshold = new TiXmlElement("coherence_thresh");
			sprintf_s(tmp2, "%.2f", threshold);
			Coherence_threshold->LinkEndChild(new TiXmlText(tmp2));
			Data_Processing_Parameters->LinkEndChild(Coherence_threshold);
		}
		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("phase-3.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

int XMLFile::XMLFile_add_dem(const char* datanode_name, const char* node_name, const char* node_path, int Row_offset, int Col_offset, const char* method, int times)
{
	if (datanode_name == NULL ||
		node_name == NULL ||
		node_path == NULL ||
		method == NULL)
	{
		fprintf(stderr, "XMLFile_add_dem(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	int ret = find_node_with_attribute(doc.RootElement(), "DataNode", "name", datanode_name, DataNode);
	string tmp;
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		DataNode->SetAttribute("name", datanode_name);
		int index = 1;
		TiXmlElement* root = doc.RootElement()->FirstChildElement()->NextSiblingElement();
		tmp = root->Value();
		for (; root != NULL; root = root->NextSiblingElement(), index++)
		{
			if (strcmp(root->Attribute("rank"), "complex-0.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "complex-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-1.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-2.0") == 0 ||
				strcmp(root->Attribute("rank"), "phase-3.0") == 0 ||
				strcmp(root->Attribute("rank"), "dem-1.0") == 0)
				continue;
			else
				break;
		}
		string index_str = int2str(index);
		DataNode->SetAttribute("index", index_str.c_str());
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "dem");
		DataNode->SetAttribute("rank", "dem-1.0");

		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("dem-1.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText("1"));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);

		TiXmlElement* Data_Processing_Parameters = new TiXmlElement("Data_Processing_Parameters");
		DataNode->LinkEndChild(Data_Processing_Parameters);
		TiXmlElement* Method = new TiXmlElement("method");
		Method->LinkEndChild(new TiXmlText(method));
		Data_Processing_Parameters->LinkEndChild(Method);

		if (strcmp(method, "Iteration") == 0)
		{
			
			TiXmlElement* Iter_times = new TiXmlElement("iter_times");
			tmp = int2str(times);
			Iter_times->LinkEndChild(new TiXmlText(tmp.c_str()));
			Data_Processing_Parameters->LinkEndChild(Iter_times);
		}
		if (!root)
		{
			doc.RootElement()->LinkEndChild(DataNode);
		}
		else
		{
			doc.RootElement()->InsertBeforeChild(root, *DataNode);
			while (root)
			{
				index_str = int2str(++index);
				root->SetAttribute("index", index_str.c_str());
				root = root->NextSiblingElement();
			}

		}
	}
	else
	{
		string str = DataNode->Attribute("data_count");
		int index = str2int(str) + 1;
		TiXmlElement* p = NULL;
		tmp = int2str(index);
		DataNode->SetAttribute("data_count", tmp.c_str());
		TiXmlElement* LastNode = DataNode->LastChild()->ToElement();

		TiXmlElement* Data = new TiXmlElement("Data");

		int data_count = 0;
		ret = get_children_count(DataNode, &data_count);
		//DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("dem-1.0"));
		Data->LinkEndChild(Data_Rank);
		TiXmlElement* Data_Index = new TiXmlElement("Data_Index");
		Data_Index->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Data_Index);
		TiXmlElement* Data_Path = new TiXmlElement("Data_Path");
		Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Data_Path);
		/*TiXmlElement* Last_Data_Path = new TiXmlElement("Last_Data_Path");
		Last_Data_Path->LinkEndChild(new TiXmlText(node_path));
		Data->LinkEndChild(Last_Data_Path);*/
		TiXmlElement* Row_Offset = new TiXmlElement("Row_Offset");
		tmp = int2str(Row_offset);
		Row_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Row_Offset);
		TiXmlElement* Col_Offset = new TiXmlElement("Col_Offset");
		tmp = int2str(Col_offset);
		Col_Offset->LinkEndChild(new TiXmlText(tmp.c_str()));
		Data->LinkEndChild(Col_Offset);
		DataNode->InsertBeforeChild(LastNode, *Data);
	}
	return 0;
}

string XMLFile::int2str(int n)
{
	string sResult;
	stringstream ssTmp;
	ssTmp << n;
	sResult = ssTmp.str();
	return sResult;
}

int XMLFile::str2int(const string& s)
{
	int n;
	stringstream ssTmp(s);
	ssTmp >> n;
	return n;
}

int XMLFile::XMLFile_save(const char* save_path)
{
	if (save_path == NULL)
	{
		fprintf(stderr, "XMLFile_save(): input check failed!\n");
		return -1;
	}
	doc.SaveFile(save_path);
	return 0;
}

int XMLFile::XMLFile_load(const char* xmlFileName)
{
	if (xmlFileName == NULL)
	{
		fprintf(stderr, "XMLFile_load(): input check failed!\n");
		return -1;
	}
	strcpy_s(m_xmlFileName, 2047, xmlFileName);
	if (!doc.LoadFile(xmlFileName))
	{
		fprintf(stderr, "XMLFile_load(): can't load XML file %s!\n", xmlFileName);
		return -1;
	}
	return 0;
}

int XMLFile::get_root(TiXmlElement*& root)
{
	root = doc.RootElement();
	return 0;
}

int XMLFile::get_children_count(TiXmlElement* pRoot, int* count)
{
	if (!pRoot)
	{
		fprintf(stderr, "get_children_count(): Node doesn't exist\n");
		return -1;
	}
	TiXmlElement* p;
	int tmp = 0;
	for (p = pRoot->FirstChildElement(); p != NULL; p = p->NextSiblingElement())
	{
		tmp++;
	}
	*count = tmp;
	return 0;
}

int XMLFile::_find_node(TiXmlElement* pRoot, const char* node_name, TiXmlElement*& pnode)
{
	if (node_name == NULL || pRoot == NULL)
	{
		fprintf(stderr, "find_node(): input check failed!\n");
		return -1;
	}
	const char* value = pRoot->Value();
	if (strcmp(pRoot->Value(), node_name) == 0)
	{
		pnode = pRoot;
		return 0;
	}

	TiXmlElement* p = pRoot;
	for (p = p->FirstChildElement(); p != NULL; p = p->NextSiblingElement())
	{
		if (0 == _find_node(p, node_name, pnode)) return 0;
	}
	
	return -1;
}

int XMLFile::find_node(const char* node_name, TiXmlElement*& pnode)
{
	if (node_name == NULL)
	{
		fprintf(stderr, "find_node(): input check failed!\n");
		return -1;
	}
	TiXmlElement* root = doc.RootElement();
	int ret = _find_node(root, node_name, pnode);
	if (ret < 0)
	{
		fprintf(stderr, "find_node(): node %s not found!\n", node_name);
		return -1;
	}
	return 0;
}

int XMLFile::find_node_with_attribute(TiXmlElement* pRoot, const char* node_name, const char* attribute_name, const char* attribute_value, TiXmlElement*& pnode)
{
	if (node_name == NULL ||
		attribute_name == NULL ||
		attribute_value == NULL)
	{
		fprintf(stderr, "find_node_with_attribute(): input check failed!\n");
		return -1;
	}
	const char* value = pRoot->Value();
	const char* attribute = pRoot->Attribute(attribute_name);
	if (strcmp(value, node_name) == 0 && strcmp(attribute, attribute_value) == 0)
	{
		pnode = pRoot;
		return 0;
	}

	TiXmlElement* p = pRoot;
	for (p = p->FirstChildElement(); p != NULL; p = p->NextSiblingElement())
	{
		if (0 == find_node_with_attribute(p, node_name, attribute_name, attribute_value, pnode)) return 0;
	}

	return -1;
}

int XMLFile::get_str_para(const char* node_name, string& value)
{
	if (node_name == NULL)
	{
		fprintf(stderr, "get_str_para(): input check failed!\n");
		return -1;
	}
	int ret;
	TiXmlElement* pnode = NULL;
	ret = find_node(node_name, pnode);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	if (pnode)
	{
		string x(pnode->GetText());
		value = x;
	}
	return 0;
}

int XMLFile::get_double_para(const char* node_name, double* value)
{
	if (node_name == NULL || value == NULL)
	{
		fprintf(stderr, "get_double_para(): input check failed!\n");
		return -1;
	}
	int ret;
	string tmp;
	ret = get_str_para(node_name, tmp);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	ret = sscanf(tmp.c_str(), "%lf", value);
	if (ret != 1)
	{
		fprintf(stderr, "get_double_para(): %s is unknown paramter format!\n", tmp.c_str());
		return -1;
	}
	return 0;
}

int XMLFile::getDoubleArray(const char* node_name, Mat& Array, TiXmlElement* rootNode)
{
	if (!node_name)
	{
		fprintf(stderr, "getDoubleArray(): input check failed!\n");
		return -1;
	}
	int ret, count;
	char* ptr = NULL;
	count = -1;
	TiXmlElement* pnode = NULL;
	if (rootNode)
	{
		ret = _find_node(rootNode, node_name, pnode);
	}
	else
	{
		ret = find_node(node_name, pnode);
	}
	if (return_check(ret, "getDoubleArray()", error_head)) return -1;
	if (pnode)
	{
		if (!pnode->FirstAttribute())
		{
			fprintf(stderr, "getDoubleArray(): no count attribute!\n");
			return -1;
		}
		ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &count);
		if (ret != 1 || count <= 0)
		{
			fprintf(stderr, "getDoubleArray(): no count attribute!\n");
			return -1;
		}
		Array.create(1, count, CV_64F);
		Array.at<double>(0, 0) = strtod(pnode->GetText(), &ptr);
		for (int i = 1; i < count; i++)
		{
			Array.at<double>(0, i) = strtod(ptr, &ptr);
		}

	}
	return 0;
}

int XMLFile::get_int_para(const char* node_name, int* value)
{
	if (node_name == NULL || value == NULL)
	{
		fprintf(stderr, "get_int_para(): input check failed!\n");
		return -1;
	}
	int ret;
	string tmp;
	ret = get_str_para(node_name, tmp);
	if (return_check(ret, "get_int_para()", error_head)) return -1;
	ret = sscanf(tmp.c_str(), "%d", value);
	if (ret != 1)
	{
		fprintf(stderr, "get_int_para(): %s is unknown paramter format!\n", tmp.c_str());
		return -1;
	}
	return 0;
}

int XMLFile::getIntArray(const char* node_name, Mat& Array, TiXmlElement* rootNode)
{
	if (!node_name)
	{
		fprintf(stderr, "getIntArray(): input check failed!\n");
		return -1;
	}
	int ret, count;
	char* ptr = NULL;
	count = -1;
	TiXmlElement* pnode = NULL;
	if (rootNode)
	{
		ret = _find_node(rootNode, node_name, pnode);
	}
	else
	{
		ret = find_node(node_name, pnode);
	}
	if (return_check(ret, "getIntArray()", error_head)) return -1;
	if (pnode)
	{
		if (!pnode->FirstAttribute())
		{
			fprintf(stderr, "getIntArray(): no count attribute!\n");
			return -1;
		}
		ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &count);
		if (ret != 1 || count <= 0)
		{
			fprintf(stderr, "getIntArray(): no count attribute!\n");
			return -1;
		}
		Array.create(1, count, CV_32S);
		Array.at<int>(0, 0) = strtol(pnode->GetText(), &ptr, 0);
		for (int i = 1; i < count; i++)
		{
			Array.at<int>(0, i) = strtol(ptr, &ptr, 0);
		}

	}
	return 0;
}

int XMLFile::get_gcps_from_TSX(Mat& gcps)
{
	/*
	* 验证根节点
	*/
	TiXmlElement* pRoot = NULL;
	pRoot = doc.RootElement();
	if (pRoot)
	{
		if (0 != strcmp(pRoot->Value(), "geoReference"))
		{
			fprintf(stderr, "get_gcps_from_TSX():  %s: unknown format!\n", this->m_xmlFileName);
			return -1;
		}
	}
	else
	{
		fprintf(stderr, "get_gcps_from_TSX():  %s: root element error!\n", this->m_xmlFileName);
		return -1;
	}
	/*
	* 确定控制点个数
	*/
	int n_gcps = 1;
	TiXmlElement* pnode = NULL;
	int ret = find_node("numberOfGridPoints", pnode);
	if (return_check(ret, "get_gcps_from_TSX()", error_head)) return -1;
	pnode = pnode->FirstChildElement();
	string tmp(pnode->GetText());
	ret = sscanf(tmp.c_str(), "%d", &n_gcps);
	if (ret != 1)
	{
		fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
		return -1;
	}
	Mat x(n_gcps, 6, CV_64F);
	/*
	* 读取控制点
	*/
	ret = find_node("gridPoint", pnode);
	if (return_check(ret, "get_gcps_from_TSX()", error_head)) return -1;
	TiXmlElement* pchild = NULL;
	double lon, lat, height, row, col, inc;
	for (int i = 1; i <= n_gcps; i++)
	{
		if (!pnode) break;
		//读取经度
		ret = _find_node(pnode, "lon", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &lon);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//读取纬度
		ret = _find_node(pnode, "lat", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &lat);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//读取高度
		ret = _find_node(pnode, "height", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &height);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//读取行数
		ret = _find_node(pnode, "row", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &row);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//读取列数
		ret = _find_node(pnode, "col", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &col);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//读取下视角
		ret = _find_node(pnode, "inc", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		tmp = pchild->GetText();
		ret = sscanf(tmp.c_str(), "%lf", &inc);
		if (ret != 1)
		{
			fprintf(stderr, "get_gcps_from_TSX(): %s: unkown format!\n", this->m_xmlFileName);
			return -1;
		}
		//赋值
		x.at<double>(i - 1, 0) = lon;
		x.at<double>(i - 1, 1) = lat;
		x.at<double>(i - 1, 2) = height;
		x.at<double>(i - 1, 3) = row;
		x.at<double>(i - 1, 4) = col;
		x.at<double>(i - 1, 5) = inc;
		//下一个节点
		pnode = pnode->NextSiblingElement();
	}
	x.copyTo(gcps);
	return 0;
}

int XMLFile::get_stateVec_from_TSX(Mat& stateVec)
{
	TiXmlElement* pRoot = NULL;
	pRoot = doc.RootElement();
	if (!pRoot)
	{
		fprintf(stderr, "get_stateVec_from_TSX(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}

	/*
	* 确定stateVec个数
	*/

	TiXmlElement* pnode = NULL;
	int ret = _find_node(pRoot, "numStateVectors", pnode);
	if (ret < 0)
	{
		fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	int numstateVec;
	ret = sscanf(pnode->GetText(), "%d", &numstateVec);
	if (ret != 1)
	{
		fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	Mat x(numstateVec, 7, CV_64F);

	/*
	* 找到第一个stateVec
	*/
	ret = _find_node(pRoot, "stateVec", pnode);
	if (ret < 0)
	{
		fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	TiXmlElement* pchild = NULL;
	double GPS_time, posX, posY, posZ, velX, velY, velZ;
	for (int i = 0; i < numstateVec; i++)
	{
		if (!pnode) break;
		//GPS时间
		ret = _find_node(pnode, "timeUTC", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		UTC2GPS(pchild->GetText(), &GPS_time);
		//ret = sscanf(pchild->GetText(), "%lf", &GPS_time);
		//if (ret != 1)
		//{
		//	fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
		//	return -1;
		//}
		//posX
		ret = _find_node(pnode, "posX", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &posX);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//posY
		ret = _find_node(pnode, "posY", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &posY);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//posZ
		ret = _find_node(pnode, "posZ", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &posZ);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//velX
		ret = _find_node(pnode, "velX", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &velX);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//velY
		ret = _find_node(pnode, "velY", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &velY);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//velZ
		ret = _find_node(pnode, "velZ", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &velZ);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//赋值
		x.at<double>(i, 0) = GPS_time;
		x.at<double>(i, 1) = posX;
		x.at<double>(i, 2) = posY;
		x.at<double>(i, 3) = posZ;
		x.at<double>(i, 4) = velX;
		x.at<double>(i, 5) = velY;
		x.at<double>(i, 6) = velZ;
		//下一个节点
		pnode = pnode->NextSiblingElement();
	}
	x.copyTo(stateVec);
	return 0;
}

int XMLFile::get_dopplerCentroid_from_TSX(Mat& doppler)
{
	TiXmlElement* pnode = NULL;
	TiXmlElement* pchild = NULL;
	TiXmlElement* pchild2 = NULL;
	int ret, numberOfDopplerRecords, polynomialDegree;
	ret = find_node("numberOfDopplerRecords", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	ret = sscanf(pnode->GetText(), "%d", &numberOfDopplerRecords);
	if (ret != 1)
	{
		fprintf(stderr, "get_dopplerCentroid_from_TSX(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	ret = find_node("dopplerEstimate", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	ret = _find_node(pnode, "combinedDoppler", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "get_dopplerCentroid_from_TSX(): node combinedDoppler not found!\n");
		return -1;
	}
	ret = _find_node(pchild, "polynomialDegree", pnode);
	if (ret < 0)
	{
		fprintf(stderr, "get_dopplerCentroid_from_TSX(): node polynomialDegree not found!\n");
		return -1;
	}
	ret = sscanf(pnode->GetText(), "%d", &polynomialDegree);
	if (ret != 1)
	{
		fprintf(stderr, "get_dopplerCentroid_from_TSX(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}

	doppler.create(numberOfDopplerRecords, polynomialDegree + 2, CV_64F);
	ret = find_node("dopplerEstimate", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	double ref_time, c;
	for (int i = 0; i < numberOfDopplerRecords; i++)
	{
		if (!pnode) break;
		ret = _find_node(pnode, "combinedDoppler", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_dopplerCentroid_from_TSX(): node combinedDoppler not found!\n");
			return -1;
		}
		ret = _find_node(pchild, "referencePoint", pchild2);
		if (ret < 0)
		{
			fprintf(stderr, "get_dopplerCentroid_from_TSX(): node referencePoint not found!\n");
			return -1;
		}
		ret = sscanf(pchild2->GetText(), "%lf", &ref_time);
		if (ret != 1)
		{
			fprintf(stderr, "get_dopplerCentroid_from_TSX(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		doppler.at<double>(i, 0) = ref_time;
		pchild = pchild2->NextSiblingElement()->NextSiblingElement();
		for (int j = 0; j < polynomialDegree + 1; j++)
		{
			if (!pchild) break;
			ret = sscanf(pchild->GetText(), "%lf", &c);
			if (ret != 1)
			{
				fprintf(stderr, "get_dopplerCentroid_from_TSX(): %s: unknown data format!\n", this->m_xmlFileName);
				return -1;
			}
			doppler.at<double>(i, j + 1) = c;
			pchild = pchild->NextSiblingElement();
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int XMLFile::get_gcps_from_sentinel(Mat& gcps)
{
	/*
	* 确定控制点个数
	*/
	int n_gcps = 1;
	TiXmlElement* pRoot = doc.RootElement();
	int ret;
	TiXmlElement* pnode = NULL;
	ret = _find_node(pRoot, "geolocationGridPointList", pnode);
	if (ret < 0)
	{
		fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node geolocationGridPointList not found!\n");
		return -1;
	}
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &n_gcps);
	if (ret != 1)
	{
		fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	gcps.create(n_gcps, 6, CV_64F);
	TiXmlElement* pchild = NULL;
	ret = _find_node(pnode, "geolocationGridPoint", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node geolocationGridPoint not found!\n");
		return -1;
	}
	double lon, lat, height, row, col, inc;
	for (int i = 0; i < n_gcps; i++)
	{
		if (!pchild) break;

		//longtitude
		ret = _find_node(pchild, "longitude", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node longitude not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &lon);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//latitude
		ret = _find_node(pchild, "latitude", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node latitude not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &lat);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//height
		ret = _find_node(pchild, "height", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node height not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &height);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//row
		ret = _find_node(pchild, "line", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node line not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &row);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//col
		ret = _find_node(pchild, "pixel", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node pixel not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &col);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//incidence angle
		ret = _find_node(pchild, "incidenceAngle", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): node pixel not found!\n");
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lf", &inc);
		if (ret != 1)
		{
			fprintf(stderr, "XMLFile::get_gcps_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//assignment
		gcps.at<double>(i, 0) = lon;
		gcps.at<double>(i, 1) = lat;
		gcps.at<double>(i, 2) = height;
		gcps.at<double>(i, 3) = row;
		gcps.at<double>(i, 4) = col;
		gcps.at<double>(i, 5) = inc;
		pchild = pchild->NextSiblingElement();
	}
	return 0;
}

int XMLFile::get_dopplerCentroid_from_sentinel(Mat& doppler)
{
	TiXmlElement* pnode, * pchild1, * pchild2;
	int ret, numOfDcEstimates, polynomialDegree;
	ret = find_node("dopplerCentroid", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	ret = sscanf(pnode->FirstChildElement()->FirstAttribute()->Value(), "%d", &numOfDcEstimates);
	if (ret != 1)
	{
		fprintf(stderr, "get_dopplerCentroid_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	ret = _find_node(pnode, "dataDcPolynomial", pchild1);
	if (ret < 0)
	{
		fprintf(stderr, "get_dopplerCentroid_from_sentinel(): node dataDcPolynomial not found!\n");
		return -1;
	}
	ret = sscanf(pchild1->FirstAttribute()->Value(), "%d", &polynomialDegree);
	if (ret != 1)
	{
		fprintf(stderr, "get_dopplerCentroid_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}
	doppler.create(numOfDcEstimates, polynomialDegree + 1, CV_64F);
	double t0, c;
	char* ptr;
	ret = find_node("dcEstimate", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	for (int i = 0; i < numOfDcEstimates; i++)
	{
		if (!pnode) break;
		ret = _find_node(pnode, "t0", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_dopplerCentroid_from_sentinel(): node t0 not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &t0);
		if (ret != 1)
		{
			fprintf(stderr, "get_dopplerCentroid_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		doppler.at<double>(i, 0) = t0;
		ret = _find_node(pnode, "dataDcPolynomial", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_dopplerCentroid_from_sentinel(): node dataDcPolynomial not found!\n");
			return -1;
		}
		c = strtod(pchild1->GetText(), &ptr);
		doppler.at<double>(i, 1) = c;
		for (int j = 1; j < polynomialDegree; j++)
		{
			c = strtod(ptr, &ptr);
			doppler.at<double>(i, j + 1) = c;
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int XMLFile::get_stateVec_from_sentinel(Mat& stateVec)
{
	TiXmlElement* pnode, * pchild, * pchild1;
	FormatConversion conversion;
	int ret, numOfstateVec;
	ret = find_node("orbitList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &numOfstateVec);
	if (ret != 1)
	{
		fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
		return -1;
	}

	ret = find_node("orbit", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	double time, x, y, z, vx, vy, vz;
	stateVec.create(numOfstateVec, 7, CV_64F);
	for (int i = 0; i < numOfstateVec; i++)
	{
		if (!pnode) break;
		//GPS时间
		ret = _find_node(pnode, "time", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node time not found!\n");
			return -1;
		}
		ret = conversion.utc2gps(pchild->GetText(), &time);
		if (return_check(ret, "utc2gps()", error_head)) return -1;
		//位置x
		ret = _find_node(pnode, "position", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node position not found!\n");
			return -1;
		}
		ret = _find_node(pchild, "x", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node x not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &x);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//位置y
		ret = _find_node(pchild, "y", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node y not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &y);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//位置z
		ret = _find_node(pchild, "z", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node z not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &z);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//速度x
		ret = _find_node(pnode, "velocity", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node velocity not found!\n");
			return -1;
		}
		ret = _find_node(pchild, "x", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node x not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &vx);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//速度y
		ret = _find_node(pchild, "y", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node y not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &vy);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		//速度z
		ret = _find_node(pchild, "z", pchild1);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): node z not found!\n");
			return -1;
		}
		ret = sscanf(pchild1->GetText(), "%lf", &vz);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_sentinel(): %s: unknown data format!\n", this->m_xmlFileName);
			return -1;
		}

		//赋值
		stateVec.at<double>(i, 0) = time;
		stateVec.at<double>(i, 1) = x;
		stateVec.at<double>(i, 2) = y;
		stateVec.at<double>(i, 3) = z;
		stateVec.at<double>(i, 4) = vx;
		stateVec.at<double>(i, 5) = vy;
		stateVec.at<double>(i, 6) = vz;
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int FormatConversion::Copy_para_from_h5_2_h5(const char* Input_file, const char* Output_file)
{
	if (Input_file == NULL ||
		Output_file == NULL)
	{
		fprintf(stderr, "Copy_para_from_h5_2_h5(): input check failed!\n");
		return -1;
	}

	string tmp_str;;
	Mat tmp_mat;
	/*图像类型*/
	if (!read_str_from_h5(Input_file, "file_type", tmp_str))
		write_str_to_h5(Output_file, "file_type", tmp_str.c_str());
	/*卫星名称*/
	if (!read_str_from_h5(Input_file, "sensor", tmp_str))
		write_str_to_h5(Output_file, "sensor", tmp_str.c_str());
	/*极化方式*/
	if (!read_str_from_h5(Input_file, "polarization", tmp_str))
		write_str_to_h5(Output_file, "polarization", tmp_str.c_str());
	/*拍摄模式*/
	if (!read_str_from_h5(Input_file, "imaging_mode", tmp_str))
		write_str_to_h5(Output_file, "imaging_mode", tmp_str.c_str());
	/*左右视*/
	if (!read_str_from_h5(Input_file, "lookside", tmp_str))
		write_str_to_h5(Output_file, "lookside", tmp_str.c_str());
	/*轨道方向*/
	if (!read_str_from_h5(Input_file, "orbit_dir", tmp_str))
		write_str_to_h5(Output_file, "orbit_dir", tmp_str.c_str());
	/*子带号（swath）*/
	if (!read_str_from_h5(Input_file, "swath", tmp_str))
		write_str_to_h5(Output_file, "swath", tmp_str.c_str());
	/*轨道高度*/
	if (!read_array_from_h5(Input_file, "orbit_altitude", tmp_mat))
		write_array_to_h5(Output_file, "orbit_altitude", tmp_mat);
	/*载频*/
	if (!read_array_from_h5(Input_file, "carrier_frequency", tmp_mat))
		write_array_to_h5(Output_file, "carrier_frequency", tmp_mat);
	/*航偏角*/
	if (!read_array_from_h5(Input_file, "heading", tmp_mat))
		write_array_to_h5(Output_file, "heading", tmp_mat);
	/*脉冲重复频率*/
	if (!read_array_from_h5(Input_file, "prf", tmp_mat))
		write_array_to_h5(Output_file, "prf", tmp_mat);
	/*gcps数据*/
	if (!read_array_from_h5(Input_file, "gcps", tmp_mat))
		write_array_to_h5(Output_file, "gcps", tmp_mat);
	/*方位向分辨率*/
	if (!read_array_from_h5(Input_file, "azimuth_resolution", tmp_mat))
		write_array_to_h5(Output_file, "azimuth_resolution", tmp_mat);
	/*距离向分辨率*/
	if (!read_array_from_h5(Input_file, "range_resolution", tmp_mat))
		write_array_to_h5(Output_file, "range_resolution", tmp_mat);
	/*方位向采样间隔*/
	if (!read_array_from_h5(Input_file, "azimuth_spacing", tmp_mat))
		write_array_to_h5(Output_file, "azimuth_spacing", tmp_mat);
	/*距离向采样间隔*/
	if (!read_array_from_h5(Input_file, "range_spacing", tmp_mat))
		write_array_to_h5(Output_file, "range_spacing", tmp_mat);
	/*卫星轨道*/
	if (!read_array_from_h5(Input_file, "state_vec", tmp_mat))
		write_array_to_h5(Output_file, "state_vec", tmp_mat);
	/*精密轨道数据*/
	if (!read_array_from_h5(Input_file, "fine_state_vec", tmp_mat))
		write_array_to_h5(Output_file, "fine_state_vec", tmp_mat);
	/*拍摄开始时间*/
	if (!read_str_from_h5(Input_file, "acquisition_start_time", tmp_str))
		write_str_to_h5(Output_file, "acquisition_start_time", tmp_str.c_str());
	/*拍摄开始时间*/
	if (!read_str_from_h5(Input_file, "acquisition_stop_time", tmp_str))
		write_str_to_h5(Output_file, "acquisition_stop_time", tmp_str.c_str());
	/*多普勒中心频率*/
	if (!read_array_from_h5(Input_file, "doppler_centroid", tmp_mat))
		write_array_to_h5(Output_file, "doppler_centroid", tmp_mat);
	/*多普勒中心频率系数a*/
	if (!read_array_from_h5(Input_file, "doppler_coefficient_a", tmp_mat))
		write_array_to_h5(Output_file, "doppler_coefficient_a", tmp_mat);
	/*多普勒中心频率系数b*/
	if (!read_array_from_h5(Input_file, "doppler_coefficient_b", tmp_mat))
		write_array_to_h5(Output_file, "doppler_coefficient_b", tmp_mat);
	/*经度拟合系数*/
	if (!read_array_from_h5(Input_file, "lon_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "lon_coefficient", tmp_mat);
	/*纬度拟合系数*/
	if (!read_array_from_h5(Input_file, "lat_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "lat_coefficient", tmp_mat);
	/*行坐标拟合系数*/
	if (!read_array_from_h5(Input_file, "row_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "row_coefficient", tmp_mat);
	/*列坐标拟合系数b*/
	if (!read_array_from_h5(Input_file, "col_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "col_coefficient", tmp_mat);
	/*下视角拟合系数*/
	if (!read_array_from_h5(Input_file, "inc_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "inc_coefficient", tmp_mat);
	/*下视角拟合系数r*/
	if (!read_array_from_h5(Input_file, "inc_coefficient_r", tmp_mat))
		write_array_to_h5(Output_file, "inc_coefficient_r", tmp_mat);
	/*行坐标拟合系数*/
	if (!read_array_from_h5(Input_file, "row_coefficient", tmp_mat))
		write_array_to_h5(Output_file, "row_coefficient", tmp_mat);
	/*最近斜距*/
	if (!read_array_from_h5(Input_file, "slant_range_first_pixel", tmp_mat))
		write_array_to_h5(Output_file, "slant_range_first_pixel", tmp_mat);
	///*行偏移量*/
	//if (!read_array_from_h5(Input_file, "offset_row", tmp_mat))
	//	write_array_to_h5(Output_file, "offset_row", tmp_mat);
	///*列偏移量*/
	//if (!read_array_from_h5(Input_file, "offset_col", tmp_mat))
	//	write_array_to_h5(Output_file, "offset_col", tmp_mat);
	return 0;
}

/*------------------------------------------------*/
/*               哨兵一号数据读取工具             */
/*------------------------------------------------*/
Sentinel1Reader::Sentinel1Reader()
{
	bXmlLoad = false;
	isDataAvailable = false;
	memset(m_xmlFileName, 0, 2048);
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	this->azimuthPixelSpacing = 0.0;
	this->azimuthSteeringRate = 0.0;
	this->azimuthTimeInterval = 0.0;
	this->numberOfLines = 0;
	this->headingAngle = 0.0;
	this->numberOfSamples = 0;
	this->burstCount = 0;
	this->linesPerBurst = 0;
	this->pass = "";
	this->polarization = "";
	this->radarFrequency = 0.0;
	this->rangePixelSpacing = 0.0;
	this->rangeSamplingRate = 0.0;
	this->sensor = "sentinel";
	this->swath = "";
	this->slantRangeTime = 0.0;

}

Sentinel1Reader::Sentinel1Reader(const char* xmlfile, const char* tiffFile, const char* PODFile)
{
	bXmlLoad = false;
	isDataAvailable = false;
	memset(m_xmlFileName, 0, 2048);
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	this->azimuthPixelSpacing = 0.0;
	this->azimuthSteeringRate = 0.0;
	this->azimuthTimeInterval = 0.0;
	this->numberOfLines = 0;
	this->headingAngle = 0.0;
	this->numberOfSamples = 0;
	this->burstCount = 0;
	this->linesPerBurst = 0;
	this->pass = "";
	this->polarization = "";
	this->radarFrequency = 0.0;
	this->rangePixelSpacing = 0.0;
	this->rangeSamplingRate = 0.0;
	this->sensor = "sentinel";
	this->swath = "";
	this->slantRangeTime = 0.0;
	this->tiffFile = tiffFile;
	if (PODFile) this->PODFile = PODFile;
	if (!xmldoc.XMLFile_load(xmlfile)) {
		bXmlLoad = true;
		std::strcpy(this->m_xmlFileName, xmlfile);
	}
}

Sentinel1Reader::~Sentinel1Reader()
{

}

int Sentinel1Reader::load(const char* xmlfile, const char* tiffFile)
{
	if (bXmlLoad) return 0;//只允许加载1次
	this->tiffFile = tiffFile;
	if (!xmldoc.XMLFile_load(xmlfile)) {
		bXmlLoad = true;
		std::strcpy(this->m_xmlFileName, xmlfile);
	}
	else
	{
		fprintf(stderr, "Sentinel1Reader::load(): can't load file %s\n", xmlfile);
		return -1;
	}
	return 0;
}

int Sentinel1Reader::getDcEstimateList()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getDcEstimateList(): input check failed!\n");
		return -1;
	}
	TiXmlElement* pnode, * pchild1, * pchild2;
	int ret, numOfDcEstimates, polynomialDegree;
	ret = xmldoc.find_node("dopplerCentroid", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	ret = sscanf(pnode->FirstChildElement()->FirstAttribute()->Value(), "%d", &numOfDcEstimates);
	ret = xmldoc._find_node(pnode, "dataDcPolynomial", pchild1);
	ret = sscanf(pchild1->FirstAttribute()->Value(), "%d", &polynomialDegree);
	DcEstimateList.create(numOfDcEstimates, polynomialDegree + 2, CV_64F);
	double t0, c;
	char* ptr;
	ret = xmldoc.find_node("dcEstimate", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	for (int i = 0; i < numOfDcEstimates; i++)
	{
		if (!pnode) break;

		ret = xmldoc._find_node(pnode, "azimuthTime", pchild1);
		UTC2GPS(pchild1->GetText(), &t0);
		DcEstimateList.at<double>(i, 0) = t0;

		ret = xmldoc._find_node(pnode, "t0", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &t0);
		DcEstimateList.at<double>(i, 1) = t0;

		ret = xmldoc._find_node(pnode, "dataDcPolynomial", pchild1);
		c = strtod(pchild1->GetText(), &ptr);
		DcEstimateList.at<double>(i, 2) = c;
		for (int j = 1; j < polynomialDegree; j++)
		{
			c = strtod(ptr, &ptr);
			DcEstimateList.at<double>(i, j + 2) = c;
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getAzimuthFmRateList()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getAzimuthFmRateList(): input check failed!\n");
		return -1;
	}
	TiXmlElement* pnode, * pchild1, * pchild2;
	int ret, numOfFmEstimates, polynomialDegree;
	ret = xmldoc.find_node("azimuthFmRateList", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &numOfFmEstimates);
	ret = xmldoc._find_node(pnode, "azimuthFmRatePolynomial", pchild1);
	ret = sscanf(pchild1->FirstAttribute()->Value(), "%d", &polynomialDegree);
	AzimuthFmRateList.create(numOfFmEstimates, polynomialDegree + 2, CV_64F);
	double t0, c;
	char* ptr;
	ret = xmldoc.find_node("azimuthFmRate", pnode);
	if (return_check(ret, "find_node", error_head)) return -1;
	for (int i = 0; i < numOfFmEstimates; i++)
	{
		if (!pnode) break;

		ret = xmldoc._find_node(pnode, "azimuthTime", pchild1);
		UTC2GPS(pchild1->GetText(), &t0);
		AzimuthFmRateList.at<double>(i, 0) = t0;

		ret = xmldoc._find_node(pnode, "t0", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &t0);
		AzimuthFmRateList.at<double>(i, 1) = t0;

		ret = xmldoc._find_node(pnode, "azimuthFmRatePolynomial", pchild1);
		c = strtod(pchild1->GetText(), &ptr);
		AzimuthFmRateList.at<double>(i, 2) = c;
		for (int j = 1; j < polynomialDegree; j++)
		{
			c = strtod(ptr, &ptr);
			AzimuthFmRateList.at<double>(i, j + 2) = c;
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getAntennaPattern()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getAntennaPattern(): input check failed!\n");
		return -1;
	}
	TiXmlElement* pnode, * pchild1, * pchild2;
	int ret = xmldoc.find_node("antennaPatternList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	int count = -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &count);
	ret = xmldoc._find_node(pnode, "slantRangeTime", pchild1);
	int count2 = -1;
	ret = sscanf(pchild1->FirstAttribute()->Value(), "%d", &count2);
	antennaPattern_slantRangeTime.create(count, count2, CV_64F);
	antennaPattern_elevationAngle.create(count, count2, CV_64F);
	Mat Array;
	ret = xmldoc._find_node(pnode, "antennaPattern", pchild1);
	for (int i = 0; i < count; i++)
	{
		if (!pchild1) break;
		ret = xmldoc.getDoubleArray("slantRangeTime", Array, pchild1);
		Array.copyTo(antennaPattern_slantRangeTime(cv::Range(i, i + 1), cv::Range(0, count2)));
		ret = xmldoc.getDoubleArray("elevationAngle", Array, pchild1);
		Array.copyTo(antennaPattern_elevationAngle(cv::Range(i, i + 1), cv::Range(0, count2)));
		pchild1 = pchild1->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getBurstCount(int* burstCount)
{
	if (!bXmlLoad || !burstCount)
	{
		fprintf(stderr, "Sentinel1Reader::getBurstCount(): input check failed!\n");
		return -1;
	}
	int ret;
	TiXmlElement* pnode;
	ret = xmldoc.find_node("burstList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", burstCount);
	return 0;
}

int Sentinel1Reader::getBurstAzimuthTime()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getBurstAzimuthTime(): input check failed!\n");
		return -1;
	}
	int ret, burstCount; double t;
	TiXmlElement* pnode, * pchild;
	ret = getBurstCount(&burstCount);
	burstAzimuthTime.create(burstCount, 1, CV_64F);
	ret = xmldoc.find_node("burst", pnode);
	for (int i = 0; i < burstCount; i++)
	{
		if (!pnode) break;
		ret = xmldoc._find_node(pnode, "azimuthTime", pchild);
		UTC2GPS(pchild->GetText(), &t);
		burstAzimuthTime.at<double>(i, 0) = t;
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getFirstValidSample()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getFirstValidSample(): input check failed!\n");
		return -1;
	}
	int ret, burstCount;
	TiXmlElement* pnode;
	Mat tmp;
	ret = getBurstCount(&burstCount);
	if (return_check(ret, "getBurstCount()", error_head)) return -1;
	firstValidSample.create(burstCount, 1, CV_32S);
	ret = xmldoc.find_node("burst", pnode);
	for (int i = 0; i < burstCount; i++)
	{
		if (!pnode) break;
		xmldoc.getIntArray("firstValidSample", tmp, pnode);
		for (int j = 0; j < tmp.cols; j++)
		{
			if (tmp.at<int>(0, j) != -1)
			{
				firstValidSample.at<int>(i, 0) = tmp.at<int>(0, j);
				break;
			}
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getLastValidSample()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getLastValidSample(): input check failed!\n");
		return -1;
	}
	int ret, burstCount;
	TiXmlElement* pnode;
	Mat tmp;
	ret = getBurstCount(&burstCount);
	if (return_check(ret, "getBurstCount()", error_head)) return -1;
	lastValidSample.create(burstCount, 1, CV_32S);
	ret = xmldoc.find_node("burst", pnode);
	for (int i = 0; i < burstCount; i++)
	{
		if (!pnode) break;
		xmldoc.getIntArray("lastValidSample", tmp, pnode);
		for (int j = tmp.cols - 1; j >= 0; j--)
		{
			if (tmp.at<int>(0, j) != -1)
			{
				lastValidSample.at<int>(i, 0) = tmp.at<int>(0, j);
				break;
			}
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getFirstValidLine()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getFirstValidLine(): input check failed!\n");
		return -1;
	}
	int ret, burstCount;
	TiXmlElement* pnode;
	Mat tmp;
	ret = getBurstCount(&burstCount);
	if (return_check(ret, "getBurstCount()", error_head)) return -1;
	firstValidLine.create(burstCount, 1, CV_32S);
	ret = xmldoc.find_node("burst", pnode);
	for (int i = 0; i < burstCount; i++)
	{
		if (!pnode) break;
		xmldoc.getIntArray("firstValidSample", tmp, pnode);
		for (int j = 0; j < tmp.cols; j++)
		{
			if (tmp.at<int>(0, j) != -1)
			{
				firstValidLine.at<int>(i, 0) = j + 1;
				break;
			}
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getLastValidLine()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getLastValidLine(): input check failed!\n");
		return -1;
	}
	int ret, burstCount;
	TiXmlElement* pnode;
	Mat tmp;
	ret = getBurstCount(&burstCount);
	if (return_check(ret, "getBurstCount()", error_head)) return -1;
	lastValidLine.create(burstCount, 1, CV_32S);
	ret = xmldoc.find_node("burst", pnode);
	for (int i = 0; i < burstCount; i++)
	{
		if (!pnode) break;
		xmldoc.getIntArray("lastValidSample", tmp, pnode);
		for (int j = tmp.cols - 1; j >= 0; j--)
		{
			if (tmp.at<int>(0, j) != -1)
			{
				lastValidLine.at<int>(i, 0) = j + 1;
				break;
			}
		}
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getGeolocationGridPoint()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getGeolocationGridPoint(): input check failed!\n");
		return -1;
	}
	/*
	* 确定控制点个数
	*/
	int n_gcps = 1;
	int ret;
	TiXmlElement* pnode = NULL;
	ret = xmldoc.find_node("geolocationGridPointList", pnode);
	if (ret < 0)
	{
		fprintf(stderr, "Sentinel1Reader::getGeolocationGridPoint(): node geolocationGridPointList not found!\n");
		return -1;
	}
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &n_gcps);
	geolocationGridPoint.create(n_gcps, 7, CV_64F);//最后一行是方位向时间
	TiXmlElement* pchild = NULL;
	ret = xmldoc._find_node(pnode, "geolocationGridPoint", pchild);
	double lon, lat, height, row, col, inc, azimuthTime;
	for (int i = 0; i < n_gcps; i++)
	{
		if (!pchild) break;

		//longtitude
		ret = xmldoc._find_node(pchild, "longitude", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &lon);

		//latitude
		ret = xmldoc._find_node(pchild, "latitude", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &lat);

		//height
		ret = xmldoc._find_node(pchild, "height", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &height);

		//row
		ret = xmldoc._find_node(pchild, "line", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &row);

		//col
		ret = xmldoc._find_node(pchild, "pixel", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &col);

		//incidence angle
		ret = xmldoc._find_node(pchild, "incidenceAngle", pnode);
		ret = sscanf(pnode->GetText(), "%lf", &inc);

		//azimuthTime
		ret = xmldoc._find_node(pchild, "azimuthTime", pnode);
		UTC2GPS(pnode->GetText(), &azimuthTime);

		//assignment
		geolocationGridPoint.at<double>(i, 0) = lon;
		geolocationGridPoint.at<double>(i, 1) = lat;
		geolocationGridPoint.at<double>(i, 2) = height;
		geolocationGridPoint.at<double>(i, 3) = row;
		geolocationGridPoint.at<double>(i, 4) = col;
		geolocationGridPoint.at<double>(i, 5) = inc;
		geolocationGridPoint.at<double>(i, 6) = azimuthTime;
		pchild = pchild->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::updateGeolocationGridPoint()
{
	if (geolocationGridPoint.empty()) return 0;
	double startTime = geolocationGridPoint.at<double>(0, 6);
	for (int i = 0; i < geolocationGridPoint.rows; i++)
	{
		double time = geolocationGridPoint.at<double>(i, 6);
		
		if (fabs(time - startTime) > azimuthTimeInterval)
		{
			geolocationGridPoint.at<double>(i, 3) = (time - startTime) / azimuthTimeInterval;
		}
	}

	return 0;
}

int Sentinel1Reader::fitCoordinateConversionCoefficient()
{
	double mean_lon, mean_lat, mean_inc, max_lon, max_lat, max_inc, min_lon, min_lat, min_inc;
	Mat lon, lat, inc, row, col, gcps;
	geolocationGridPoint(cv::Range(0, geolocationGridPoint.rows), cv::Range(0, 6)).copyTo(gcps);
	gcps(cv::Range(0, gcps.rows), cv::Range(0, 1)).copyTo(lon);
	gcps(cv::Range(0, gcps.rows), cv::Range(1, 2)).copyTo(lat);
	gcps(cv::Range(0, gcps.rows), cv::Range(3, 4)).copyTo(row);
	gcps(cv::Range(0, gcps.rows), cv::Range(4, 5)).copyTo(col);
	gcps(cv::Range(0, gcps.rows), cv::Range(5, 6)).copyTo(inc);
	mean_lon = cv::mean(lon)[0];
	mean_lat = cv::mean(lat)[0];
	mean_inc = cv::mean(inc)[0];
	cv::minMaxLoc(lon, &min_lon, &max_lon);
	cv::minMaxLoc(lat, &min_lat, &max_lat);
	cv::minMaxLoc(inc, &min_inc, &max_inc);
	lon = (lon - mean_lon) / (max_lon - min_lon + 1e-10);
	lat = (lat - mean_lat) / (max_lat - min_lat + 1e-10);
	inc = (inc - mean_inc) / (max_inc - min_inc + 1e-10);
	row = (row + 1 - double(numberOfSamples) * 0.5) / (double(numberOfSamples) + 1e-10);//sentinel行列起点为0，+1统一为1.
	col = (col + 1 - double(numberOfSamples) * 0.5) / (double(numberOfSamples) + 1e-10);

	//拟合经度

	Mat A, B, b, temp, coefficient, error, eye, b_t, a, a_t;
	double rms;
	lon.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lon;
		temp.at<double>(0, 1) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 4) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 5) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		temp.copyTo(lon_coefficient);
	}

	//拟合纬度

	lat.copyTo(b);

	A = Mat::ones(lon.rows, 25, CV_64F);
	row.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = row.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	col.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	col.copyTo(temp);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(row);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = mean_lat;
		temp.at<double>(0, 1) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 4) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 5) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		temp.copyTo(lat_coefficient);
	}

	//拟合下视角

	inc.copyTo(b);
	A = Mat::ones(inc.rows, 6, CV_64F);
	col.copyTo(A(cv::Range(0, inc.rows), cv::Range(1, 2)));
	temp = col.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(2, 3)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(3, 4)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(4, 5)));
	temp = temp.mul(col);
	temp.copyTo(A(cv::Range(0, inc.rows), cv::Range(5, 6)));
	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 11, CV_64F);
		temp.at<double>(0, 0) = mean_inc;
		temp.at<double>(0, 1) = max_inc - min_inc + 1e-10;
		temp.at<double>(0, 2) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 3) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 10) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(4, 10)));
		temp.copyTo(inc_coefficient);
	}

	//拟合行坐标

	row.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 1) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		temp.copyTo(row_coefficient);
	}

	//拟合列坐标

	col.copyTo(b);
	A = Mat::ones(lon.rows, 25, CV_64F);
	lon.copyTo(A(cv::Range(0, lon.rows), cv::Range(1, 2)));
	temp = lon.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(2, 3)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(3, 4)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(4, 5)));

	lat.copyTo(temp);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(5, 6)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(6, 7)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(7, 8)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(8, 9)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(9, 10)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(10, 11)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(11, 12)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(12, 13)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(13, 14)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(14, 15)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(15, 16)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(16, 17)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(17, 18)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(18, 19)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(19, 20)));

	lat.copyTo(temp);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp = temp.mul(lat);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(20, 21)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(21, 22)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(22, 23)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(23, 24)));
	temp = temp.mul(lon);
	temp.copyTo(A(cv::Range(0, lon.rows), cv::Range(24, 25)));

	cv::transpose(A, temp);
	B = temp * b;
	A.copyTo(a);
	cv::transpose(a, a_t);
	A = temp * A;
	rms = -1.0;
	if (cv::invert(A, error, cv::DECOMP_LU) > 0)
	{
		cv::transpose(b, b_t);
		error = b_t * b - (b_t * a) * error * (a_t * b);
		//error = b_t * (eye - a * error * a_t) * b;
		rms = sqrt(error.at<double>(0, 0) / double(b.rows));
	}
	if (cv::solve(A, B, coefficient, cv::DECOMP_NORMAL))
	{
		temp.create(1, 32, CV_64F);
		temp.at<double>(0, 0) = double(numberOfSamples) * 0.5;
		temp.at<double>(0, 1) = double(numberOfSamples) + 1e-10;
		temp.at<double>(0, 2) = mean_lon;
		temp.at<double>(0, 3) = max_lon - min_lon + 1e-10;
		temp.at<double>(0, 4) = mean_lat;
		temp.at<double>(0, 5) = max_lat - min_lat + 1e-10;
		temp.at<double>(0, 31) = rms;
		cv::transpose(coefficient, coefficient);
		coefficient.copyTo(temp(cv::Range(0, 1), cv::Range(6, 31)));
		temp.copyTo(col_coefficient);
	}
	return 0;
}

int Sentinel1Reader::getOrbitList()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getOrbitList(): input check failed!\n");
		return -1;
	}
	TiXmlElement* pnode, * pchild, * pchild1;
	int ret, numOfstateVec;
	ret = xmldoc.find_node("orbitList", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &numOfstateVec);

	ret = xmldoc.find_node("orbit", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	double time, x, y, z, vx, vy, vz;
	orbitList.create(numOfstateVec, 7, CV_64F);
	for (int i = 0; i < numOfstateVec; i++)
	{
		if (!pnode) break;
		//GPS时间
		ret = xmldoc._find_node(pnode, "time", pchild);
		ret = UTC2GPS(pchild->GetText(), &time);
		//位置x
		ret = xmldoc._find_node(pnode, "position", pchild);
		ret = xmldoc._find_node(pchild, "x", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &x);
		//位置y
		ret = xmldoc._find_node(pchild, "y", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &y);
		//位置z
		ret = xmldoc._find_node(pchild, "z", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &z);
		//速度x
		ret = xmldoc._find_node(pnode, "velocity", pchild);
		ret = xmldoc._find_node(pchild, "x", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &vx);
		//速度y
		ret = xmldoc._find_node(pchild, "y", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &vy);
		//速度z
		ret = xmldoc._find_node(pchild, "z", pchild1);
		ret = sscanf(pchild1->GetText(), "%lf", &vz);

		//赋值
		orbitList.at<double>(i, 0) = time;
		orbitList.at<double>(i, 1) = x;
		orbitList.at<double>(i, 2) = y;
		orbitList.at<double>(i, 3) = z;
		orbitList.at<double>(i, 4) = vx;
		orbitList.at<double>(i, 5) = vy;
		orbitList.at<double>(i, 6) = vz;
		pnode = pnode->NextSiblingElement();
	}
	return 0;
}

int Sentinel1Reader::getPOD(const char* POD_file)
{
	if (!bXmlLoad || !POD_file)
	{
		fprintf(stderr, "Sentinel1Reader::getPOD(): input check failed!\n");
		return -1;
	}
	/*
	* 读取精密轨道数据
	*/

	int ret, numOfstateVec;
	double start_time, stop_time;
	string start_time_str, stop_time_str;
	ret = xmldoc.get_str_para("startTime", start_time_str);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = xmldoc.get_str_para("stopTime", stop_time_str);
	UTC2GPS(start_time_str.c_str(), &start_time);
	UTC2GPS(stop_time_str.c_str(), &stop_time);
	TiXmlElement* pnode, * pchild;
	XMLFile doc;
	ret = doc.XMLFile_load(POD_file);
	if (return_check(ret, "XMLFile_load()", error_head)) return -1;
	ret = doc.find_node("List_of_OSVs", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	ret = sscanf(pnode->FirstAttribute()->Value(), "%d", &numOfstateVec);


	Mat tmp = Mat::zeros(numOfstateVec, 7, CV_64F);
	ret = doc.find_node("OSV", pnode);
	if (return_check(ret, "find_node()", error_head)) return -1;
	string str;
	double gps_time, x, y, z, vx, vy, vz;
	bool start = false; bool stop = false;
	int count = 0;
	for (int i = 0; i < numOfstateVec; i++)
	{
		if (!pnode || stop) break;
		ret = doc._find_node(pnode, "UTC", pchild);
		str = pchild->GetText();
		str = str.substr(4);
		ret = UTC2GPS(str.c_str(), &gps_time);
		if (gps_time <= start_time && fabs(gps_time - start_time) <= 100.0) start = true;
		if (gps_time >= stop_time && fabs(gps_time - stop_time) >= 100.0) stop = true;

		if (start && !stop)//开始记录
		{
			ret = doc._find_node(pnode, "X", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &x);
			ret = doc._find_node(pnode, "Y", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &y);
			ret = doc._find_node(pnode, "Z", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &z);
			ret = doc._find_node(pnode, "VX", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &vx);
			ret = doc._find_node(pnode, "VY", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &vy);
			ret = doc._find_node(pnode, "VZ", pchild);
			ret = sscanf(pchild->GetText(), "%lf", &vz);

			tmp.at<double>(count, 0) = gps_time;
			tmp.at<double>(count, 1) = x;
			tmp.at<double>(count, 2) = y;
			tmp.at<double>(count, 3) = z;
			tmp.at<double>(count, 4) = vx;
			tmp.at<double>(count, 5) = vy;
			tmp.at<double>(count, 6) = vz;
			count++;
		}

		pnode = pnode->NextSiblingElement();
	}
	if (count < 1)
	{
		fprintf(stderr, "Sentinel1Reader::getPOD(): orbit mismatch! please check if POD file!\n");
		return -1;
	}
	tmp(cv::Range(0, count), cv::Range(0, 7)).copyTo(preciseOrbitList);
	return 0;
}

int Sentinel1Reader::getOtherParameters()
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getOtherParameters(): input check failed!\n");
		return -1;
	}
	int ret = xmldoc.get_double_para("azimuthPixelSpacing", &this->azimuthPixelSpacing);
	if (return_check(ret, "get_double_para()", error_head)) return -1;
	ret = xmldoc.get_double_para("azimuthSteeringRate", &this->azimuthSteeringRate);
	ret = xmldoc.get_double_para("azimuthTimeInterval", &this->azimuthTimeInterval);
	ret = xmldoc.get_double_para("platformHeading", &this->headingAngle);
	ret = xmldoc.get_double_para("radarFrequency", &this->radarFrequency);
	ret = xmldoc.get_double_para("rangePixelSpacing", &this->rangePixelSpacing);
	ret = xmldoc.get_double_para("rangeSamplingRate", &this->rangeSamplingRate);
	ret = xmldoc.get_double_para("slantRangeTime", &this->slantRangeTime);
	ret = xmldoc.get_double_para("incidenceAngleMidSwath", &this->incidence_center);

	ret = xmldoc.get_int_para("numberOfLines", &this->numberOfLines);
	ret = xmldoc.get_int_para("numberOfSamples", &this->numberOfSamples);
	ret = xmldoc.get_int_para("linesPerBurst", &this->linesPerBurst);
	ret = getBurstCount(&this->burstCount);

	ret = xmldoc.get_str_para("polarisation", this->polarization);
	ret = xmldoc.get_str_para("swath", this->swath);
	ret = xmldoc.get_str_para("pass", this->pass);
	xmldoc.get_str_para("startTime", this->startTime);
	xmldoc.get_str_para("stopTime", this->stopTime);
	return 0;
}

int Sentinel1Reader::prepareData(const char* PODFile)
{
	if (!bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::prepareData(): input check failed!\n");
		return -1;
	}
	int ret;
	ret = getOtherParameters();
	if (return_check(ret, "getOtherParameters()", error_head)) return -1;
	this->getAntennaPattern();
	this->getAzimuthFmRateList();
	this->getBurstAzimuthTime();
	this->getDcEstimateList();
	this->getFirstValidLine();
	this->getFirstValidSample();
	this->getGeolocationGridPoint();
	this->getLastValidLine();
	this->getLastValidSample();
	this->getOrbitList();
	if (PODFile)
	{
		this->getPOD(PODFile);
	}
	//更新控制点
	updateGeolocationGridPoint();
	//拟合坐标转换系数
	fitCoordinateConversionCoefficient();
	isDataAvailable = true;
	return 0;
}

int Sentinel1Reader::getSLC(ComplexMat& slc)
{
	int ret;
	if (this->tiffFile.empty() || !bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::getSLC(): input check failed!\n");
		return -1;
	}
	FILE* fp = fopen(this->tiffFile.c_str(), "rb");
	if (!fp) {
		fprintf(stderr, "getSLC(): can't open %s\n", tiffFile.c_str());
		return -1;
	}
	int byteOffset;
	ret = xmldoc.get_int_para("byteOffset", &byteOffset);
	if (return_check(ret, "get_int_para()", error_head)) {  
		if (fp) fclose(fp);
		return -1;
	}
	size_t size = numberOfLines * numberOfSamples * 2 * sizeof(short);
	short* buf = (short*)malloc(size);
	if (!buf) {
		fprintf(stderr, "getSLC(): out of memory!\n");
		if (fp) fclose(fp);
		return -1;
	}
	fseek(fp, byteOffset, SEEK_SET);
	fread(buf, sizeof(short), numberOfLines * numberOfSamples * 2, fp);
	if (fp) fclose(fp);
	size_t offset = 0;
	slc.re.create(numberOfLines, numberOfSamples, CV_16S);
	slc.im.create(numberOfLines, numberOfSamples, CV_16S);
	for (int j = 0; j < numberOfLines; j++)
	{
		for (int k = 0; k < numberOfSamples; k++)
		{
			slc.re.at<short>(j, k) = buf[offset];
			offset++;
			slc.im.at<short>(j, k) = buf[offset];
			offset++;
		}
	}
	if (buf) free(buf);
	return 0;
}

int Sentinel1Reader::writeToh5(const char* h5File)
{
	int ret;
	if (!h5File || !bXmlLoad)
	{
		fprintf(stderr, "Sentinel1Reader::writeToh5(): input check failed!\n");
		return -1;
	}
	if (!isDataAvailable)
	{
		if (PODFile.empty()) ret = prepareData();
		else ret = prepareData(PODFile.c_str());
		if (return_check(ret, "prepareData()", error_head)) return -1;
	}
	
	FormatConversion conversion;
	ret = conversion.creat_new_h5(h5File);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_str_to_h5(h5File, "polarization", this->polarization.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	conversion.write_str_to_h5(h5File, "orbit_dir", pass.c_str());
	conversion.write_str_to_h5(h5File, "swath", swath.c_str());
	conversion.write_str_to_h5(h5File, "imaging_mode", "TOPS");
	conversion.write_str_to_h5(h5File, "sensor", sensor.c_str());
	conversion.write_str_to_h5(h5File, "acquisition_start_time", startTime.c_str());
	conversion.write_str_to_h5(h5File, "acquisition_stop_time", stopTime.c_str());

	conversion.write_double_to_h5(h5File, "azimuth_spacing", this->azimuthPixelSpacing);
	conversion.write_double_to_h5(h5File, "range_spacing", this->rangePixelSpacing);
	conversion.write_double_to_h5(h5File, "azimuthSteeringRate", this->azimuthSteeringRate);
	conversion.write_double_to_h5(h5File, "prf", 1.0 / this->azimuthTimeInterval);
	conversion.write_double_to_h5(h5File, "heading", this->headingAngle);
	conversion.write_double_to_h5(h5File, "incidence_center", this->incidence_center);
	conversion.write_double_to_h5(h5File, "carrier_frequency", this->radarFrequency);
	conversion.write_double_to_h5(h5File, "slant_range_first_pixel", this->slantRangeTime * VEL_C / 2.0);

	conversion.write_int_to_h5(h5File, "burstCount", this->burstCount);
	conversion.write_int_to_h5(h5File, "linesPerBurst", this->linesPerBurst);
	conversion.write_int_to_h5(h5File, "samplesPerBurst", this->numberOfSamples);
	conversion.write_int_to_h5(h5File, "range_len", this->numberOfSamples);
	conversion.write_int_to_h5(h5File, "azimuth_len", this->numberOfLines);

	conversion.write_array_to_h5(h5File, "antennaPattern_elevationAngle", this->antennaPattern_elevationAngle);
	conversion.write_array_to_h5(h5File, "antennaPattern_slantRangeTime", this->antennaPattern_slantRangeTime);
	conversion.write_array_to_h5(h5File, "azimuthFmRateList", this->AzimuthFmRateList);
	conversion.write_array_to_h5(h5File, "burstAzimuthTime", this->burstAzimuthTime);
	conversion.write_array_to_h5(h5File, "dcEstimateList", this->DcEstimateList);
	conversion.write_array_to_h5(h5File, "firstValidLine", this->firstValidLine);
	conversion.write_array_to_h5(h5File, "firstValidSample", this->firstValidSample);
	conversion.write_array_to_h5(h5File, "lon_coefficient", this->lon_coefficient);
	conversion.write_array_to_h5(h5File, "lat_coefficient", this->lat_coefficient);
	conversion.write_array_to_h5(h5File, "row_coefficient", this->row_coefficient);
	conversion.write_array_to_h5(h5File, "col_coefficient", this->col_coefficient);
	conversion.write_array_to_h5(h5File, "inc_coefficient", this->inc_coefficient);


	Mat gcps;
	if (geolocationGridPoint.cols > 6) geolocationGridPoint(cv::Range(0, geolocationGridPoint.rows), cv::Range(0, 6)).copyTo(gcps);
	else geolocationGridPoint.copyTo(gcps);
	conversion.write_array_to_h5(h5File, "gcps", gcps);
	conversion.write_array_to_h5(h5File, "lastValidLine", this->lastValidLine);
	conversion.write_array_to_h5(h5File, "lastValidSample", this->lastValidSample);
	conversion.write_array_to_h5(h5File, "state_vec", this->orbitList);
	if(!preciseOrbitList.empty())
		conversion.write_array_to_h5(h5File, "fine_state_vec", this->preciseOrbitList);

	//写入图像数据
	ComplexMat slc;
	ret = getSLC(slc);
	if (return_check(ret, "getSLC()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(h5File, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	return 0;
}

/*------------------------------------------------*/
/*             哨兵一号数据/计算读取工具          */
/*------------------------------------------------*/

Sentinel1Utils::Sentinel1Utils(const char* h5File)
{
	bInitialized = false;
	memset(m_xmlFileName, 0, 2048);
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	this->azimuthPixelSpacing = 0.0;
	this->azimuthSteeringRate = 0.0;
	this->azimuthTimeInterval = 0.0;
	this->numberOfLines = 0;
	this->headingAngle = 0.0;
	this->numberOfSamples = 0;
	this->samplesPerBurst = 0;
	this->burstCount = 0;
	this->linesPerBurst = 0;
	this->pass = "";
	this->polarization = "";
	this->radarFrequency = 0.0;
	this->rangePixelSpacing = 0.0;
	this->rangeSamplingRate = 0.0;
	this->sensor = "sentinel";
	this->swath = "";
	this->slantRangeTime = 0.0;
	this->h5File = h5File;
	stateVectors = NULL;

	this->isDopplerCentroidAvailable = false;
	this->isDopplerRateAvailable = false;
	this->isRangeDependDopplerRateAvailiable = false;
	this->isReferenceTimeAvailable = false;

	this->burstOffset = -9999;

}

Sentinel1Utils::~Sentinel1Utils()
{
	if (stateVectors)
	{
		delete stateVectors;
		stateVectors = NULL;
	}
}

int Sentinel1Utils::init()
{
	if (h5File.empty())
	{
		fprintf(stderr, "init(): input check failed!\n");
		return -1;
	}
	int ret; FormatConversion conversion;
	//读取数据
	ret = conversion.read_double_from_h5(h5File.c_str(), "azimuth_spacing", &this->azimuthPixelSpacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5File.c_str(), "azimuthSteeringRate", &this->azimuthSteeringRate);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5File.c_str(), "prf", &this->azimuthTimeInterval);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	azimuthTimeInterval = 1.0 / azimuthTimeInterval;
	ret = conversion.read_double_from_h5(h5File.c_str(), "carrier_frequency", &this->radarFrequency);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5File.c_str(), "range_spacing", &this->rangePixelSpacing);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(h5File.c_str(), "slant_range_first_pixel", &this->slantRangeTime);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	slantRangeTime = 2.0 * slantRangeTime / VEL_C;

	ret = conversion.read_int_from_h5(h5File.c_str(), "burstCount", &this->burstCount);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(h5File.c_str(), "linesPerBurst", &this->linesPerBurst);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(h5File.c_str(), "range_len", &this->numberOfSamples);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	this->samplesPerBurst = this->numberOfSamples;
	ret = conversion.read_int_from_h5(h5File.c_str(), "azimuth_len", &this->numberOfLines);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_array_from_h5(h5File.c_str(), "azimuthFmRateList", this->AzimuthFmRateList);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "burstAzimuthTime", this->burstAzimuthTime);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "dcEstimateList", this->DcEstimateList);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "firstValidSample", this->firstValidSample);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "firstValidLine", this->firstValidLine);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "lastValidLine", this->lastValidLine);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "lastValidSample", this->lastValidSample);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "state_vec", this->orbitList);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "fine_state_vec", this->preciseOrbitList);
	ret = conversion.read_array_from_h5(h5File.c_str(), "antennaPattern_elevationAngle", this->antennaPattern_elevationAngle);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "antennaPattern_slantRangeTime", this->antennaPattern_slantRangeTime);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(h5File.c_str(), "gcps", this->geolocationGridPoint);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	string str;
	double startTime, stopTime;
	ret = conversion.read_str_from_h5(h5File.c_str(), "acquisition_start_time", str);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	UTC2GPS(str.c_str(), &startTime);
	ret = conversion.read_str_from_h5(h5File.c_str(), "acquisition_stop_time", str);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	UTC2GPS(str.c_str(), &stopTime);

	if (!preciseOrbitList.empty())
	{
		this->stateVectors = new orbitStateVectors(preciseOrbitList, startTime, stopTime);
	}
	else
	{
		this->stateVectors = new orbitStateVectors(orbitList, startTime, stopTime);
	}
	ret = this->stateVectors->applyOrbit();
	if (return_check(ret, "applyOrbit()", error_head)) return -1;
	bInitialized = true;
	computeDopplerCentroid();
	return 0;
}

int Sentinel1Utils::computeReferenceTime()
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::computeReferenceTime(): input check failed!\n");
		return -1;
	}
	int ret;
	if (!isDopplerCentroidAvailable) {
		ret = computeDopplerCentroid();
		if (return_check(ret, "computeDopplerCentroid()", error_head)) return -1;
	}

	if (!isRangeDependDopplerRateAvailiable) {
		ret = computeRangeDependDopplerRate();
		if (return_check(ret, "computeRangeDependDopplerRate()", error_head)) return -1;
	}
	double tmp;
	tmp = (double)linesPerBurst * azimuthTimeInterval / 2.0;
	referenceTime.create(burstCount, samplesPerBurst, CV_64F);
	for (int i = 0; i < burstCount; i++)
	{
		double tmp2 = tmp + dopplerCentroid.at<double>(i, firstValidSample.at<int>(i, 0)) /
			rangeDependDopplerRate.at<double>(i, firstValidSample.at<int>(i, 0));
		for (int j = 0; j < samplesPerBurst; j++)
		{
			referenceTime.at<double>(i, j) = tmp2 - dopplerCentroid.at<double>(i, j) / rangeDependDopplerRate.at<double>(i, j);
		}
	}
	isReferenceTimeAvailable = true;
	return 0;
}

int Sentinel1Utils::computeRangeDependDopplerRate()
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::computeRangeDependDopplerRate(): input check failed!\n");
		return -1;
	}
	int ret;
	rangeDependDopplerRate.create(burstCount, samplesPerBurst, CV_64F);
	for (int i = 0; i < burstCount; i++)
	{
		for (int j = 0; j < samplesPerBurst; j++)
		{
			double slrt = 2 * (slantRangeTime / 2 + j * rangePixelSpacing / VEL_C);
			double dt; int k;
			if (burstAzimuthTime.at<double>(i, 0) <= AzimuthFmRateList.at<double>(0, 0))
			{
				k = 0;
				dt = slrt - AzimuthFmRateList.at<double>(k, 1);
			}
			else if (burstAzimuthTime.at<double>(i, 0) > AzimuthFmRateList.at<double>(AzimuthFmRateList.rows - 1, 0))
			{
				k = AzimuthFmRateList.rows - 1;
				dt = slrt - AzimuthFmRateList.at<double>(k, 1);
				
			}
			else
			{
				for (k = 1; k < AzimuthFmRateList.rows; k++)
				{
					if (AzimuthFmRateList.at<double>(k, 0) >= burstAzimuthTime.at<double>(i, 0) &&
						AzimuthFmRateList.at<double>(k - 1, 0) < burstAzimuthTime.at<double>(i, 0))
					{
						dt = slrt - AzimuthFmRateList.at<double>(k, 1);
						break;
					}
				}
			}
			double c0, c1, c2;
			c0 = AzimuthFmRateList.at<double>(k, 2);
			c1 = AzimuthFmRateList.at<double>(k, 3);
			c2 = AzimuthFmRateList.at<double>(k, 4);
			rangeDependDopplerRate.at<double>(i, j) = c0 + c1 * dt + c2 * dt * dt;
			
		}
	}
	isRangeDependDopplerRateAvailiable = true;
	return 0;
}

int Sentinel1Utils::computeDopplerCentroid()
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::computeDopplerCentroid(): input check failed!\n");
		return -1;
	}
	int ret;
	dopplerCentroid.create(burstCount, samplesPerBurst, CV_64F);
	for (int i = 0; i < burstCount; i++)
	{
		for (int j = 0; j < samplesPerBurst; j++)
		{
			double slrt = 2 * (slantRangeTime / 2 + j * rangePixelSpacing / VEL_C);
			double dt; int k;
			if (burstAzimuthTime.at<double>(i, 0) <= DcEstimateList.at<double>(0, 0))
			{
				k = 0;
				dt = slrt - DcEstimateList.at<double>(k, 1);
			}
			else if (burstAzimuthTime.at<double>(i, 0) > DcEstimateList.at<double>(DcEstimateList.rows - 1, 0))
			{
				k = DcEstimateList.rows - 1;
				dt = slrt - DcEstimateList.at<double>(k, 1);

			}
			else
			{
				for (k = 1; k < DcEstimateList.rows; k++)
				{
					if (DcEstimateList.at<double>(k, 0) >= burstAzimuthTime.at<double>(i, 0) &&
						DcEstimateList.at<double>(k - 1, 0) < burstAzimuthTime.at<double>(i, 0))
					{
						dt = slrt - DcEstimateList.at<double>(k, 1);
						break;
					}
				}
			}
			double c0, c1, c2;
			c0 = DcEstimateList.at<double>(k, 2);
			c1 = DcEstimateList.at<double>(k, 3);
			c2 = DcEstimateList.at<double>(k, 4);
			dopplerCentroid.at<double>(i, j) = c0 + c1 * dt + c2 * dt * dt;

		}
	}
	isDopplerCentroidAvailable = true;
	return 0;
}

int Sentinel1Utils::computeDopplerRate()
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::computeDopplerRate(): input check failed!\n");
		return -1;
	}
	int ret;
	if (!isRangeDependDopplerRateAvailiable)
	{
		ret = computeRangeDependDopplerRate();
		if (return_check(ret, "computeRangeDependDopplerRate()", error_head)) return -1;
	}
	double waveLength = VEL_C / radarFrequency;
	dopplerRate.create(burstCount, samplesPerBurst, CV_64F);
	for (int i = 0; i < burstCount; i++)
	{
		double v = sqrt(orbitList.at<double>(0, 6) * orbitList.at<double>(0, 6) +
			orbitList.at<double>(0, 5) * orbitList.at<double>(0, 5) +
			orbitList.at<double>(0, 4) * orbitList.at<double>(0, 4));
		double krot = 2 * v * azimuthSteeringRate * PI / 180.0 / waveLength;
		for (int j = 0; j < samplesPerBurst; j++)
		{
			dopplerRate.at<double>(i, j) = rangeDependDopplerRate.at<double>(i, j) * krot /
				(rangeDependDopplerRate.at<double>(i, j) - krot);
		}
	}
	isDopplerRateAvailable = true;
	return 0;
}

int Sentinel1Utils::computeDerampDemodPhase(
	int burstIndex, 
	Mat& derampDemodPhase
)
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::computeDerampDemodPhase(): input check failed!\n");
		return -1;
	}
	if (burstIndex < 1 || burstIndex > burstCount)
	{
		fprintf(stderr, "Sentinel1Utils::computeDerampDemodPhase(): input check failed!\n");
		return -1;
	}
	int ret;
	if (!isDopplerRateAvailable)
	{
		ret = computeDopplerRate();
		if (return_check(ret, "computeDopplerRate()", error_head)) return -1;
	}
	if (!isDopplerCentroidAvailable)
	{
		ret = computeDopplerCentroid();
		if (return_check(ret, "computeDopplerCentroid()", error_head)) return -1;
	}
	if (!isReferenceTimeAvailable)
	{
		ret = computeReferenceTime();
		if (return_check(ret, "computeReferenceTime()", error_head)) return -1;
	}
	derampDemodPhase.create(linesPerBurst, samplesPerBurst, CV_64F);
	int firstLineInBurst = (burstIndex - 1) * linesPerBurst;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < linesPerBurst; i++)
	{
		//double ta = double(i - firstLineInBurst) * azimuthTimeInterval;
		double ta = (double)i * azimuthTimeInterval;
		for (int j = 0; j < samplesPerBurst; j++)
		{
			double kt = dopplerRate.at<double>(burstIndex - 1, j);
			double deramp = -PI * kt * pow(ta - referenceTime.at<double>(burstIndex - 1, j), 2.0);
			double demod = -2 * PI * ta * dopplerCentroid.at<double>(burstIndex - 1, j);
			derampDemodPhase.at<double>(i, j) = deramp + demod;
		}
	}
	return 0;
}

int Sentinel1Utils::getBurst(int burstIndex, ComplexMat& burstSLC)
{
	if (!bInitialized)
	{
		fprintf(stderr, "Sentinel1Utils::getBurst(): input check failed!\n");
		return -1;
	}
	if (burstIndex < 1 || burstIndex > burstCount)
	{
		fprintf(stderr, "Sentinel1Utils::getBurst(): burstIndex out of legal range!\n");
		return -1;
	}
	int ret;
	FormatConversion conversion;
	ret = conversion.read_subarray_from_h5(this->h5File.c_str(), "s_re", linesPerBurst * (burstIndex - 1),
		0, linesPerBurst, samplesPerBurst, burstSLC.re);
	if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
	ret = conversion.read_subarray_from_h5(this->h5File.c_str(), "s_im", linesPerBurst * (burstIndex - 1),
		0, linesPerBurst, samplesPerBurst, burstSLC.im);
	if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;

	return 0;
}

int Sentinel1Utils::getDopplerFrequency(
	Position groundPosition,
	Position satellitePosition, 
	Velocity satelliteVelocity,
	double* dopplerFrequency
)
{
	if (!bInitialized || !dopplerFrequency)
	{
		fprintf(stderr, "getDopplerFrequency(): input check failed!");
		return -1;
	}
	double waveLength = VEL_C / radarFrequency;
	double xdiff = groundPosition.x - satellitePosition.x;
	double ydiff = groundPosition.y - satellitePosition.y;
	double zdiff = groundPosition.z - satellitePosition.z;
	double distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
	*dopplerFrequency = 2.0 * (xdiff * satelliteVelocity.vx + ydiff * satelliteVelocity.vy + zdiff * satelliteVelocity.vz) / (waveLength * distance);
	return 0;
}

int Sentinel1Utils::getZeroDopplerTime(Position groundPosition, double* zeroDopplerTime, double dopplerFrequency)
{
	int ret;
	if (!bInitialized || !zeroDopplerTime)
	{
		fprintf(stderr, "getZeroDopplerTime(): input check failed!");
		return -1;
	}

	int numOrbitVec = stateVectors->newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double wavelength = VEL_C / radarFrequency;
	for (int i = 0; i < numOrbitVec; i++) {
		Position orb_pos(stateVectors->newStateVectors.at<double>(i, 1), stateVectors->newStateVectors.at<double>(i, 2),
			stateVectors->newStateVectors.at<double>(i, 3));
		Velocity orb_vel(stateVectors->newStateVectors.at<double>(i, 4), stateVectors->newStateVectors.at<double>(i, 5),
			stateVectors->newStateVectors.at<double>(i, 6));
		double currentFreq = 0;
		getDopplerFrequency(groundPosition, orb_pos, orb_vel, &currentFreq);
		//currentFreq = currentFreq - dopplerFrequency;
		if (i == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors->newStateVectors.at<double>(i, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors->newStateVectors.at<double>(i, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "getZeroDopplerTime(): zeroDopplerTime out of legal range!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = azimuthTimeInterval;

	int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	int numIterations = 0; Position pos; Velocity vel;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors->getPosition(midTime, pos);
		stateVectors->getVelocity(midTime, vel);
		getDopplerFrequency(groundPosition, pos, vel, &midFreq);
		//midFreq = midFreq - dopplerFrequency;
		if ((midFreq - dopplerFrequency) * (lowerBoundFreq - dopplerFrequency) > 0.0) {
			lowerBoundTime = midTime;
			lowerBoundFreq = midFreq;
		}
		else if ((midFreq - dopplerFrequency) * (upperBoundFreq - dopplerFrequency) > 0.0) {
			upperBoundTime = midTime;
			upperBoundFreq = midFreq;
		}
		else if (fabs(midFreq - dopplerFrequency) < 0.01) {
			*zeroDopplerTime =  midTime;
			return 0;
		}

		diffTime = fabs(upperBoundTime - lowerBoundTime);
		numIterations++;
	}


	*zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);

	return 0;
}

int Sentinel1Utils::getRgAzPosition(
	int burstIndex,
	Position groundPosition,
	double* rangeIndex,
	double* azimuthIndex
)
{
	int ret;
	if (!bInitialized || !rangeIndex || !azimuthIndex)
	{
		fprintf(stderr, "getRgAzPosition(): input check failed!");
		return -1;
	}
	double zeroDopplerTime, slantRange; 
	ret = getZeroDopplerTime(groundPosition, &zeroDopplerTime, dopplerCentroid.at<double>(burstIndex - 1, (int)samplesPerBurst / 2));
	if (return_check(ret, "getZeroDopplerTime()", error_head)) return -1;
	*azimuthIndex = (zeroDopplerTime - burstAzimuthTime.at<double>(burstIndex - 1)) / azimuthTimeInterval;
	ret = getSlantRange(zeroDopplerTime, groundPosition, &slantRange);
	if (return_check(ret, "getSlantRange()", error_head)) return -1;
	*rangeIndex = (slantRange - slantRangeTime * VEL_C * 0.5) / rangePixelSpacing;

	if (*azimuthIndex < 0.0 || *rangeIndex < 0.0 || *rangeIndex >= samplesPerBurst || *azimuthIndex >= linesPerBurst) return -1;
	int x = *rangeIndex - 1; x = x < 0 ? 0 : x;
	ret = getZeroDopplerTime(groundPosition, &zeroDopplerTime, dopplerCentroid.at<double>(burstIndex - 1, (int)x));
	if (return_check(ret, "getZeroDopplerTime()", error_head)) return -1;
	*azimuthIndex = (zeroDopplerTime - burstAzimuthTime.at<double>(burstIndex - 1)) / azimuthTimeInterval;
	ret = getSlantRange(zeroDopplerTime, groundPosition, &slantRange);
	if (return_check(ret, "getSlantRange()", error_head)) return -1;
	*rangeIndex = (slantRange - slantRangeTime * VEL_C * 0.5) / rangePixelSpacing;
	if (*azimuthIndex < 0.0 || *rangeIndex < 0.0 || *rangeIndex >= samplesPerBurst || *azimuthIndex >= linesPerBurst) return -1;

	return 0;
}

int Sentinel1Utils::getSlantRange(double azimuthTime, Position groundPosition, double* slantRange)
{
	int ret;
	if (!bInitialized || !slantRange)
	{
		fprintf(stderr, "getSlantRange(): input check failed!");
		return -1;
	}
	Position satellitePosition;
	stateVectors->getPosition(azimuthTime, satellitePosition);
	double xdiff = groundPosition.x - satellitePosition.x;
	double ydiff = groundPosition.y - satellitePosition.y;
	double zdiff = groundPosition.z - satellitePosition.z;
	*slantRange = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
	return 0;
}

int Sentinel1Utils::getBurstIndice(Position groundPosition, BurstIndices& burstIndice)
{
	if (!bInitialized)
	{
		fprintf(stderr, "getBurstIndice(): input check failed!");
		return -1;
	}
	int ret;
	double zeroDopplerTime;
	ret = getZeroDopplerTime(groundPosition, &zeroDopplerTime);
	if (return_check(ret, "getZeroDopplerTime()", error_head)) return -1;
	int k = 0;
	double burstFirstLineTime, burstLastLineTime;
	for (int i = 0; i < burstCount; i++) {
		burstFirstLineTime = burstAzimuthTime.at<double>(i, 0);
		burstLastLineTime = burstFirstLineTime + azimuthTimeInterval * (linesPerBurst - 1);
		if (zeroDopplerTime >= burstFirstLineTime && zeroDopplerTime < burstLastLineTime) {
			bool inUpperPartOfBurst = (zeroDopplerTime >= (burstFirstLineTime + burstLastLineTime) / 2.0);

			if (k == 0) {
				burstIndice.firstBurstIndex = i + 1;
				burstIndice.inUpperPartOfFirstBurst = inUpperPartOfBurst;
			}
			else {
				burstIndice.secondBurstIndex = i + 1;
				burstIndice.inUpperPartOfSecondBurst = inUpperPartOfBurst;
				break;
			}
			++k;
		}
	}
	if (k == 0) return -1;
	return 0;
}

int Sentinel1Utils::computeImageGeoBoundry(
	double* lonMin,
	double* lonMax,
	double* latMin,
	double* latMax
)
{
	if (!bInitialized || !lonMin || !lonMax || !latMin || !latMax)
	{
		fprintf(stderr, "computeImageGeoBoundry(): input check failed!");
		return -1;
	}
	*lonMin = 181.0;
	*lonMax = -181.0;
	*latMin = 91.0;
	*latMax = -91.0;
	for (int i = 0; i < geolocationGridPoint.rows; i++)
	{
		double lon = geolocationGridPoint.at<double>(i, 0);
		double lat = geolocationGridPoint.at<double>(i, 1);
		*lonMin = *lonMin > lon ? lon : *lonMin;
		*lonMax = *lonMax < lon ? lon : *lonMax;
		*latMin = *latMin > lat ? lat : *latMin;
		*latMax = *latMax < lat ? lat : *latMax;
	}
	double extra = 5.0 / 6000;
	*lonMin = *lonMin - extra * 100;
	*lonMax = *lonMax + extra * 100;
	*latMin = *latMin - extra * 100;
	*latMax = *latMax + extra * 100;
	return 0;
}

int Sentinel1Utils::deburst(const char* outFile)
{
	if (!bInitialized || !outFile)
	{
		fprintf(stderr, "deburst(): input check failed!\n");
		return -1;
	}

	//创建新的deburst文件

	FormatConversion conversion;
	int ret = conversion.creat_new_h5(outFile);
	Mat start(this->burstCount, 1, CV_32S), end(this->burstCount, 1, CV_32S);
	start.at<int>(0, 0) = 1;
	end.at<int>(0, 0) = this->lastValidLine.at<int>(0, 0);
	double lastValidTime = this->burstAzimuthTime.at<double>(0, 0) +
		(this->lastValidLine.at<int>(0, 0) - 1) * this->azimuthTimeInterval;
	double firstValidTime;
	int deburstLines = 0;
	int overlap;
	for (int i = 1; i < this->burstCount; i++)
	{
		firstValidTime = this->burstAzimuthTime.at<double>(i, 0) + (this->firstValidLine.at<int>(i, 0) - 1) *
			this->azimuthTimeInterval;

		overlap = round((lastValidTime - firstValidTime) / this->azimuthTimeInterval + 1);

		end.at<int>(i - 1, 0) = end.at<int>(i - 1, 0) - int(overlap / 2);

		start.at<int>(i, 0) = this->linesPerBurst * i + this->firstValidLine.at<int>(i, 0) + overlap - int(overlap / 2);

		end.at<int>(i, 0) = this->linesPerBurst * i + this->lastValidLine.at<int>(i, 0);

		lastValidTime = this->burstAzimuthTime.at<double>(i, 0) +
			(this->lastValidLine.at<int>(i, 0) - 1) * this->azimuthTimeInterval;
	}
	end.at<int>(this->burstCount - 1, 0) = this->linesPerBurst * this->burstCount;
	start -= 1;
	//deburst
	ComplexMat tmp, tmp2, slc;
	ret = conversion.read_slc_from_h5(this->h5File.c_str(), tmp);
	tmp.convertTo(tmp, CV_32F);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	slc = tmp(cv::Range(start.at<int>(0, 0), end.at<int>(0, 0)), cv::Range(0, this->samplesPerBurst));
	for (int i = 1; i < this->burstCount; i++)
	{
		tmp2 = tmp(cv::Range(start.at<int>(i, 0), end.at<int>(i, 0)), cv::Range(0, this->samplesPerBurst));
		cv::vconcat(slc.re, tmp2.re, slc.re);
		cv::vconcat(slc.im, tmp2.im, slc.im);
	}
	//往文件里写入卫星数据和辅助数据
	ret = conversion.write_slc_to_h5(outFile, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(this->h5File.c_str(), outFile);
	conversion.write_str_to_h5(outFile, "process_state", "deburst");
	conversion.write_str_to_h5(outFile, "comment", "complex-1.0");
	conversion.write_int_to_h5(outFile, "offset_row", 0);
	conversion.write_int_to_h5(outFile, "offset_col", 0);
	conversion.write_int_to_h5(outFile, "azimuth_len", slc.GetRows());
	conversion.write_int_to_h5(outFile, "range_len", slc.GetCols());

	return 0;
}













DigitalElevationModel::DigitalElevationModel()
{
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	this->lonSpacing = 5.0 / 6000.0;
	this->latSpacing = 5.0 / 6000.0;
}

DigitalElevationModel::~DigitalElevationModel()
{
}

int DigitalElevationModel::getSRTMFileName(
	double lonMin,
	double lonMax,
	double latMin,
	double latMax,
	vector<string>& name
)
{
	if (fabs(lonMin) > 180.0 ||
		fabs(lonMax) > 180.0 ||
		fabs(latMin) >= 60.0 ||
		fabs(latMax) >= 60.0
		)
	{
		fprintf(stderr, "getSRTMFileName(): input check failed!\n");
		return -1;
	}
	name.clear();
	char tmp[512];
	int maxRows = 24; int maxCols = 72; int startRow, endRow, startCol, endCol;
	double spacing = 5.0;
	startRow = (int)((60.0 - latMax) / spacing) + 1;
	endRow = (int)((60.0 - latMin) / spacing) + 1;
	startCol = (int)((lonMin + 180.0) / spacing) + 1;
	endCol = (int)((lonMax + 180.0) / spacing) + 1;
	if (startRow == endRow)
	{
		if (startCol == endCol)
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));
		}
		else
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, startRow);
			}
			name.push_back(string(tmp));
		}
	}
	else
	{
		if (startCol == endCol)
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));

			memset(tmp, 0, 512);
			if (startCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, endRow);
			}
			name.push_back(string(tmp));
		}
		else
		{
			memset(tmp, 0, 512);
			if (startCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, startRow);
			}
			else if (startCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, startRow);
			}
			else if (endCol >= 10 && startRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, startRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, startRow);
			}
			name.push_back(string(tmp));


			memset(tmp, 0, 512);
			if (endCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", endCol, endRow);
			}
			else if (endCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", endCol, endRow);
			}
			else if (endCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", endCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", endCol, endRow);
			}
			name.push_back(string(tmp));

			memset(tmp, 0, 512);
			if (startCol < 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_0%d_0%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow >= 10)
			{
				sprintf(tmp, "srtm_%d_%d.zip", startCol, endRow);
			}
			else if (startCol >= 10 && endRow < 10)
			{
				sprintf(tmp, "srtm_%d_0%d.zip", startCol, endRow);
			}
			else
			{
				sprintf(tmp, "srtm_0%d_%d.zip", startCol, endRow);
			}
			name.push_back(string(tmp));
		}
	}
	
	return 0;
}

int DigitalElevationModel::downloadSRTM(const char* name)
{
	int ret;
	string url = this->SRTMURL + name;
	string savefile = this->DEMPath + string("\\") + name;
	std::replace(savefile.begin(), savefile.end(), '/', '\\');
	HRESULT Result = URLDownloadToFileA(NULL, url.c_str(), savefile.c_str(), 0, NULL);
	if (Result != S_OK)
	{
		fprintf(stderr, "downloadSRTM(): download failded!\n");
		return -1;
	}
	return 0;
}

int DigitalElevationModel::getRawDEM(
	const char* filepath,
	double lonMin,
	double lonMax,
	double latMin,
	double latMax
)
{
	if (!filepath) return -1;
	this->DEMPath = filepath;
	if (GetFileAttributesA(filepath) == -1)
	{
		if (_mkdir(filepath) != 0) return -1;
	}
	vector<string> srtmFileName;
	vector<bool> bAlreadyExist;
	int ret = getSRTMFileName(lonMin, lonMax, latMin, latMax, srtmFileName);
	if (return_check(ret, "getSRTMFileName()", error_head)) return -1;
	//判断文件是否已经存在
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string tmp = this->DEMPath + "\\" + srtmFileName[i];
		std::replace(tmp.begin(), tmp.end(), '/', '\\');
		if (-1 != GetFileAttributesA(tmp.c_str()))bAlreadyExist.push_back(true);
		else bAlreadyExist.push_back(false);
	}
	//不存在则下载
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		if (!bAlreadyExist[i])
		{
			ret = downloadSRTM(srtmFileName[i].c_str());
			if (ret < 0)//未下载到DEM数据,则以0填充
			{
				int rows = (latMax - latMin) / latSpacing;
				int cols = (lonMax - lonMin) / lonSpacing;
				Mat temp = Mat::zeros(rows, cols, CV_16S);
				temp.copyTo(this->rawDEM);
				this->lonUpperLeft = lonMin;
				this->latUpperLeft = latMax;
				return 0;
			}
		}
	}
	//解压文件
	for (int i = 0; i < srtmFileName.size(); i++)
	{
		string folderName = srtmFileName[i];
		folderName = folderName.substr(0, folderName.length() - 4);
		string path = this->DEMPath + string("\\") + folderName;
		std::replace(path.begin(), path.end(), '/', '\\');
		if (-1 != GetFileAttributesA(path.c_str())) continue;
		string srcFile = this->DEMPath + "\\" + srtmFileName[i];
		std::replace(srcFile.begin(), srcFile.end(), '/', '\\');
		if (GetFileAttributesA(srcFile.c_str()) == -1) continue;
		ret = unzip(srcFile.c_str(), path.c_str());
		if (return_check(ret, "unzip()", error_head)) return -1;
	}


	int startRow, startCol, endRow, endCol;
	double lonUpperLeft, lonLowerRight, latUpperLeft, latLowerRight;
	int total_rows, total_cols;

	//DEM在一个SRTM方格内
	if (srtmFileName.size() == 1)
	{
		total_rows = 6000, total_cols = 6000;
		int xx, yy;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		latUpperLeft = 60.0 - (yy - 1) * 5.0;
		latLowerRight = latUpperLeft - 5.0;
		lonUpperLeft = -180.0 + (xx - 1) * 5.0;
		lonLowerRight = lonUpperLeft + 5.0;

		startRow = (latUpperLeft - latMax) / this->latSpacing;
		startRow = startRow < 1 ? 1 : startRow;
		startRow = startRow > total_rows ? total_rows : startRow;
		endRow = (latUpperLeft - latMin) / this->latSpacing;
		endRow = endRow < 1 ? 1 : endRow;
		endRow = endRow > total_rows ? total_rows : endRow;
		startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
		startCol = startCol < 1 ? 1 : startCol;
		startCol = startCol > total_cols ? total_cols : startCol;
		endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
		endCol = endCol < 1 ? 1 : endCol;
		endCol = endCol > total_cols ? total_cols : endCol;

		string folderName = srtmFileName[0];
		folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
		string path = this->DEMPath + string("\\") + folderName;
		path = path + string("\\") + folderName + string(".tif");
		Mat outDEM = Mat::zeros(6000, 6000, CV_16S);
		std::replace(path.begin(), path.end(), '/', '\\');
		ret = geotiffread(path.c_str(), outDEM);
		//if (return_check(ret, "geotiffread()", error_head)) return -1;
		outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
		this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
		this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
	}
	//DEM在2个方格内
	else if (srtmFileName.size() == 2)
	{
		int xx, yy, xx2, yy2;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		sscanf(srtmFileName[1].c_str(), "srtm_%d_%d.zip", &xx2, &yy2);
		//同一列
		if (xx == xx2)
		{
			total_rows = 6000 * 2; total_cols = 6000;
			latUpperLeft = 60.0 - ((yy < yy2 ? yy : yy2) - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = -180.0 + (xx - 1) * 5.0;
			lonLowerRight = lonUpperLeft + 5.0;

			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;


			Mat outDEM, outDEM2;
			
			if (yy < yy2)
			{
				string folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				string path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[1];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM2 = Mat::zeros(6000, 6000, CV_16S);
				ret = geotiffread(path.c_str(), outDEM2);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;
				cv::vconcat(outDEM, outDEM2, outDEM);
			}
			else
			{
				string folderName = srtmFileName[1];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				string path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM = Mat::zeros(6000, 6000, CV_16S);
				ret = geotiffread(path.c_str(), outDEM);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;

				folderName = srtmFileName[0];
				folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
				path = this->DEMPath + string("\\") + folderName;
				path = path + string("\\") + folderName + string(".tif");
				std::replace(path.begin(), path.end(), '/', '\\');
				outDEM2 = Mat::zeros(6000, 6000, CV_16S);
				ret = geotiffread(path.c_str(), outDEM2);
				//if (return_check(ret, "geotiffread()", error_head)) return -1;
				cv::vconcat(outDEM, outDEM2, outDEM);
			}
			
			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
			this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
			this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
		}
		//同一行
		else if(yy == yy2)
		{
			total_cols = 6000 * 2; total_rows = 6000;
			//跨越-180.0/180.0线
			if ((xx == 1 && xx2 == 72) || (xx == 72 && xx2 == 1))
			{
				latUpperLeft = 60.0 - (yy - 1) * 5.0;
				latLowerRight = latUpperLeft - 5.0;
				lonUpperLeft = 175.0;
				lonLowerRight = -175.0;
				startRow = (latUpperLeft - latMax) / this->latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / this->latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMax - lonUpperLeft) / this->lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMin - lonUpperLeft + 360.0) / this->lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;

				Mat outDEM, outDEM2;

				if (xx > xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}
				else
				{
					string folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}

				outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
				this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
				this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
			}
			else
			{
				latUpperLeft = 60.0 - (yy - 1) * 5.0;
				latLowerRight = latUpperLeft - 5.0;
				lonUpperLeft = -180.0 + ((xx < xx2 ? xx : xx2) - 1) * 5.0;
				lonLowerRight = lonUpperLeft + 10.0;

				startRow = (latUpperLeft - latMax) / this->latSpacing;
				startRow = startRow < 1 ? 1 : startRow;
				startRow = startRow > total_rows ? total_rows : startRow;
				endRow = (latUpperLeft - latMin) / this->latSpacing;
				endRow = endRow < 1 ? 1 : endRow;
				endRow = endRow > total_rows ? total_rows : endRow;
				startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
				startCol = startCol < 1 ? 1 : startCol;
				startCol = startCol > total_cols ? total_cols : startCol;
				endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
				endCol = endCol < 1 ? 1 : endCol;
				endCol = endCol > total_cols ? total_cols : endCol;


				Mat outDEM, outDEM2;

				if (xx < xx2)
				{
					string folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}
				else
				{
					string folderName = srtmFileName[1];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					string path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;

					folderName = srtmFileName[0];
					folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
					path = this->DEMPath + string("\\") + folderName;
					path = path + string("\\") + folderName + string(".tif");
					std::replace(path.begin(), path.end(), '/', '\\');
					outDEM2 = Mat::zeros(6000, 6000, CV_16S);
					ret = geotiffread(path.c_str(), outDEM2);
					//if (return_check(ret, "geotiffread()", error_head)) return -1;
					cv::hconcat(outDEM, outDEM2, outDEM);
				}

				outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
				this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
				this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
			}			
		}
		else
		{
			return -1;
		}
		

		
	}
	//DEM在4个方格内
	else if (srtmFileName.size() == 4)
	{
		int xx, yy, xx2, yy2, xx3, yy3, xx4, yy4, temp;
		sscanf(srtmFileName[0].c_str(), "srtm_%d_%d.zip", &xx, &yy);
		sscanf(srtmFileName[1].c_str(), "srtm_%d_%d.zip", &xx2, &yy2);
		sscanf(srtmFileName[2].c_str(), "srtm_%d_%d.zip", &xx3, &yy3);
		sscanf(srtmFileName[3].c_str(), "srtm_%d_%d.zip", &xx4, &yy4);
		total_rows = 6000 * 2; total_cols = 6000 * 2;
		//跨越-180.0/180.0线
		if (lonMax * lonMin < 0 && (fabs(lonMin) + fabs(lonMax)) > 180.0)
		{
			startRow = (int)((60.0 - latMax) / 5.0) + 1;
			endRow = (int)((60.0 - latMin) / 5.0) + 1;
			endCol = (int)((lonMin + 180.0) / 5.0) + 1;
			startCol = (int)((lonMax + 180.0) / 5.0) + 1;
			latUpperLeft = 60.0 - (startRow - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = 175.0;
			lonLowerRight = -175.0;

			

			Mat outDEM, outDEM2, outDEM3;

			char tmpstr[512];
			const char* format = NULL;
			if (startCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, startRow);
			string folderName(tmpstr);
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			string path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, startRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM, outDEM2, outDEM);

			if (startCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);


			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMin - lonUpperLeft + 360.0) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
			this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
			this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
		}
		else
		{
			startRow = (int)((60.0 - latMax) / 5.0) + 1;
			endRow = (int)((60.0 - latMin) / 5.0) + 1;
			startCol = (int)((lonMin + 180.0) / 5.0) + 1;
			endCol = (int)((lonMax + 180.0) / 5.0) + 1;
			latUpperLeft = 60.0 - (startRow - 1) * 5.0;
			latLowerRight = latUpperLeft - 10.0;
			lonUpperLeft = -180.0 + (startCol - 1) * 5.0;
			lonLowerRight = lonUpperLeft + 10.0;

			

			Mat outDEM, outDEM2, outDEM3;

			char tmpstr[512];
			const char* format = NULL;
			if (startCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, startRow);
			string folderName(tmpstr);
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			string path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && startRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && startRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && startRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, startRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM, outDEM2, outDEM);

			if (startCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (startCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (startCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, startCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM2 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM2);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;

			if (endCol < 10 && endRow < 10) format = "srtm_0%d_0%d.zip";
			else if (endCol >= 10 && endRow < 10) format = "srtm_%d_0%d.zip";
			else if (endCol < 10 && endRow >= 10) format = "srtm_0%d_%d.zip";
			else format = "srtm_%d_%d.zip";
			sprintf(tmpstr, format, endCol, endRow);
			folderName = tmpstr;
			folderName = folderName.substr(0, folderName.length() - 4);//去掉.zip后缀
			path = this->DEMPath + string("\\") + folderName;
			path = path + string("\\") + folderName + string(".tif");
			std::replace(path.begin(), path.end(), '/', '\\');
			outDEM3 = Mat::zeros(6000, 6000, CV_16S);
			ret = geotiffread(path.c_str(), outDEM3);
			//if (return_check(ret, "geotiffread()", error_head)) return -1;
			cv::hconcat(outDEM2, outDEM3, outDEM2);

			cv::vconcat(outDEM, outDEM2, outDEM);

			startRow = (latUpperLeft - latMax) / this->latSpacing;
			startRow = startRow < 1 ? 1 : startRow;
			startRow = startRow > total_rows ? total_rows : startRow;
			endRow = (latUpperLeft - latMin) / this->latSpacing;
			endRow = endRow < 1 ? 1 : endRow;
			endRow = endRow > total_rows ? total_rows : endRow;
			startCol = (lonMin - lonUpperLeft) / this->lonSpacing;
			startCol = startCol < 1 ? 1 : startCol;
			startCol = startCol > total_cols ? total_cols : startCol;
			endCol = (lonMax - lonUpperLeft) / this->lonSpacing;
			endCol = endCol < 1 ? 1 : endCol;
			endCol = endCol > total_cols ? total_cols : endCol;

			outDEM(cv::Range(startRow - 1, endRow), cv::Range(startCol - 1, endCol)).copyTo(this->rawDEM);
			this->lonUpperLeft = lonUpperLeft + (startCol - 1) * lonSpacing;
			this->latUpperLeft = latUpperLeft - (startRow - 1) * latSpacing;
		}

	}
	else return -1;
	this->rows = this->rawDEM.rows;
	this->cols = this->rawDEM.cols;
	return 0;
}

int DigitalElevationModel::getElevation(double lon, double lat, double* elevation)
{
	if (!elevation) return -1;
	int row, col;
	row = (this->latUpperLeft - lat) / this->latSpacing;
	col = (lon - this->lonUpperLeft) / this->lonSpacing;
	if (row < 0 || col < 0 || row >= this->rows || col >= this->cols) return -1;
	double elevationUL, elevationUR, elevationLL, elevationLR, upper, lower;
	int r, c, r1, c1;
	r = row; c = col; r1 = r + 1; c1 = c + 1;
	r1 = r1 > this->rows - 1 ? this->rows - 1 : r1;
	c1 = c1 > this->cols - 1 ? this->cols - 1 : c1;
	elevationUL = this->rawDEM.at<short>(r, c);
	elevationUR = this->rawDEM.at<short>(r, c1);
	elevationLL = this->rawDEM.at<short>(r1, c);
	elevationLR = this->rawDEM.at<short>(r1, c1);
	*elevation = (elevationLL + elevationLR + elevationUL + elevationUR) / 4.0;
	return 0;
}

int DigitalElevationModel::geotiffread(const char* filename, Mat& outDEM)
{
	if (!filename) return -1;
	GDALAllRegister();	//注册已知驱动
	GDALDataset* poDataset = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);	//打开geotiff文件
	if (poDataset == NULL)
	{
		fprintf(stderr, "geotiffread(): failed to open %s!\n", filename);
		GDALDestroyDriverManager();
		return -1;
	}
	int nBand = poDataset->GetRasterCount();	//获取波段数（geotiff应为1）
	int xsize = 0;
	int ysize = 0;
	if (nBand == 1)
	{
		GDALRasterBand* poBand = poDataset->GetRasterBand(1);	//获取指向波段1的指针
		xsize = poBand->GetXSize();		//cols
		ysize = poBand->GetYSize();		//rows
		if (xsize < 0 || ysize < 0)
		{
			fprintf(stderr, "geotiffread(): band rows and cols error!\n");
			GDALClose(poDataset);
			GDALDestroyDriverManager();
			return -1;
		}
		GDALDataType dataType = poBand->GetRasterDataType();	//数据存储类型，geotiff应为16位整型
		short* pbuf = NULL;
		pbuf = (short*)malloc(sizeof(short) * xsize * ysize);		//分配数据指针空间
		if (!pbuf)
		{
			fprintf(stderr, "geotiffread(): out of memory!\n");
			GDALClose(poDataset);
			GDALDestroyDriverManager();
			return -1;
		}
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize, pbuf, xsize, ysize, dataType, 0, 0);		//读取复图像数据到pbuf中
		int i, j;
		outDEM.create(ysize, xsize, CV_16S);
		memcpy(outDEM.data, pbuf, sizeof(short) * xsize * ysize);
		if (pbuf)
		{
			free(pbuf);
			pbuf = NULL;
		}
		GDALClose(poDataset);
		GDALDestroyDriverManager();
	}
	else
	{
		fprintf(stderr, "geotiffread(): number of Bands != 1\n");
		GDALClose(poDataset);
		GDALDestroyDriverManager();
		return -1;
	}
	int rows = outDEM.rows;
	int cols = outDEM.cols;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (outDEM.at<short>(i, j) < 0) outDEM.at<short>(i, j) = 0;
		}
	}
	return 0;
}

int DigitalElevationModel::unzip(const char* srcFile, const char* dstPath)
{
	if (!srcFile || !dstPath)
	{
		fprintf(stderr, "unzip(): input check failed!\n");
		return -1;
	}
	int ret;
	//如果目标文件夹不存在，则创建
	if (-1 == GetFileAttributesA(dstPath))
	{
		ret = _mkdir(dstPath);
		if (ret < 0)return -1;
	}
	//////////////////////////创建并调用unzip.exe进程///////////////////////////////
	char szFilePath[MAX_PATH + 1] = { 0 };
	GetModuleFileNameA(NULL, szFilePath, MAX_PATH);
	string str(szFilePath);
	str = str.substr(0, str.rfind("\\"));
	string commandline = str + string("\\unzip.exe ") + string(srcFile) + string(" ") + string(dstPath);
	char szCommandLine[1024];
	strcpy(szCommandLine, commandline.c_str());
	STARTUPINFOA si;
	PROCESS_INFORMATION p_i;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&p_i, sizeof(p_i));
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = FALSE;
	BOOL bRet = ::CreateProcessA(
		NULL,           // 不在此指定可执行文件的文件名
		szCommandLine,      // 命令行参数
		NULL,           // 默认进程安全性
		NULL,           // 默认线程安全性
		FALSE,          // 指定当前进程内的句柄不可以被子进程继承
		CREATE_NEW_CONSOLE, // 为新进程创建一个新的控制台窗口
		NULL,           // 使用本进程的环境变量
		NULL,           // 使用本进程的驱动器和目录
		&si,
		&p_i);
	if (bRet)
	{
		char snaphu_job_name[512]; snaphu_job_name[0] = 0;
		time_t tt = std::time(0);
		sprintf(snaphu_job_name, "UNZIP_%lld", tt);
		string snaphu_job_name_string(snaphu_job_name);
		HANDLE hd = CreateJobObjectA(NULL, snaphu_job_name_string.c_str());
		if (hd)
		{
			JOBOBJECT_EXTENDED_LIMIT_INFORMATION extLimitInfo;
			extLimitInfo.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
			BOOL retval = SetInformationJobObject(hd, JobObjectExtendedLimitInformation, &extLimitInfo, sizeof(extLimitInfo));
			if (retval)
			{
				if (p_i.hProcess)
				{
					retval = AssignProcessToJobObject(hd, p_i.hProcess);
				}
			}
		}
		WaitForSingleObject(p_i.hProcess, INFINITE);
		::CloseHandle(p_i.hThread);
		::CloseHandle(p_i.hProcess);
	}
	else
	{
		fprintf(stderr, "unzip(): create unzip.exe process failed!\n\n");
		return -1;
	}
	return 0;
}






















Sentinel1BackGeocoding::Sentinel1BackGeocoding()
{
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "FORMATCONVERSION_DLL_ERROR: error happens when using ");
	this->dem = NULL;
	this->burstOffsetComputed = false;
	this->isMasterRgAzComputed = false;
	this->isdeBurstConfig = false;
	this->masterIndex = 1;
	this->numOfImages = 0;
}

Sentinel1BackGeocoding::~Sentinel1BackGeocoding()
{
	if (dem)
	{
		delete dem; dem = NULL;
	}
	for (int i = 0; i < su.size(); i++)
	{
		if (su[i])
		{
			delete su[i]; su[i] = NULL;
		}
	}
}

int Sentinel1BackGeocoding::init(
	vector<string>& h5Files,
	vector<string>& outFiles,
	const char* DEMPath,
	int masterIndex
)
{
	int ret;
	ret = loadData(h5Files);
	if (return_check(ret, "loadData()", error_head)) return -1;
	ret = setDEMPath(DEMPath);
	if (return_check(ret, "setDEMPath()", error_head)) return -1;
	ret = loadOutFiles(outFiles);
	if (return_check(ret, "loadOutFiles()", error_head)) return -1;
	ret = setMasterIndex(masterIndex);
	if (return_check(ret, "setMasterIndex()", error_head)) return -1;
	deBurstConfig();
	ret = prepareOutFiles();
	if (return_check(ret, "prepareOutFiles()", error_head)) return -1;

	

	return 0;
}

int Sentinel1BackGeocoding::loadData(vector<string>& h5Files)
{
	if (h5Files.size() < 2)
	{
		fprintf(stderr, "loadData(): input check failed!\n");
		return -1;
	}
	int ret;
	this->numOfImages = h5Files.size();
	//清空已有数据
	for (int i = 0; i < su.size(); i++)
	{
		if (su[i])
		{
			delete su[i]; su[i] = NULL;
		}
	}
	su.clear();
	for (int i = 0; i < numOfImages; i++)
	{
		su.push_back(new Sentinel1Utils(h5Files[i].c_str()));
		if (su[i])
		{
			ret = su[i]->init();
			if (return_check(ret, "init()", error_head)) return -1;
		}
	}
	return 0;
}

int Sentinel1BackGeocoding::setDEMPath(const char* DEMPath)
{
	if (!DEMPath)
	{
		fprintf(stderr, "setDEMPath(): input check failed!\n");
		return -1;
	}
	this->DEMPath = DEMPath;
	return 0;
}

int Sentinel1BackGeocoding::loadDEM(
	const char* filepath,
	double lonMin, 
	double lonMax,
	double latMin,
	double latMax
)
{
	int ret;
	if (dem) {
		delete dem; dem = NULL;
	}
	dem = new DigitalElevationModel();
	ret = dem->getRawDEM(filepath, lonMin, lonMax, latMin, latMax);
	if (return_check(ret, "getRawDEM()", error_head)) return -1;
	return 0;
}

int Sentinel1BackGeocoding::loadOutFiles(vector<string>& outFiles)
{
	if (outFiles.size() != su.size())
	{
		fprintf(stderr, "loadOutFiles(): input check failed!\n");
		return -1;
	}
	this->outFiles.clear();
	for (int i = 0; i < outFiles.size(); i++)
	{
		this->outFiles.push_back(outFiles[i]);
	}
	return 0;
}

int Sentinel1BackGeocoding::prepareOutFiles()
{
	if (!isdeBurstConfig) deBurstConfig();
	int ret; FormatConversion conversion;
	if (numOfImages < 2) return -1;
	ComplexMat tmp, tmp2, slc;
	ret = conversion.read_slc_from_h5(su[masterIndex - 1]->h5File.c_str(), tmp);
	tmp.convertTo(tmp, CV_32F);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	//deburst
	slc = tmp(cv::Range(start.at<int>(0, 0), end.at<int>(0, 0)), cv::Range(0, su[masterIndex - 1]->samplesPerBurst));
	for (int i = 1; i < su[masterIndex - 1]->burstCount; i++)
	{
		tmp2 = tmp(cv::Range(start.at<int>(i, 0), end.at<int>(i, 0)), cv::Range(0, su[masterIndex - 1]->samplesPerBurst));
		cv::vconcat(slc.re, tmp2.re, slc.re);
		cv::vconcat(slc.im, tmp2.im, slc.im);
	}
	for (int i = 0; i < numOfImages; i++)
	{
		conversion.creat_new_h5(this->outFiles[i].c_str());
		ret = conversion.write_slc_to_h5(this->outFiles[i].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	}
	this->deburstLines = slc.GetRows();
	return 0;
}

int Sentinel1BackGeocoding::setMasterIndex(int masterIndex)
{
	if (masterIndex < 1 || masterIndex > numOfImages)
	{
		fprintf(stderr, "setMasterIndex(): input check failed!\n");
		return -1;
	}
	this->masterIndex = masterIndex;
	return 0;
}

int Sentinel1BackGeocoding::computeBurstOffset()
{
	if (burstOffsetComputed) return 0;
	if (!dem || su.size() < 2)
	{
		fprintf(stderr, "computeBurstOffset(): input check failed!\n");
		return -1;
	}
	int ret, numOfGeoLocationPoints;
	Position earthPoint;
	double lon, lat, elevation;
	numOfGeoLocationPoints = su[masterIndex - 1]->geolocationGridPoint.rows;
	for (int i = 0; i < numOfGeoLocationPoints; i++)
	{
		lon = su[masterIndex - 1]->geolocationGridPoint.at<double>(i, 0);
		lat = su[masterIndex - 1]->geolocationGridPoint.at<double>(i, 1);
		dem->getElevation(lon, lat, &elevation);
		Utils::ell2xyz(lon, lat, elevation, earthPoint);
		BurstIndices mBurstIndices, sBurstIndices;
		ret = su[masterIndex - 1]->getBurstIndice(earthPoint, mBurstIndices);
		if (ret < 0) continue;
		for (int j = 0; j < numOfImages; j++)
		{
			if (j == masterIndex - 1) continue;
			ret = su[j]->getBurstIndice(earthPoint, sBurstIndices);
			if (ret < 0 || (mBurstIndices.firstBurstIndex == -1 && mBurstIndices.secondBurstIndex == -1) ||
				(sBurstIndices.firstBurstIndex == -1 && sBurstIndices.secondBurstIndex == -1)) {
				continue;
			}
			if (mBurstIndices.inUpperPartOfFirstBurst == sBurstIndices.inUpperPartOfFirstBurst) {
				su[j]->burstOffset = sBurstIndices.firstBurstIndex - mBurstIndices.firstBurstIndex;
			}
			else if (sBurstIndices.secondBurstIndex != -1 &&
				mBurstIndices.inUpperPartOfFirstBurst == sBurstIndices.inUpperPartOfSecondBurst) {
				su[j]->burstOffset = sBurstIndices.secondBurstIndex - mBurstIndices.firstBurstIndex;
			}
			else if (mBurstIndices.secondBurstIndex != -1 &&
				mBurstIndices.inUpperPartOfSecondBurst == sBurstIndices.inUpperPartOfFirstBurst) {
				su[j]->burstOffset = sBurstIndices.firstBurstIndex - mBurstIndices.secondBurstIndex;
			}
			else if (mBurstIndices.secondBurstIndex != -1 && sBurstIndices.secondBurstIndex != -1 &&
				mBurstIndices.inUpperPartOfSecondBurst == sBurstIndices.inUpperPartOfSecondBurst) {
				su[j]->burstOffset = sBurstIndices.secondBurstIndex - mBurstIndices.secondBurstIndex;
			}
		}
		bool allComputed = true;
		for (int j = 0; j < numOfImages; j++) {
			if (j == masterIndex - 1) continue;
			if (su[j]->burstOffset == -9999) {
				allComputed = false;
				break;
			}
		}
		if (!allComputed)
			continue;

		burstOffsetComputed = true;
		return 0;
	}
	for (int j = 0; j < numOfImages; j++) {
		if (j == masterIndex - 1) continue;
		su[j]->burstOffset = 0;
	}
	burstOffsetComputed = true;
	return 0;
}

int Sentinel1BackGeocoding::performDerampDemod(Mat& derampDemodPhase, ComplexMat& slc)
{
	if (derampDemodPhase.size != slc.re.size || derampDemodPhase.size != slc.im.size || derampDemodPhase.empty())
	{
		fprintf(stderr, "performDerampDemod(): input check failed!\n");
		return -1;
	}
	ComplexMat tmp;
	Utils util;
	if (derampDemodPhase.type() != CV_64F) derampDemodPhase.convertTo(derampDemodPhase, CV_64F);
	util.phase2cos(derampDemodPhase, tmp.re, tmp.im);
	if (slc.type() != CV_64F) slc.convertTo(slc, CV_64F);
	slc.Mul(tmp, slc, false);

	return 0;
}

int Sentinel1BackGeocoding::computeSlavePosition(int slaveImagesIndex, int mBurstIndex)
{
	if (slaveImagesIndex < 1 || slaveImagesIndex > numOfImages)
	{
		fprintf(stderr, "computeSlavePosition(): input check failed!\n");
		return -1;
	}
	int sBurstIndex = mBurstIndex + su[slaveImagesIndex - 1]->burstOffset;
	if (sBurstIndex < 1 || sBurstIndex > su[slaveImagesIndex - 1]->burstCount) {
		return -1;
	}
	int ret;
	double lonMin, lonMax, latMin, latMax;
	if (!isMasterRgAzComputed)
	{
		masterAzimuth.create(dem->rows, dem->cols, CV_64F);
		slaveAzimuth.create(dem->rows, dem->cols, CV_64F);
		masterRange.create(dem->rows, dem->cols, CV_64F);
		slaveRange.create(dem->rows, dem->cols, CV_64F);
	}
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < dem->rows; i++)
	{
		double lon, lat, elevation, rangeIndex, azimuthIndex;
		Position earthPoint;
		for (int j = 0; j < dem->cols; j++)
		{
			lat = dem->latUpperLeft - i * dem->latSpacing;
			lon = dem->lonUpperLeft + j * dem->lonSpacing;
			lon = lon > 180.0 ? lon - 360.0 : lon;
			elevation = dem->rawDEM.at<short>(i, j);
			Utils::ell2xyz(lon, lat, elevation, earthPoint);
			if (!isMasterRgAzComputed)
			{
				if (su[masterIndex - 1]->getRgAzPosition(mBurstIndex, earthPoint, &rangeIndex, &azimuthIndex) == 0)
				{
					masterAzimuth.at<double>(i, j) = azimuthIndex;
					masterRange.at<double>(i, j) = rangeIndex;
				}
				else
				{
					masterAzimuth.at<double>(i, j) = invalidRgAzIndex;
					masterRange.at<double>(i, j) = invalidRgAzIndex;
				}
			}
			if (su[slaveImagesIndex - 1]->getRgAzPosition(sBurstIndex, earthPoint, &rangeIndex, &azimuthIndex) == 0)
			{
				slaveAzimuth.at<double>(i, j) = azimuthIndex;
				slaveRange.at<double>(i, j) = rangeIndex;
			}
			else
			{
				slaveAzimuth.at<double>(i, j) = invalidRgAzIndex;
				slaveRange.at<double>(i, j) = invalidRgAzIndex;
			}
		}
	}
	
	if (!isMasterRgAzComputed)isMasterRgAzComputed = true;
	return 0;
}

int Sentinel1BackGeocoding::computeSlaveOffset(Mat& slaveAzimuthOffset, Mat& slaveRangeOffset)
{
	if (!isMasterRgAzComputed)
	{
		fprintf(stderr, "computeSlaveOffset(): input check failed!\n");
		return -1;
	}
	slaveAzimuthOffset.create(masterRange.size(), CV_64F);
	slaveRangeOffset.create(masterRange.size(), CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < slaveAzimuthOffset.rows; i++)
	{
		for (int j = 0; j < slaveAzimuthOffset.cols; j++)
		{
			//方位向偏移量计算
			if (masterAzimuth.at<double>(i, j) < -0.5 || slaveAzimuth.at<double>(i, j) < -0.5)
			{
				slaveAzimuthOffset.at<double>(i, j) = invalidOffset;
			}
			else
			{
				slaveAzimuthOffset.at<double>(i, j) = slaveAzimuth.at<double>(i, j) - masterAzimuth.at<double>(i, j);
			}
			//距离向偏移量计算
			if (masterRange.at<double>(i, j) < -0.5 || slaveRange.at<double>(i, j) < -0.5)
			{
				slaveRangeOffset.at<double>(i, j) = invalidOffset;
			}
			else
			{
				slaveRangeOffset.at<double>(i, j) = slaveRange.at<double>(i, j) - masterRange.at<double>(i, j);
			}
		}
	}
	return 0;
}

int Sentinel1BackGeocoding::fitSlaveOffset(
	Mat& slaveOffset,
	double* a0,
	double* a1, 
	double* a2
)
{
	if (slaveOffset.empty() || slaveOffset.type() != CV_64F)
	{
		fprintf(stderr, "fitSlaveOffset(): input check failed!\n");
		return -1;
	}
	int count = 0, nr, nc;
	nr = slaveOffset.rows;
	nc = slaveOffset.cols;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(slaveOffset.at<double>(i, j) - invalidOffset) > 0.0001) count++;
		}
	}
	if (count < 4)
	{
		fprintf(stderr, "fitSlaveOffset(): not enough valid offsetpoints!\n");
		return -1;
	}
	Mat offset(count, 1, CV_64F);
	Mat range(count, 1, CV_64F);
	Mat azimuth(count, 1, CV_64F);
	Mat A = Mat::ones(count, 3, CV_64F);
	count = 0;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(slaveOffset.at<double>(i, j) - invalidOffset) > 0.0001)
			{
				offset.at<double>(count, 0) = slaveOffset.at<double>(i, j);
				range.at<double>(count, 0) = masterRange.at<double>(i, j);
				azimuth.at<double>(count++, 0) = masterAzimuth.at<double>(i, j);
			}
		}
	}
	range.copyTo(A(cv::Range(0, count), cv::Range(1, 2)));
	azimuth.copyTo(A(cv::Range(0, count), cv::Range(2, 3)));
	Mat A_t, b, coef;
	cv::transpose(A, A_t);
	A = A_t * A;
	b = A_t * offset;
	if (!cv::solve(A, b, coef, cv::DECOMP_NORMAL))
	{
		fprintf(stderr, "fitSlaveOffset(): matrix defficiency!\n");
		return -1;
	}
	if (a0) *a0 = coef.at<double>(0, 0);
	if (a1) *a1 = coef.at<double>(1, 0);
	if (a2) *a2 = coef.at<double>(2, 0);
	return 0;
}

int Sentinel1BackGeocoding::performBilinearResampling(
	ComplexMat& slave,
	int dstHeight, 
	int dstWidth, 
	double a0Rg, double a1Rg, double a2Rg,
	double a0Az, double a1Az, double a2Az
)
{
	if (slave.isempty() || dstHeight < 2 || dstWidth < 2)
	{
		fprintf(stderr, "performBilinearResampling(): input check failed!\n");
		return -1;
	}
	ComplexMat slcResampled;
	if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
	slcResampled.re.create(dstHeight, dstWidth, CV_64F);
	slcResampled.im.create(dstHeight, dstWidth, CV_64F);
	Mat coef_r(3, 1, CV_64F), coef_c(3, 1, CV_64F);
	coef_r.at<double>(0, 0) = a0Az;
	coef_r.at<double>(1, 0) = a1Az;
	coef_r.at<double>(2, 0) = a2Az;
	coef_c.at<double>(0, 0) = a0Rg;
	coef_c.at<double>(1, 0) = a1Rg;
	coef_c.at<double>(2, 0) = a2Rg;
	int rows = dstHeight; int cols = dstWidth;
	int cols_slave = slave.GetCols(); int rows_slave = slave.GetRows();
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < rows; i++)
	{
		double x, y, ii, jj; Mat tmp(1, 3, CV_64F); Mat result;
		int mm, nn, mm1, nn1;
		double offset_rows, offset_cols, upper, lower;
		for (int j = 0; j < cols; j++)
		{
			jj = (double)j;
			ii = (double)i;
			tmp.at<double>(0, 0) = 1.0;
			tmp.at<double>(0, 1) = jj;
			tmp.at<double>(0, 2) = ii;

			result = tmp * coef_r;
			offset_rows = result.at<double>(0, 0);
			result = tmp * coef_c;
			offset_cols = result.at<double>(0, 0);

			ii += offset_rows;
			jj += offset_cols;

			mm = (int)floor(ii); nn = (int)floor(jj);
			if (mm < 0 || nn < 0 || mm > rows_slave - 1 || nn > cols_slave - 1)
			{
				slcResampled.re.at<double>(i, j) = 0.0;
				slcResampled.im.at<double>(i, j) = 0.0;
			}
			else
			{
				mm1 = mm + 1; nn1 = nn + 1;
				mm1 = mm1 >= rows_slave - 1 ? rows_slave - 1 : mm1;
				nn1 = nn1 >= cols_slave - 1 ? cols_slave - 1 : nn1;
				//实部插值
				upper = slave.re.at<double>(mm, nn) + (slave.re.at<double>(mm, nn1) - slave.re.at<double>(mm, nn)) * (jj - (double)nn);
				lower = slave.re.at<double>(mm1, nn) + (slave.re.at<double>(mm1, nn1) - slave.re.at<double>(mm1, nn)) * (jj - (double)nn);
				slcResampled.re.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
				//虚部插值
				upper = slave.im.at<double>(mm, nn) + (slave.im.at<double>(mm, nn1) - slave.im.at<double>(mm, nn)) * (jj - (double)nn);
				lower = slave.im.at<double>(mm1, nn) + (slave.im.at<double>(mm1, nn1) - slave.im.at<double>(mm1, nn)) * (jj - (double)nn);
				slcResampled.im.at<double>(i, j) = upper + (lower - upper) * (ii - (double)mm);
			}

		}
	}
	slave = slcResampled;
	return 0;
}

int Sentinel1BackGeocoding::slaveBilinearInterpolation(
	int mBurstIndex,
	int slaveImageIndex,
	ComplexMat& slave
)
{
	if (slaveImageIndex < 1 || 
		slaveImageIndex > numOfImages ||
		mBurstIndex < 1 || 
		mBurstIndex > su[masterIndex - 1]->burstCount)
	{
		fprintf(stderr, "slaveBilinearInterpolation(): input check failed!\n");
		return -1;
	}
	int ret;
	int sBurstIndex = mBurstIndex + su[slaveImageIndex - 1]->burstOffset;
	if (sBurstIndex < 1 || sBurstIndex > su[slaveImageIndex - 1]->burstCount) {
		return -1;
	}
	ComplexMat tmp;
	double a0Rg, a1Rg, a2Rg, a0Az, a1Az, a2Az;
	ret = su[slaveImageIndex - 1]->getBurst(sBurstIndex, slave);
	if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
	if (return_check(ret, "getBurst()", error_head)) return -1;
	Mat derampDemodPhase;
	ret = su[slaveImageIndex - 1]->computeDerampDemodPhase(sBurstIndex, derampDemodPhase);
	if (return_check(ret, "computeDerampDemodPhase()", error_head)) return -1;
	ret = performDerampDemod(derampDemodPhase, slave);
	if (return_check(ret, "performDerampDemod()", error_head)) return -1;
	ret = computeSlavePosition(slaveImageIndex, mBurstIndex);
	if (return_check(ret, "computeSlavePosition()", error_head)) return -1;
	Mat slaveAzimuthOffset, slaveRangeOffset;
	ret = computeSlaveOffset(slaveAzimuthOffset, slaveRangeOffset);
	if (return_check(ret, "computeSlaveOffset()", error_head)) return -1;
	Utils util;
	ret = fitSlaveOffset(slaveAzimuthOffset, &a0Az, &a1Az, &a2Az);
	if (return_check(ret, "fitSlaveOffset()", error_head)) return -1;
	ret = fitSlaveOffset(slaveRangeOffset, &a0Rg, &a1Rg, &a2Rg);
	if (return_check(ret, "fitSlaveOffset()", error_head)) return -1;
	ret = performBilinearResampling(slave, su[masterIndex - 1]->linesPerBurst, su[masterIndex - 1]->samplesPerBurst,
		a0Rg, a1Rg, a2Rg, a0Az, a1Az, a2Az);
	if (return_check(ret, "performBilinearResampling()", error_head)) return -1;
	tmp.SetRe(derampDemodPhase); tmp.SetIm(derampDemodPhase);
	ret = performBilinearResampling(tmp, su[masterIndex - 1]->linesPerBurst, su[masterIndex - 1]->samplesPerBurst,
		a0Rg, a1Rg, a2Rg, a0Az, a1Az, a2Az);
	if (return_check(ret, "performBilinearResampling()", error_head)) return -1;
	tmp.re.copyTo(derampDemodPhase);
	util.phase2cos(derampDemodPhase, tmp.re, tmp.im);
	slave.Mul(tmp, slave, true);//reramp
	slave.convertTo(slave, CV_32F);
	return 0;
}

int Sentinel1BackGeocoding::deBurstConfig()
{
	if (isdeBurstConfig) return 0;
	Mat start(su[masterIndex - 1]->burstCount, 1, CV_32S), end(su[masterIndex - 1]->burstCount, 1, CV_32S);
	start.at<int>(0, 0) = 1;
	end.at<int>(0, 0) = su[masterIndex - 1]->lastValidLine.at<int>(0, 0);
	double lastValidTime = su[masterIndex - 1]->burstAzimuthTime.at<double>(0, 0) +
		(su[masterIndex - 1]->lastValidLine.at<int>(0, 0) - 1) * su[masterIndex - 1]->azimuthTimeInterval;
	double firstValidTime;
	int deburstLines = 0;
	int overlap;
	for (int i = 1; i < su[masterIndex - 1]->burstCount; i++)
	{
		firstValidTime = su[masterIndex - 1]->burstAzimuthTime.at<double>(i, 0) + (su[masterIndex - 1]->firstValidLine.at<int>(i, 0) - 1) *
			su[masterIndex - 1]->azimuthTimeInterval;

		overlap = round((lastValidTime - firstValidTime) / su[masterIndex - 1]->azimuthTimeInterval + 1);

		end.at<int>(i - 1, 0) = end.at<int>(i - 1, 0) - int(overlap / 2);

		start.at<int>(i, 0) = su[masterIndex - 1]->linesPerBurst * i + su[masterIndex - 1]->firstValidLine.at<int>(i, 0) + overlap - int(overlap / 2);

		end.at<int>(i, 0) = su[masterIndex - 1]->linesPerBurst * i + su[masterIndex - 1]->lastValidLine.at<int>(i, 0);

		lastValidTime = su[masterIndex - 1]->burstAzimuthTime.at<double>(i, 0) +
			(su[masterIndex - 1]->lastValidLine.at<int>(i, 0) - 1) * su[masterIndex - 1]->azimuthTimeInterval;
	}
	end.at<int>(su[masterIndex - 1]->burstCount - 1, 0) = su[masterIndex - 1]->linesPerBurst * su[masterIndex - 1]->burstCount;
	start -= 1;
	//end -= 1;
	start.copyTo(this->start);
	end.copyTo(this->end);
	//for (int i = 0; i < su[masterIndex - 1]->burstCount; i++)
	//{
	//	deburstLines += end.at<int>(i, 0) - start.at<int>(i, 0);
	//}
	//this->deburstLines = deburstLines;
	isdeBurstConfig = true;
	return 0;
}

int Sentinel1BackGeocoding::backGeoCodingCoregistration()
{
	if (!isdeBurstConfig) deBurstConfig();
	FormatConversion conversion;
	ComplexMat slaveSLC, tmp;
	Utils util;
	int linesPerBurst, ret;
	int samplesPerBurst = su[masterIndex - 1]->samplesPerBurst;
	int offset_row = 0, lines = 0;
	double lonMin, lonMax, latMin, latMax;
	ret = su[masterIndex - 1]->computeImageGeoBoundry(&lonMin, &lonMax, &latMin, &latMax);
	if (return_check(ret, "computeImageGeoBoundry()", error_head)) return -1;
	ret = loadDEM(this->DEMPath.c_str(), lonMin, lonMax, latMin, latMax);
	if (return_check(ret, "loadDEM()", error_head)) return -1;
	for (int i = 0; i < su[masterIndex - 1]->burstCount; i++)
	{
		if (!burstOffsetComputed)
		{
			int ret = computeBurstOffset();
			if (return_check(ret, "computeBurstOffset()", error_head)) return -1;
		}
		lines = i * su[masterIndex - 1]->linesPerBurst;
		linesPerBurst = end.at<int>(i, 0) - start.at<int>(i, 0);
		for (int j = 0; j < numOfImages; j++)
		{
			if (j == masterIndex - 1) continue;
			ret = slaveBilinearInterpolation(i + 1, j + 1, slaveSLC);
			if (return_check(ret, "slaveBilinearInterpolation()", error_head)) return -1;
			tmp = slaveSLC(cv::Range(start.at<int>(i, 0) - lines, end.at<int>(i, 0) - lines), 
				cv::Range(0, su[masterIndex - 1]->samplesPerBurst));
			ret = conversion.write_subarray_to_h5(this->outFiles[j].c_str(), "s_re", tmp.re,
				offset_row, 0, linesPerBurst, samplesPerBurst);
			if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
			ret = conversion.write_subarray_to_h5(this->outFiles[j].c_str(), "s_im", tmp.im,
				offset_row, 0, linesPerBurst, samplesPerBurst);
			if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
		}
		offset_row += linesPerBurst;
		isMasterRgAzComputed = false;
	}
	return 0;
}


























orbitStateVectors::orbitStateVectors(Mat& stateVectors, double startTime, double stopTime)
{
	this->startTime = startTime;
	this->stopTime = stopTime;
	this->isOrbitUpdated = false;
	stateVectors.copyTo(this->stateVectors);
	if (this->stateVectors.type() != CV_64F) this->stateVectors.convertTo(this->stateVectors, CV_64F);
	this->dt = 0.00000001;
	setSceneStartStopTime(startTime, stopTime);

}

orbitStateVectors::~orbitStateVectors()
{
}

int orbitStateVectors::setSceneStartStopTime(double startTime, double stopTime)
{
	this->startTime = startTime;
	this->stopTime = stopTime;
	return 0;
}

int orbitStateVectors::getPosition(double azimuthTime, Position& position)
{
	if (newStateVectors.cols != 7 || newStateVectors.rows < 2 || !isOrbitUpdated)
	{
		fprintf(stderr, "getPosition(): input check failed!\n");
		return -1;
	}
	if (azimuthTime < newStateVectors.at<double>(0, 0) || azimuthTime > newStateVectors.at<double>(newStateVectors.rows - 1, 0))
	{
		fprintf(stderr, "getPosition(): azimuthTime out of legal range\n");
		return -1;
	}
	int i0, iN;
	if (newStateVectors.rows <= nv) {
		i0 = 0;
		iN = newStateVectors.rows - 1;
	}
	else {
		i0 = max((int)((azimuthTime - newStateVectors.at<double>(0, 0)) / dt) - nv / 2 + 1, 0);
		iN = min(i0 + nv - 1, newStateVectors.rows - 1);
		i0 = (iN < newStateVectors.rows - 1 ? i0 : iN - nv + 1);
	}
	position.x = 0.0;
	position.y = 0.0;
	position.z = 0.0;
	for (int i = i0; i <= iN; ++i) {
		double weight = 1;
		for (int j = i0; j <= iN; ++j) {
			if (j != i) {
				double time2 = newStateVectors.at<double>(j, 0);
				weight *= (azimuthTime - time2) / (newStateVectors.at<double>(i, 0) - time2);
			}
		}
		position.x += weight * newStateVectors.at<double>(i, 1);
		position.y += weight * newStateVectors.at<double>(i, 2);
		position.z += weight * newStateVectors.at<double>(i, 3);
	}
	return 0;
}

int orbitStateVectors::getVelocity(double azimuthTime, Velocity& velocity)
{
	if (newStateVectors.cols != 7 || newStateVectors.rows < 2 || !isOrbitUpdated)
	{
		fprintf(stderr, "getVelocity(): input check failed!\n");
		return -1;
	}
	if (azimuthTime < newStateVectors.at<double>(0, 0) || azimuthTime > newStateVectors.at<double>(newStateVectors.rows - 1, 0))
	{
		fprintf(stderr, "getVelocity(): azimuthTime out of legal range\n");
		return -1;
	}
	int i0, iN;
	if (newStateVectors.rows <= nv) {
		i0 = 0;
		iN = newStateVectors.rows - 1;
	}
	else {
		i0 = max((int)((azimuthTime - newStateVectors.at<double>(0, 0)) / dt) - nv / 2 + 1, 0);
		iN = min(i0 + nv - 1, newStateVectors.rows - 1);
		i0 = (iN < newStateVectors.rows - 1 ? i0 : iN - nv + 1);
	}
	velocity.vx = 0.0;
	velocity.vy = 0.0;
	velocity.vz = 0.0;
	for (int i = i0; i <= iN; ++i) {
		double weight = 1.0;
		for (int j = i0; j <= iN; ++j) {
			if (j != i) {
				double time2 = newStateVectors.at<double>(j, 0);
				weight *= (azimuthTime - time2) / (newStateVectors.at<double>(i, 0) - time2);
			}
		}
		velocity.vx += weight * newStateVectors.at<double>(i, 4);
		velocity.vy += weight * newStateVectors.at<double>(i, 5);
		velocity.vz += weight * newStateVectors.at<double>(i, 6);
	}
	return 0;
}

int orbitStateVectors::getOrbitData(double time, OSV* osv)
{
	if (!osv || time < 0 || stateVectors.empty()|| isOrbitUpdated)
	{
		fprintf(stderr, "getOrbitData(): input check failed!\n");
		return -1;
	}
	int ret;
	int numVectors = stateVectors.rows;
	double t0 = stateVectors.at<double>(0, 0);
	double tN = stateVectors.at<double>(numVectors - 1, 0);

	int numVecPolyFit = polyDegree + 1; //4;
	int halfNumVecPolyFit = numVecPolyFit / 2;
	Mat vectorIndices = Mat::zeros(1, numVecPolyFit, CV_32S);
	int vecIdx = (int)((time - t0) / (tN - t0) * (numVectors - 1));
	if (vecIdx <= halfNumVecPolyFit - 1) {
		for (int i = 0; i < numVecPolyFit; i++) {
			vectorIndices.at<int>(0, i) = i;
		}
	}
	else if (vecIdx >= numVectors - halfNumVecPolyFit) {
		for (int i = 0; i < numVecPolyFit; i++) {
			vectorIndices.at<int>(0, i) = numVectors - numVecPolyFit + i;
		}
	}
	else {
		for (int i = 0; i < numVecPolyFit; i++) {
			vectorIndices.at<int>(0, i) = vecIdx - halfNumVecPolyFit + 1 + i;
		}
	}

	Mat timeArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat xPosArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat yPosArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat zPosArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat xVelArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat yVelArray = Mat::zeros(numVecPolyFit, 1, CV_64F);
	Mat zVelArray = Mat::zeros(numVecPolyFit, 1, CV_64F);


	for (int i = 0; i < numVecPolyFit; i++) {
		timeArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 0) - t0;
		xPosArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 1);
		yPosArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 2);
		zPosArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 3);
		xVelArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 4);
		yVelArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 5);
		zVelArray.at<double>(i, 0) = stateVectors.at<double>(vectorIndices.at<int>(0, i), 6);

	}
	Mat A, xPosCoeff, yPosCoeff, zPosCoeff, xVelCoeff, yVelCoeff, zVelCoeff;
	Utils::createVandermondeMatrix(timeArray, A, polyDegree);
	ret = Utils::ployFit(A, xPosArray, xPosCoeff);
	if (ret < 0) return -1;
	ret = Utils::ployFit(A, yPosArray, yPosCoeff); 
	if (ret < 0) return -1;
	ret = Utils::ployFit(A, zPosArray, zPosCoeff);
	if (ret < 0) return -1;
	ret = Utils::ployFit(A, xVelArray, xVelCoeff); 
	if (ret < 0) return -1;
	ret = Utils::ployFit(A, yVelArray, yVelCoeff);
	if (ret < 0) return -1;
	ret = Utils::ployFit(A, zVelArray, zVelCoeff);
	if (ret < 0) return -1;
	double normalizedTime = time - t0;

	osv->time = time;
	ret = Utils::polyVal(xPosCoeff, normalizedTime, &osv->x);
	if (ret < 0) return -1;
	ret = Utils::polyVal(yPosCoeff, normalizedTime, &osv->y);
	if (ret < 0) return -1;
	ret = Utils::polyVal(zPosCoeff, normalizedTime, &osv->z);
	if (ret < 0) return -1;
	ret = Utils::polyVal(xVelCoeff, normalizedTime, &osv->vx);
	if (ret < 0) return -1;
	ret = Utils::polyVal(yVelCoeff, normalizedTime, &osv->vy);
	if (ret < 0) return -1;
	ret = Utils::polyVal(zVelCoeff, normalizedTime, &osv->vz);
	if (ret < 0) return -1;
	return 0;
}

int orbitStateVectors::applyOrbit()
{
	if (isOrbitUpdated) return 0;
	double delta_t = 1.0;//1.0s
	this->dt = delta_t;
	double extra = 10.0;//10.0s
	double start = startTime - extra;
	double stop = stopTime + extra;
	int numVectors = (int)((stop - start) / delta_t);
	OSV osv;
	Mat newStateVectors = Mat::zeros(numVectors, 7, CV_64F);
	int ret;
	for (int i = 0; i < numVectors; i++)
	{
		ret = getOrbitData(start + (double)i * delta_t, &osv);
		if (ret < 0) return -1;
		newStateVectors.at<double>(i, 0) = start + (double)i * delta_t;
		newStateVectors.at<double>(i, 1) = osv.x;
		newStateVectors.at<double>(i, 2) = osv.y;
		newStateVectors.at<double>(i, 3) = osv.z;
		newStateVectors.at<double>(i, 4) = osv.vx;
		newStateVectors.at<double>(i, 5) = osv.vy;
		newStateVectors.at<double>(i, 6) = osv.vz;


	}
	newStateVectors.copyTo(this->newStateVectors);
	isOrbitUpdated = true;
	return 0;
}
