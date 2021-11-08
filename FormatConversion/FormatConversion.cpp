#include"pch.h"
#include"..\include\FormatConversion.h"
#include"gdal_priv.h"
#ifdef _DEBUG
#pragma comment(lib,"ComplexMat_d.lib")
#else
#pragma comment(lib,"ComplexMat.lib")
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
	//拍摄起始时间
	ret = xmldoc.get_str_para("startTimeUTC", acquisition_start_time);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_start_time", acquisition_start_time.c_str());
	if (return_check(ret, "write_str_to_h5()", error_head)) return -1;
	//拍摄结束时间
	ret = xmldoc.get_str_para("stopTimeUTC", acquisition_stop_time);
	if (return_check(ret, "get_str_para()", error_head)) return -1;
	ret = write_str_to_h5(dst_h5_filename, "acquisition_stop_time", acquisition_stop_time.c_str());
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
	//载频
	TiXmlElement* pnode, * pchild;
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
	ret = xmldoc._find_node(pnode, "slantRange", pchild);
	if (return_check(ret, "_find_node()", error_head)) return -1;
	ret = sscanf(pchild->GetText(), "%lf", &range_spacing);
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
		if (gps_time >= stop_time && fabs(gps_time - start_time) >= 100.0) stop = true;

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

int FormatConversion::read_slc_from_Sentinel(const char* filename, const char* xml_filename, ComplexMat& slc)
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

	/*
	* 为读取slc预分配内存
	*/
	slc.re.create(linesPerBurst * burst_count, samplesPerBurst, CV_16S);
	slc.im.create(linesPerBurst * burst_count, samplesPerBurst, CV_16S);

	INT16* buf = NULL;
	buf = (INT16*)malloc(linesPerBurst * samplesPerBurst * 2 * sizeof(INT16));
	if (!buf)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): out of memory!\n");
		return -1;
	}
	//逐个burst读取内容
	TiXmlElement* pchild = NULL;
	ret = xmldoc._find_node(pnode, "burst", pchild);
	if (ret < 0)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): node 'burst' not found!\n");
		if (buf) free(buf);
		return -1;
	}
	int64 bytesoffset;
	FILE* fp = NULL;
	fopen_s(&fp, filename, "rb");
	if (!fp)
	{
		fprintf(stderr, "read_slc_from_Sentinel(): failed to open %s!\n", filename);
		if (buf) free(buf);
		return -1;
	}
	for (int i = 0; i < burst_count; i++)
	{
		if (!pchild) break;
		ret = xmldoc._find_node(pchild, "byteOffset", pnode);
		if (ret < 0)
		{
			fprintf(stderr, "read_slc_from_Sentinel(): node 'byteOffset' not found!\n");
			if (buf) free(buf);
			if (fp) fclose(fp);
			return -1;
		}
		ret = sscanf(pnode->GetText(), "%lld", &bytesoffset);
		if (ret != 1)
		{
			fprintf(stderr, "read_slc_from_Sentinel(): %s: unknown data format!\n", xml_filename);
			if (buf) free(buf);
			if (fp) fclose(fp);
			return -1;
		}
		fseek(fp, bytesoffset, SEEK_SET);
		fread(buf, sizeof(INT16), linesPerBurst * samplesPerBurst * 2, fp);
		size_t offset = 0;
		for (int j = 0; j < linesPerBurst; j++)
		{
			for (int k = 0; k < samplesPerBurst; k++)
			{
				slc.re.at<short>(i * linesPerBurst + j, k) = buf[offset];
				offset++;
				slc.im.at<short>(i * linesPerBurst + j, k) = buf[offset];
				offset++;
			}
		}
		pchild = pchild->NextSiblingElement();
	}
	if (buf)free(buf);
	if (fp)fclose(fp);
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

	ComplexMat slc;Mat sentinel;
	int rows, cols;
	ret = read_slc_from_Sentinel(tiff_filename, xml_filename, slc);//需要deburst
	if (return_check(ret, "read_slc_from_Sentinel()", error_head)) return -1;
	ret = sentinel_deburst(xml_filename, slc, sentinel);
	if (return_check(ret, "sentinel_deburst()", error_head)) return -1;
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
		if (xx > 0)
		{
			gcps.at<double>(i, 3) = gcps.at<double>(i, 3) - (double)sentinel.at<int>(xx - 1, 0);
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
		int lon_topleft, lat_topleft, lon_bottomright, lat_bottomright, center_lat, center_lon;
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

		////中心经纬度
		//center_lat = (lat_topleft + lat_bottomright) / 2;
		//center_lon = (lon_topleft + lon_bottomright) / 2;

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
	//最近斜距
	//fseek(fp, 720 + 727 - 1, SEEK_SET);
	//memset(str, 0, 2048);
	//fread(str, 1, 16, fp);
	//tmp.at<double>(0, 0) = strtod(str, &ptr) / 1000000 * 3e8 * 0.5;
	//write_array_to_h5(dst_h5, "slant_range_first_pixel", tmp);
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
	this->doc.SaveFile(filename.c_str());
	return 0;
}

int XMLFile::XMLFile_add_origin(const char* node_name, const char* node_path)
{
	if (node_name == NULL ||
		node_path == NULL
		)
	{
		fprintf(stderr, "XMLFile_add_origin(): input check failed!\n");
		return -1;
	}
	TiXmlElement* DataNode = NULL;
	TiXmlElement* p = doc.RootElement();
	int ret = find_node_with_attribute(p, "DataNode", "name", "Origin_Images", DataNode);
	if (!DataNode)
	{
		DataNode = new TiXmlElement("DataNode");
		doc.RootElement()->LinkEndChild(DataNode);
		DataNode->SetAttribute("name", "Origin_Images");
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
		Sensor->LinkEndChild(new TiXmlText("TerraSAR-X"));
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











int XMLFile::XMLFile_add_cut(const char* datanode_name, const char* node_name, const char* node_path, int Row_offset, int Col_offset, double lon, double lat, double width, double height)
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
		doc.RootElement()->LinkEndChild(DataNode);
		DataNode->SetAttribute("name", datanode_name);
		DataNode->SetAttribute("index", "2");
		DataNode->SetAttribute("data_count", "1");
		DataNode->SetAttribute("data_processing", "cut");
		DataNode->SetAttribute("rank", "complex-1.0");
		TiXmlElement* Data = new TiXmlElement("Data");
		DataNode->LinkEndChild(Data);
		TiXmlElement* Data_Name = new TiXmlElement("Data_Name");
		Data_Name->LinkEndChild(new TiXmlText(node_name));
		Data->LinkEndChild(Data_Name);
		TiXmlElement* Data_Rank = new TiXmlElement("Data_Rank");
		Data_Rank->LinkEndChild(new TiXmlText("complex-1.0"));
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
		Data_Rank->LinkEndChild(new TiXmlText("complex-1.0"));
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
		ret = _find_node(pnode, "timeGPS", pchild);
		if (ret < 0)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
		ret = sscanf(pchild->GetText(), "%lf", &GPS_time);
		if (ret != 1)
		{
			fprintf(stderr, "get_stateVec_from_TSX(): %s : unknown data format!\n", this->m_xmlFileName);
			return -1;
		}
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
	return 0;
}