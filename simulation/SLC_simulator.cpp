#include"..\include\SLC_simulator.h"
#include"..\include\Utils.h"
#include"..\include\FormatConversion.h"
#include"..\include\ComplexMat.h"
#include"..\include\Filter.h"
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
SLC_simulator::SLC_simulator()
{
	memset(this->error_head, 0, 256);
	strcpy(this->error_head, "SLC_SIMULATOR_ERROR: error happens when using ");
}

SLC_simulator::~SLC_simulator()
{

}

int SLC_simulator::reflectivity(Mat& incidenceAngle, Mat& sigma)
{
	if (incidenceAngle.empty() || (incidenceAngle.type() != CV_64F && incidenceAngle.type() != CV_32F))
	{
		fprintf(stderr, "reflectivity(): input check failed!\n");
		return -1;
	}
	int rows = incidenceAngle.rows;
	int cols = incidenceAngle.cols;
	incidenceAngle.copyTo(sigma);
	if (incidenceAngle.type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double temp = incidenceAngle.at<double>(i, j);
				if (temp <= 17.0) sigma.at<double>(i, j) = 2 - 0.7058823529412 * temp;
				else if (temp > 17.0 && temp <= 30.0) sigma.at<double>(i, j) = -10 - 0.3846153846154 * (temp - 17.0);
				else if (temp > 30.0 && temp <= 80.0) sigma.at<double>(i, j) = -15 - 0.2 * (temp - 30);
				else if (temp > 80.0 && temp <= 90.0) sigma.at<double>(i, j) = -25 - 7.5 * (temp - 80);
				else sigma.at<double>(i, j) = -200.0;
				sigma.at<double>(i, j) = pow(10.0, sigma.at<double>(i, j) / 10.0);
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				double temp = incidenceAngle.at<float>(i, j);
				if (temp <= 17.0) sigma.at<float>(i, j) = 2 - 0.7058823529412 * temp;
				else if (temp > 17.0 && temp <= 30.0) sigma.at<float>(i, j) = -10 - 0.3846153846154 * (temp - 17.0);
				else if (temp > 30.0 && temp <= 80.0) sigma.at<float>(i, j) = -15 - 0.2 * (temp - 30);
				else if (temp > 80.0 && temp <= 90.0) sigma.at<float>(i, j) = -25 - 7.5 * (temp - 80);
				else sigma.at<float>(i, j) = -200.0;
				sigma.at<float>(i, j) = pow(10.0, sigma.at<float>(i, j) / 10.0);
			}
		}
	}

	return 0;
}

int SLC_simulator::computeIncidenceAngle(
	Mat& DEM_x,
	Mat& DEM_y,
	Mat& DEM_z,
	Mat& satellitePos_x,
	Mat& satellitePos_y,
	Mat& satellitePos_z,
	Mat& sigma
)
{
	if (DEM_x.size() != DEM_y.size() ||
		DEM_x.size() != DEM_z.size() ||
		DEM_x.size() != satellitePos_x.size() ||
		DEM_x.size() != satellitePos_y.size() ||
		DEM_x.size() != satellitePos_z.size() ||
		DEM_x.type() != DEM_y.type() ||
		DEM_x.type() != DEM_z.type() ||
		DEM_x.type() != satellitePos_x.type() ||
		DEM_x.type() != satellitePos_y.type() ||
		DEM_x.type() != satellitePos_z.type() ||
		DEM_x.rows < 1 || DEM_x.cols < 1 ||
		(DEM_x.type() != CV_64F && DEM_x.type() != CV_32F)
		)
	{
		fprintf(stderr, "computeIncidenceAngle(): input check failed!\n");
		return -1;
	}
	Utils util; int ret;
	int rows = DEM_x.rows;
	int cols = DEM_x.cols;
	Mat vec1, vec2, angle;
	if (DEM_x.type() == CV_64F)
	{
		vec1.create(rows * cols, 3, CV_64F), vec2.create(rows * cols, 3, CV_64F), angle.create(rows - 1, cols - 1, CV_64F);
	}
	else
	{
		vec1.create(rows * cols, 3, CV_32F), vec2.create(rows * cols, 3, CV_32F), angle.create(rows - 1, cols - 1, CV_32F);
	}
	vec2 = 0.0; vec1 = 0.0; angle = 0.0;
	if (DEM_x.type() == CV_64F)
	{
		for (int i = 1; i < rows; i++)
		{
			for (int j = 1; j < cols; j++)
			{
				double x0 = DEM_x.at<double>(i, j);
				double y0 = DEM_y.at<double>(i, j);
				double z0 = DEM_z.at<double>(i, j);

				double x1 = DEM_x.at<double>(i - 1, j);
				double y1 = DEM_y.at<double>(i - 1, j);
				double z1 = DEM_z.at<double>(i - 1, j);
				double x2 = DEM_x.at<double>(i - 1, j - 1);
				double y2 = DEM_y.at<double>(i - 1, j - 1);
				double z2 = DEM_z.at<double>(i - 1, j - 1);
				vec1.at<double>(i * cols + j, 0) = x0 - x1;
				vec1.at<double>(i * cols + j, 1) = y0 - y1;
				vec1.at<double>(i * cols + j, 2) = z0 - z1;

				vec2.at<double>(i * cols + j, 0) = x0 - x2;
				vec2.at<double>(i * cols + j, 1) = y0 - y2;
				vec2.at<double>(i * cols + j, 2) = z0 - z2;
			}
		}
	}
	else
	{
		for (int i = 1; i < rows; i++)
		{
			for (int j = 1; j < cols; j++)
			{
				double x0 = DEM_x.at<float>(i, j);
				double y0 = DEM_y.at<float>(i, j);
				double z0 = DEM_z.at<float>(i, j);

				double x1 = DEM_x.at<float>(i - 1, j);
				double y1 = DEM_y.at<float>(i - 1, j);
				double z1 = DEM_z.at<float>(i - 1, j);
				double x2 = DEM_x.at<float>(i - 1, j - 1);
				double y2 = DEM_y.at<float>(i - 1, j - 1);
				double z2 = DEM_z.at<float>(i - 1, j - 1);
				vec1.at<float>(i * cols + j, 0) = x0 - x1;
				vec1.at<float>(i * cols + j, 1) = y0 - y1;
				vec1.at<float>(i * cols + j, 2) = z0 - z1;

				vec2.at<float>(i * cols + j, 0) = x0 - x2;
				vec2.at<float>(i * cols + j, 1) = y0 - y2;
				vec2.at<float>(i * cols + j, 2) = z0 - z2;
			}
		}
	}
	ret = util.cross(vec1, vec2, vec1);
	if (return_check(ret, "cross()", error_head)) return -1;
	if (DEM_x.type() == CV_64F)
	{
		for (int i = 1; i < rows; i++)
		{
			for (int j = 1; j < cols; j++)
			{
				double sate_x = satellitePos_x.at<double>(i, j);
				double sate_y = satellitePos_y.at<double>(i, j);
				double sate_z = satellitePos_z.at<double>(i, j);
				double r1 = sate_x - DEM_x.at<double>(i, j);
				double r2 = sate_y - DEM_y.at<double>(i, j);
				double r3 = sate_z - DEM_z.at<double>(i, j);
				double r = sqrt(r1 * r1 + r2 * r2 + r3 * r3);
				double v1 = vec1.at<double>(i * cols + j, 0);
				double v2 = vec1.at<double>(i * cols + j, 1);
				double v3 = vec1.at<double>(i * cols + j, 2);
				double xx = v1 * DEM_x.at<double>(i, j) + v2 * DEM_y.at<double>(i, j) + v3 * DEM_z.at<double>(i, j);
				xx = xx / sqrt(xx * xx);
				double v = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
				double a = xx * (r1 * v1 + r2 * v2 + r3 * v3) / (r * v);
				angle.at<double>(i - 1, j - 1) = acos(a) / PI * 180.0;
			}
		}
	}
	else
	{
		for (int i = 1; i < rows; i++)
		{
			for (int j = 1; j < cols; j++)
			{
				double sate_x = satellitePos_x.at<float>(i, j);
				double sate_y = satellitePos_y.at<float>(i, j);
				double sate_z = satellitePos_z.at<float>(i, j);
				double r1 = sate_x - DEM_x.at<float>(i, j);
				double r2 = sate_y - DEM_y.at<float>(i, j);
				double r3 = sate_z - DEM_z.at<float>(i, j);
				double r = sqrt(r1 * r1 + r2 * r2 + r3 * r3);
				double v1 = vec1.at<float>(i * cols + j, 0);
				double v2 = vec1.at<float>(i * cols + j, 1);
				double v3 = vec1.at<float>(i * cols + j, 2);
				double xx = v1 * DEM_x.at<float>(i, j) + v2 * DEM_y.at<float>(i, j) + v3 * DEM_z.at<float>(i, j);
				xx = xx / sqrt(xx * xx);
				double v = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
				double a = xx * (r1 * v1 + r2 * v2 + r3 * v3) / (r * v);
				angle.at<float>(i - 1, j - 1) = acos(a) / PI * 180.0;
			}
		}
	}
	
	ret = reflectivity(angle, sigma);
	if (return_check(ret, "reflectivity", error_head)) return -1;
	cv::copyMakeBorder(sigma, sigma, 1, 0, 1, 0, cv::BORDER_REFLECT);
	return 0;
}

int SLC_simulator::generateSLC(
	Mat& stateVec,
	Mat& dem, 
	double lon_upperleft,
	double lat_upperleft,
	int sceneHeight,
	int sceneWidth,
	double nearRange,
	double prf, 
	double wavelength,
	double rangeSpacing,
	double azimuthSpacing,
	double acquisitionStartTime,
	double acquisitionStopTime,
	double SNR,
	ComplexMat& slc
)
{
	if (stateVec.cols != 7 ||
		stateVec.rows < 7 ||
		stateVec.type() != CV_64F ||
		dem.type() != CV_16S ||
		dem.empty() ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		sceneHeight < 1 ||
		sceneWidth < 1 ||
		nearRange <= 0.0 ||
		prf <= 0.0 ||
		wavelength <= 0.0 ||
		rangeSpacing <= 0.0 ||
		azimuthSpacing <= 0.0 ||
		acquisitionStartTime <= 0.0 ||
		acquisitionStartTime >= acquisitionStopTime
		)
	{
		fprintf(stderr, "generateSLC(): input check failed!\n");
		return -1;
	}
	slc.re.create(sceneHeight, sceneWidth, CV_32F);
	slc.im.create(sceneHeight, sceneWidth, CV_32F);
	slc.re = 0.0; slc.im = 0.0;
	//Mat pixel_count(sceneHeight, sceneWidth, CV_8U); pixel_count = 0;
	//分块计算，确定DEM划分大小与方式

	int interp_times_row = 90.0 / azimuthSpacing;
	int interp_times_col = 90.0 / rangeSpacing;
	int interp_cell = 4;
	double lon_spacing_old = 5.0 / 6000.0 / (double)interp_times_col;
	double lat_spacing_old = 5.0 / 6000.0 / (double)interp_times_row;
	double lon_spacing = lon_spacing_old / (double)interp_cell;
	double lat_spacing = lat_spacing_old / (double)interp_cell;
	int rows = dem.rows * interp_times_row;
	int cols = dem.cols * interp_times_col;
	int block_rows = 1000;//分块大小
	int block_cols = 1000;
	Mat dem_interp, dem_temp;
	dem.convertTo(dem_interp, CV_32F);
	cv::resize(dem_interp, dem_interp, cv::Size(cols, rows), 0, 0, cv::INTER_CUBIC);
	int num_block_row = rows / block_rows;
	int num_block_col = cols / block_cols;
	block_rows = rows / num_block_row;
	block_cols = cols / num_block_col;
	Utils util;
	num_block_row = 2;
	num_block_col = num_block_col;
	//初始化轨道类
	orbitStateVectors stateVectors(stateVec, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	for (int i = 10; i < 12; i++)
	{
		for (int j = 0; j < num_block_col; j++)
		{
			seed++;
			int row_start = i * block_rows - 1;
			row_start = row_start < 0 ? 0 : row_start;
			int row_end = (i + 1) * block_rows;
			row_end = row_end > rows ? rows : row_end;
			int col_start = j * block_cols - 1;
			col_start = col_start < 0 ? 0 : col_start;
			int col_end = (j + 1) * block_cols;
			col_end = col_end > cols ? cols : col_end;
			Mat dem_temp2;
			dem_interp(cv::Range(row_start, row_end), cv::Range(col_start, col_end)).copyTo(dem_temp);
			//dem_temp.convertTo(dem_temp, CV_32F);
			cv::resize(dem_temp, dem_temp2, cv::Size(dem_temp.cols * interp_cell, dem_temp.rows * interp_cell),
				0, 0, cv::INTER_CUBIC);
			
			//dem_temp2.convertTo(dem_temp2, CV_64F);
			//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\dem.bin", dem_temp2); return 0;
			//dem_temp2.convertTo(dem_temp2, CV_32F);
			double upper_left_lon = lon_upperleft + col_start * lon_spacing_old;
			upper_left_lon = upper_left_lon > 180.0 ? upper_left_lon - 360.0 : upper_left_lon;
			double upper_left_lat = lat_upperleft - row_start * lat_spacing_old;

			//生成随机相位
			Mat randomAngle(dem_temp2.rows, dem_temp2.cols, CV_32F);
			cv::RNG rng(seed);
			rng.fill(randomAngle, cv::RNG::UNIFORM, 0.0, 2.0 * PI);

			//加入热噪声项
			Mat noise_real(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary(dem_temp2.rows, dem_temp2.cols, CV_32F);
			double noise_sigma = sqrt(pow(10.0, -SNR / 10.0) / 2.0);
			cv::RNG rng2;
			rng2.fill(noise_real, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng3;
			rng3.fill(noise_imaginary, cv::RNG::NORMAL, 0, noise_sigma);

			//计算DEM点的成像卫星位置
			Mat /*satepos_x(dem_temp2.rows, dem_temp2.cols, CV_64F), satepos_y(dem_temp2.rows, dem_temp2.cols, CV_64F),
				satepos_z(dem_temp2.rows, dem_temp2.cols, CV_64F), DEM_x(dem_temp2.rows, dem_temp2.cols, CV_64F),
				DEM_y(dem_temp2.rows, dem_temp2.cols, CV_64F), DEM_z(dem_temp2.rows, dem_temp2.cols, CV_64F),*/
				imaging_time(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range(dem_temp2.rows, dem_temp2.cols, CV_64F);
			int DEM_rows = dem_temp2.rows;
			int DEM_cols = dem_temp2.cols;
#pragma omp parallel for schedule(guided)
			for (int ii = 0; ii < DEM_rows; ii++)
			{
				for (int jj = 0; jj < DEM_cols; jj++)
				{
					Position groundPosition;
					double lat, lon, height;
					lat = upper_left_lat - (double)ii * lat_spacing;
					lon = upper_left_lon + (double)jj * lon_spacing;
					lon = lon > 180.0 ? (lon - 360.0) : lon;
					height = dem_temp2.at<float>(ii, jj);
					Utils::ell2xyz(lon, lat, height, groundPosition);
					//DEM_x.at<double>(ii, jj) = groundPosition.x;
					//DEM_y.at<double>(ii, jj) = groundPosition.y;
					//DEM_z.at<double>(ii, jj) = groundPosition.z;
					int numOrbitVec = stateVectors.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors.newStateVectors.at<double>(iii, 1), stateVectors.newStateVectors.at<double>(iii, 2),
							stateVectors.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors.newStateVectors.at<double>(iii, 4), stateVectors.newStateVectors.at<double>(iii, 5),
							stateVectors.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors.newStateVectors.at<double>(iii, 0);
							secondVecFreq = currentFreq;
							break;
						}
					}

					if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
						continue;
					}

					double lowerBoundTime = firstVecTime;
					double upperBoundTime = secondVecTime;
					double lowerBoundFreq = firstVecFreq;
					double upperBoundFreq = secondVecFreq;
					double midTime, midFreq;
					double diffTime = fabs(upperBoundTime - lowerBoundTime);
					double absLineTimeInterval = time_interval;

					int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
					int numIterations = 0; Position pos; Velocity vel;
					while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

						midTime = (upperBoundTime + lowerBoundTime) / 2.0;
						stateVectors.getPosition(midTime, pos);
						stateVectors.getVelocity(midTime, vel);
						xdiff = groundPosition.x - pos.x;
						ydiff = groundPosition.y - pos.y;
						zdiff = groundPosition.z - pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength * distance);
						if ((midFreq - dopplerFrequency) * (lowerBoundFreq - dopplerFrequency) > 0.0) {
							lowerBoundTime = midTime;
							lowerBoundFreq = midFreq;
						}
						else if ((midFreq - dopplerFrequency) * (upperBoundFreq - dopplerFrequency) > 0.0) {
							upperBoundTime = midTime;
							upperBoundFreq = midFreq;
						}
						else if (fabs(midFreq - dopplerFrequency) < 0.01) {
							zeroDopplerTime = midTime;
							break;
						}

						diffTime = fabs(upperBoundTime - lowerBoundTime);
						numIterations++;
					}
					zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
					//satepos_x.at<double>(ii, jj) = pos.x;
					//satepos_y.at<double>(ii, jj) = pos.y;
					//satepos_z.at<double>(ii, jj) = pos.z;
					imaging_time.at<double>(ii, jj) = zeroDopplerTime;
					slant_range.at<double>(ii, jj) = distance;
				}
			}

			//计算RCS
			//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\slant_range.bin", slant_range);
			//return 0;
			//Mat sigma;
			//sigma.create(dem_temp2.rows, dem_temp2.cols, CV_64F); sigma = 0.5;
			//ret = computeIncidenceAngle(DEM_x, DEM_y, DEM_z, satepos_x, satepos_y, satepos_z, sigma);
			//if (return_check(ret, "computeIncidenceAngle()", error_head)) return -1;
			
			//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\sigma.bin", sigma); return 0;
			//逐点计算DEM像素点的实部虚部
			for (int ii = 0; ii < DEM_rows; ii++)
			{
				for (int jj = 0; jj < DEM_cols; jj++)
				{
					double zeroDopplerTime = imaging_time.at<double>(ii, jj);
					double distance = slant_range.at<double>(ii, jj);
					int azimuthIndex = round((zeroDopplerTime - acquisitionStartTime) / time_interval);
					int rangeIndex = round((distance - nearRange) / rangeSpacing);
					if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
					{

					}
					else
					{
						//pixel_count.at<uchar>(azimuthIndex, rangeIndex) += 1;
						double theta = -4.0 * PI * distance / wavelength + randomAngle.at<float>(ii, jj);
						double real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real.at<float>(ii, jj));
						double imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary.at<float>(ii, jj));
						slc.re.at<float>(azimuthIndex, rangeIndex) += real;
						slc.im.at<float>(azimuthIndex, rangeIndex) += imaginary;
					}
				}
			}

			fprintf(stdout, "process %lf: %d / %d\n", (double)seed / (double)(num_block_row * num_block_col) * 100.0, seed,
				num_block_row * num_block_col);
		}
	}

//#pragma omp parallel for schedule(guided)
//	for (int i = 0; i < sceneHeight; i++)
//	{
//		for (int j = 0; j < sceneWidth; j++)
//		{
//			if (pixel_count.at<uchar>(i, j) > 0)
//			{
//				slc.re.at<float>(i, j) = slc.re.at<float>(i, j) / (double)pixel_count.at<uchar>(i, j);
//				slc.im.at<float>(i, j) = slc.im.at<float>(i, j) / (double)pixel_count.at<uchar>(i, j);
//			}
//		}
//	}
	//pixel_count.convertTo(pixel_count, CV_64F);
	//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\pixel_count.bin", pixel_count);
	return 0;
}


