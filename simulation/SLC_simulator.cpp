#include"..\include\SLC_simulator.h"
#include"..\include\Utils.h"
#include"..\include\FormatConversion.h"
#include"..\include\ComplexMat.h"
#include"..\include\Filter.h"
#include"..\include\Deflat.h"
#ifdef _DEBUG
#pragma comment(lib, "ComplexMat_d.lib")
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "FormatConversion_d.lib")
#pragma comment(lib, "Deflat_d.lib")
#else
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "FormatConversion.lib")
#pragma comment(lib, "Deflat.lib")
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
enum ConvolutionType {
	/* Return the full convolution, including border */
	CONVOLUTION_FULL,

	/* Return only the part that corresponds to the original image */
	CONVOLUTION_SAME,
	/* Return only the submatrix containing elements that were not influenced by the border */
	CONVOLUTION_VALID
};
Mat conv2(const Mat& img, const Mat& ikernel, ConvolutionType type)
{
	Mat dest;
	Mat kernel;
	flip(ikernel, kernel, -1);
	Mat source = img;
	if (CONVOLUTION_FULL == type)
	{
		source = Mat();
		const int additionalRows = kernel.rows - 1, additionalCols = kernel.cols - 1;
		copyMakeBorder(img, source, (additionalRows + 1) / 2, additionalRows / 2, (additionalCols + 1) / 2, additionalCols / 2, cv::BORDER_CONSTANT, cv::Scalar(0));
	}
	cv::Point anchor(kernel.cols - kernel.cols / 2 - 1, kernel.rows - kernel.rows / 2 - 1);
	int borderMode = cv::BORDER_CONSTANT;
	filter2D(source, dest, img.depth(), kernel, anchor, 0, borderMode);

	if (CONVOLUTION_VALID == type)
	{
		dest = dest.colRange((kernel.cols - 1) / 2, dest.cols - kernel.cols / 2).rowRange((kernel.rows - 1) / 2, dest.rows - kernel.rows / 2);
	}
	return dest;
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
	ComplexMat& slc,
	Mat& GCP
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
	//分块计算，确定DEM划分大小与方式

	int interp_times_row = 90.0 / azimuthSpacing;
	int interp_times_col = 90.0 / rangeSpacing;
	int interp_cell = 2;
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
	num_block_row = num_block_row;
	num_block_col = num_block_col;
	vector<double> GCPs;//控制点信息
	//初始化轨道类
	orbitStateVectors stateVectors(stateVec, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	for (int i = 0; i < num_block_row; i++)
	{
		for (int j = 0; j < num_block_col; j++)
		{
			bool flag = (i % 3 == 0) && (j % 3 == 0) && i != 0 && j != 0;
			if (!flag) continue;
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
			cv::resize(dem_temp, dem_temp2, cv::Size(dem_temp.cols * interp_cell, dem_temp.rows * interp_cell),
				0, 0, cv::INTER_CUBIC);
			
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
			cv::RNG rng2(seed + 10000);
			rng2.fill(noise_real, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng3(seed + 10001);
			rng3.fill(noise_imaginary, cv::RNG::NORMAL, 0, noise_sigma);

			//计算DEM点的成像卫星位置
			Mat imaging_time(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range(dem_temp2.rows, dem_temp2.cols, CV_64F);
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
					while (diffTime > absLineTimeInterval * 0.01 && numIterations <= totalIterations) {

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
						else if (fabs(midFreq - dopplerFrequency) < 0.0001) {
							zeroDopplerTime = midTime;
							break;
						}

						diffTime = fabs(upperBoundTime - lowerBoundTime);
						numIterations++;
					}
					zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
					imaging_time.at<double>(ii, jj) = zeroDopplerTime;
					slant_range.at<double>(ii, jj) = distance;
				}
			}
			//控制点信息
			if ((i % 3 == 0) && (j % 3 == 0) && i != 0 && j != 0)
			{
				int gcp_row = DEM_rows / 2;
				int gcp_col = DEM_cols / 2;
				double gcp_sigma = 1000.0;
				double zeroDopplerTime = imaging_time.at<double>(gcp_row, gcp_col);
				double distance = slant_range.at<double>(gcp_row, gcp_col);
				int azimuthIndex = floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
				int rangeIndex = floor((distance - nearRange) / rangeSpacing);
				if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
				{

				}
				else
				{
					double lat, lon, height, lat_total = 0.0, lon_total = 0.0, height_total = 0.0, distance_total = 0.0;
					int azimuthIndex_tmp, rangeIndex_tmp, pixel_count = 0;
					for (int ii = 0; ii < DEM_rows; ii++)
					{
						for (int jj = 0; jj < DEM_cols; jj++)
						{
							zeroDopplerTime = imaging_time.at<double>(ii, jj);
							distance = slant_range.at<double>(ii, jj);
							azimuthIndex_tmp = floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
							rangeIndex_tmp = floor((distance - nearRange) / rangeSpacing);

							if (azimuthIndex_tmp == azimuthIndex && rangeIndex_tmp == (rangeIndex - 50))
							{
								lat = upper_left_lat - (double)ii * lat_spacing;
								lon = upper_left_lon + (double)jj * lon_spacing;
								lon = lon > 180.0 ? (lon - 360.0) : lon;
								height = dem_temp2.at<float>(ii, jj);
								pixel_count++;
								lon_total += lon;
								lat_total += lat;
								height_total += height;
								distance_total += distance;
							}
						}
					}
					if (pixel_count > 0)
					{
						lon_total /= (double)pixel_count;
						lat_total /= (double)pixel_count;
						height_total /= (double)pixel_count;
						distance_total /= (double)pixel_count;
					}
					GCPs.push_back(azimuthIndex + 1);
					GCPs.push_back(rangeIndex - 50 + 1);
					GCPs.push_back(lon_total);
					GCPs.push_back(lat_total);
					GCPs.push_back(height_total);
					GCPs.push_back(distance_total);
					GCPs.push_back(distance_total);
				}
			}
			

			printf("\rprocess %lf", (double)(i * num_block_col + j + 1) / (double)(num_block_row * num_block_col) * 100.0);
			fflush(stdout);
		}
	}

	int total_rows = GCPs.size() / 7;
	if (total_rows > 0)
	{
		GCP.create(total_rows, 7, CV_64F);

		for (int i = 0; i < total_rows; i++)
		{
			GCP.at<double>(i, 0) = GCPs[i * 7];
			GCP.at<double>(i, 1) = GCPs[i * 7 + 1];
			GCP.at<double>(i, 2) = GCPs[i * 7 + 2];
			GCP.at<double>(i, 3) = GCPs[i * 7 + 3];
			GCP.at<double>(i, 4) = GCPs[i * 7 + 4];
			GCP.at<double>(i, 5) = GCPs[i * 7 + 5];
			GCP.at<double>(i, 6) = GCPs[i * 7 + 6];
		}
	}
	return 0;
}

int SLC_simulator::generateSLC(
	Mat& stateVec1,
	Mat& stateVec2,
	Mat& dem,
	double lon_upperleft,
	double lat_upperleft,
	int sceneHeight1,
	int sceneWidth1,
	int sceneHeight2,
	int sceneWidth2,
	double nearRange1,
	double nearRange2,
	double prf,
	double wavelength,
	double rangeSpacing,
	double azimuthSpacing,
	double acquisitionStartTime1,
	double acquisitionStopTime1,
	double acquisitionStartTime2,
	double acquisitionStopTime2,
	double SNR,
	ComplexMat& slc1,
	ComplexMat& slc2,
	Mat& GCP,
	int mode
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec2.type() != CV_64F ||
		dem.type() != CV_16S ||
		dem.empty() ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		sceneHeight1 < 1 ||
		sceneWidth1 < 1 ||
		nearRange1 <= 0.0 ||
		sceneHeight2 < 1 ||
		sceneWidth2 < 1 ||
		nearRange2 <= 0.0 ||
		prf <= 0.0 ||
		wavelength <= 0.0 ||
		rangeSpacing <= 0.0 ||
		azimuthSpacing <= 0.0 ||
		acquisitionStartTime1 <= 0.0 ||
		acquisitionStartTime1 >= acquisitionStopTime1 ||
		acquisitionStartTime2 <= 0.0 ||
		acquisitionStartTime2 >= acquisitionStopTime2
		)
	{
		fprintf(stderr, "generateSLC(): input check failed!\n");
		return -1;
	}
	slc1.re.create(sceneHeight1, sceneWidth1, CV_32F);
	slc1.im.create(sceneHeight1, sceneWidth1, CV_32F);
	slc1.re = 0.0; slc1.im = 0.0;

	slc2.re.create(sceneHeight2, sceneWidth2, CV_32F);
	slc2.im.create(sceneHeight2, sceneWidth2, CV_32F);
	slc2.re = 0.0; slc2.im = 0.0;
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
	//int block_rows = 1000;//分块大小
	//int block_cols = 1000;
	int block_rows = 500;//分块大小
	int block_cols = 500;
	Mat dem_interp, dem_temp;
	dem.convertTo(dem_interp, CV_32F);
	cv::resize(dem_interp, dem_interp, cv::Size(cols, rows), 0, 0, cv::INTER_CUBIC);
	int num_block_row = rows / block_rows;
	int num_block_col = cols / block_cols;
	block_rows = rows / num_block_row;
	block_cols = cols / num_block_col;
	Utils util;
	num_block_row = num_block_row;
	num_block_col = num_block_col;
	vector<double> GCPs, GCPs2;//控制点信息
	//初始化轨道类
	orbitStateVectors stateVectors1(stateVec1, acquisitionStartTime1, acquisitionStopTime1);
	stateVectors1.applyOrbit();
	orbitStateVectors stateVectors2(stateVec2, acquisitionStartTime2, acquisitionStopTime2);
	stateVectors2.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	for (int i = 0; i < num_block_row; i++)
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
			cv::resize(dem_temp, dem_temp2, cv::Size(dem_temp.cols * interp_cell, dem_temp.rows * interp_cell),
				0, 0, cv::INTER_CUBIC);
			double upper_left_lon = lon_upperleft + col_start * lon_spacing_old;
			upper_left_lon = upper_left_lon > 180.0 ? upper_left_lon - 360.0 : upper_left_lon;
			double upper_left_lat = lat_upperleft - row_start * lat_spacing_old;

			//生成随机相位
			Mat randomAngle(dem_temp2.rows, dem_temp2.cols, CV_32F);
			cv::RNG rng(seed);
			rng.fill(randomAngle, cv::RNG::UNIFORM, 0.0, 2.0 * PI);

			//加入热噪声项
			Mat noise_real(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary(dem_temp2.rows, dem_temp2.cols, CV_32F);
			Mat noise_real2(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary2(dem_temp2.rows, dem_temp2.cols, CV_32F);
			double noise_sigma = sqrt(pow(10.0, -SNR / 10.0) / 2.0);
			cv::RNG rng2(seed + num_block_row * num_block_col + 1);
			rng2.fill(noise_real, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng3(seed + num_block_row * num_block_col + 2);
			rng3.fill(noise_imaginary, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng4(seed + num_block_row * num_block_col + 3);
			rng4.fill(noise_real2, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng5(seed + num_block_row * num_block_col + 4);
			rng5.fill(noise_imaginary2, cv::RNG::NORMAL, 0, noise_sigma);

			//计算DEM点的主图像成像卫星位置
			Mat imaging_time1(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range1(dem_temp2.rows, dem_temp2.cols, CV_64F);
			Mat imaging_time2(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range2(dem_temp2.rows, dem_temp2.cols, CV_64F);
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
					int numOrbitVec = stateVectors1.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors1.newStateVectors.at<double>(iii, 1), stateVectors1.newStateVectors.at<double>(iii, 2),
							stateVectors1.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors1.newStateVectors.at<double>(iii, 4), stateVectors1.newStateVectors.at<double>(iii, 5),
							stateVectors1.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
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
						stateVectors1.getPosition(midTime, pos);
						stateVectors1.getVelocity(midTime, vel);
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
					imaging_time1.at<double>(ii, jj) = zeroDopplerTime;
					slant_range1.at<double>(ii, jj) = distance;
				}
			}

			//计算DEM点的辅图像成像卫星位置
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
					int numOrbitVec = stateVectors2.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors2.newStateVectors.at<double>(iii, 1), stateVectors2.newStateVectors.at<double>(iii, 2),
							stateVectors2.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors2.newStateVectors.at<double>(iii, 4), stateVectors2.newStateVectors.at<double>(iii, 5),
							stateVectors2.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
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
						stateVectors2.getPosition(midTime, pos);
						stateVectors2.getVelocity(midTime, vel);
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
					imaging_time2.at<double>(ii, jj) = zeroDopplerTime;
					slant_range2.at<double>(ii, jj) = distance;
				}
			}

			//逐点计算DEM像素点的实部虚部
			if (mode == 2)//单发双收
			{
				for (int ii = 0; ii < DEM_rows; ii++)
				{
					for (int jj = 0; jj < DEM_cols; jj++)
					{
						double real, imaginary, theta;
						double zeroDopplerTime1 = imaging_time1.at<double>(ii, jj);
						double distance1 = slant_range1.at<double>(ii, jj);
						int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
						int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);
						if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1)
						{

						}
						else
						{
							theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(ii, jj);
							real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real.at<float>(ii, jj));
							imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary.at<float>(ii, jj));
							slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
							slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;
						}

						double zeroDopplerTime2 = imaging_time2.at<double>(ii, jj);
						double distance2 = slant_range2.at<double>(ii, jj);
						int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
						int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);
						if (azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
						{

						}
						else
						{
							theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(ii, jj);
							real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real2.at<float>(ii, jj));
							imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary2.at<float>(ii, jj));
							slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
							slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;
						}
					}
				}
				//控制点信息
				if ((i % 3 == 0) && (j % 3 == 0) && i != 0 && j != 0)
				{
					int gcp_row = DEM_rows / 2;
					int gcp_col = DEM_cols / 2;
					double gcp_sigma = 1000.0;
					double zeroDopplerTime1 = imaging_time1.at<double>(gcp_row, gcp_col);
					double distance1 = slant_range1.at<double>(gcp_row, gcp_col);
					int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
					int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);

					double zeroDopplerTime2 = imaging_time2.at<double>(gcp_row, gcp_col);
					double distance2 = slant_range2.at<double>(gcp_row, gcp_col);
					int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
					int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);

					if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1 ||
						azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
					{

					}
					else
					{
						double lat, lon, height;
						lat = upper_left_lat - (double)gcp_row * lat_spacing;
						lon = upper_left_lon + (double)gcp_col * lon_spacing;
						lon = lon > 180.0 ? (lon - 360.0) : lon;
						height = dem_temp2.at<float>(gcp_row, gcp_col);
						GCPs.push_back(azimuthIndex1 + 1);
						GCPs.push_back(rangeIndex1 + 1);
						GCPs.push_back(lon);
						GCPs.push_back(lat);
						GCPs.push_back(height);
						GCPs.push_back(distance1);
						GCPs.push_back(distance2);
						double theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
						double real = gcp_sigma * (cos(theta) + noise_real.at<float>(gcp_row, gcp_col));
						double imaginary = gcp_sigma * (sin(theta) + noise_imaginary.at<float>(gcp_row, gcp_col));
						slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
						slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;

						theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
						real = gcp_sigma * (cos(theta) + noise_real2.at<float>(gcp_row, gcp_col));
						imaginary = gcp_sigma * (sin(theta) + noise_imaginary2.at<float>(gcp_row, gcp_col));
						slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
						slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;
					}
				}
			}

			else//单发单收
			{
				for (int ii = 0; ii < DEM_rows; ii++)
				{
					for (int jj = 0; jj < DEM_cols; jj++)
					{
						double real, imaginary, theta;
						double zeroDopplerTime1 = imaging_time1.at<double>(ii, jj);
						double distance1 = slant_range1.at<double>(ii, jj);
						int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
						int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);
						if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1)
						{

						}
						else
						{
							theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(ii, jj);
							real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real.at<float>(ii, jj));
							imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary.at<float>(ii, jj));
							slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
							slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;
						}

						double zeroDopplerTime2 = imaging_time2.at<double>(ii, jj);
						double distance2 = slant_range2.at<double>(ii, jj);
						int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
						int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);
						if (azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
						{

						}
						else
						{
							theta = -4.0 * PI * distance2 / wavelength + randomAngle.at<float>(ii, jj);
							real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real2.at<float>(ii, jj));
							imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary2.at<float>(ii, jj));
							slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
							slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;
						}
					}
				}
				//控制点信息
				if ((i % 3 == 0) && (j % 3 == 0) && i != 0 && j != 0)
				{
					int gcp_row = DEM_rows / 2;
					int gcp_col = DEM_cols / 2;
					double gcp_sigma = 1000.0;
					double zeroDopplerTime1 = imaging_time1.at<double>(gcp_row, gcp_col);
					double distance1 = slant_range1.at<double>(gcp_row, gcp_col);
					int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
					int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);

					double zeroDopplerTime2 = imaging_time2.at<double>(gcp_row, gcp_col);
					double distance2 = slant_range2.at<double>(gcp_row, gcp_col);
					int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
					int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);

					if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1 ||
						azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
					{

					}
					else
					{
						double lat, lon, height;
						lat = upper_left_lat - (double)gcp_row * lat_spacing;
						lon = upper_left_lon + (double)gcp_col * lon_spacing;
						lon = lon > 180.0 ? (lon - 360.0) : lon;
						height = dem_temp2.at<float>(gcp_row, gcp_col);
						GCPs.push_back(azimuthIndex1 + 1);
						GCPs.push_back(rangeIndex1 + 1);
						GCPs.push_back(lon);
						GCPs.push_back(lat);
						GCPs.push_back(height);
						GCPs.push_back(distance1);
						GCPs.push_back(distance2);
						double theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
						double real = gcp_sigma * (cos(theta) + noise_real.at<float>(gcp_row, gcp_col));
						double imaginary = gcp_sigma * (sin(theta) + noise_imaginary.at<float>(gcp_row, gcp_col));
						slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
						slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;

						theta = -4.0 * PI * distance2 / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
						real = gcp_sigma * (cos(theta) + noise_real2.at<float>(gcp_row, gcp_col));
						imaginary = gcp_sigma * (sin(theta) + noise_imaginary2.at<float>(gcp_row, gcp_col));
						slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
						slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;
					}
				}
			}
			

			fprintf(stdout, "process %lf: %d / %d\n", (double)seed / (double)(num_block_row * num_block_col) * 100.0, seed,
				num_block_row * num_block_col);
		}
	}

	////二维复数卷积
	int wright = 16;
	double Drange = rangeSpacing;
	double Dazimuth = azimuthSpacing;
	double res_range = rangeSpacing * 3.0;
	double res_azimuth = azimuthSpacing * 3.0;

	Mat sinc(wright * 2 + 1, wright * 2 + 1, CV_32F);
	sinc = 0.0;
	for (int i = -wright; i <= wright; i++)
	{
		double bb;
		double r = Drange * i;
		if (i == 0) bb = 1.0;
		else
		{
			bb = sin(PI * r / res_range) / (PI * r / res_range);
		}
		for (int j = -wright; j <= wright; j++)
		{
			double aa;
			double x = Dazimuth * j;
			if (j == 0) aa = 1.0;
			else
			{
				aa = sin(PI * x / res_azimuth) / (PI * x / res_azimuth);
			}
			sinc.at<float>(i + wright, j + wright) = aa * bb;
		}
	}

	Mat x = conv2(slc1.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc1.re);
	x = conv2(slc1.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc1.im);

	x = conv2(slc2.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc2.re);
	x = conv2(slc2.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc2.im);

	int total_rows = GCPs.size() / 7;
	if (total_rows > 0)
	{
		GCP.create(total_rows, 7, CV_64F);

		for (int i = 0; i < total_rows; i++)
		{
			GCP.at<double>(i, 0) = GCPs[i * 7];
			GCP.at<double>(i, 1) = GCPs[i * 7 + 1];
			GCP.at<double>(i, 2) = GCPs[i * 7 + 2];
			GCP.at<double>(i, 3) = GCPs[i * 7 + 3];
			GCP.at<double>(i, 4) = GCPs[i * 7 + 4];
			GCP.at<double>(i, 5) = GCPs[i * 7 + 5];
			GCP.at<double>(i, 6) = GCPs[i * 7 + 6];
		}
	}
	return 0;
}

int SLC_simulator::generateSLC(
	Mat& stateVec1,
	Mat& stateVec2,
	Mat& dem,
	double lon_upperleft,
	double lat_upperleft,
	int sceneHeight1,
	int sceneWidth1,
	int sceneHeight2, 
	int sceneWidth2,
	double nearRange1, 
	double nearRange2,
	double prf,
	double wavelength, 
	double rangeSpacing, 
	double azimuthSpacing, 
	double acquisitionStartTime1, 
	double acquisitionStopTime1,
	double acquisitionStartTime2, 
	double acquisitionStopTime2,
	double SNR, 
	ComplexMat& slc1,
	ComplexMat& slc2,
	ComplexMat& slc3, 
	ComplexMat& slc4, 
	Mat& GCP1,
	Mat& GCP2
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec2.type() != CV_64F ||
		dem.type() != CV_16S ||
		dem.empty() ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		sceneHeight1 < 1 ||
		sceneWidth1 < 1 ||
		nearRange1 <= 0.0 ||
		sceneHeight2 < 1 ||
		sceneWidth2 < 1 ||
		nearRange2 <= 0.0 ||
		prf <= 0.0 ||
		wavelength <= 0.0 ||
		rangeSpacing <= 0.0 ||
		azimuthSpacing <= 0.0 ||
		acquisitionStartTime1 <= 0.0 ||
		acquisitionStartTime1 >= acquisitionStopTime1 ||
		acquisitionStartTime2 <= 0.0 ||
		acquisitionStartTime2 >= acquisitionStopTime2
		)
	{
		fprintf(stderr, "generateSLC(): input check failed!\n");
		return -1;
	}
	slc1.re.create(sceneHeight1, sceneWidth1, CV_32F);
	slc1.im.create(sceneHeight1, sceneWidth1, CV_32F);
	slc1.re = 0.0; slc1.im = 0.0;

	slc2.re.create(sceneHeight1, sceneWidth1, CV_32F);
	slc2.im.create(sceneHeight1, sceneWidth1, CV_32F);
	slc2.re = 0.0; slc2.im = 0.0;

	slc3.re.create(sceneHeight2, sceneWidth2, CV_32F);
	slc3.im.create(sceneHeight2, sceneWidth2, CV_32F);
	slc3.re = 0.0; slc3.im = 0.0;

	slc4.re.create(sceneHeight2, sceneWidth2, CV_32F);
	slc4.im.create(sceneHeight2, sceneWidth2, CV_32F);
	slc4.re = 0.0; slc4.im = 0.0;

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
	//int block_rows = 1000;//分块大小
	//int block_cols = 1000;
	int block_rows = 500;//分块大小
	int block_cols = 500;
	Mat dem_interp, dem_temp;
	dem.convertTo(dem_interp, CV_32F);
	cv::resize(dem_interp, dem_interp, cv::Size(cols, rows), 0, 0, cv::INTER_CUBIC);
	int num_block_row = rows / block_rows;
	int num_block_col = cols / block_cols;
	block_rows = rows / num_block_row;
	block_cols = cols / num_block_col;
	Utils util;
	num_block_row = num_block_row;
	num_block_col = num_block_col;
	vector<double> GCPs, GCPs2;//控制点信息
	//初始化轨道类
	orbitStateVectors stateVectors1(stateVec1, acquisitionStartTime1, acquisitionStopTime1);
	stateVectors1.applyOrbit();
	orbitStateVectors stateVectors2(stateVec2, acquisitionStartTime2, acquisitionStopTime2);
	stateVectors2.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	char process[512];
	for (int i = 0; i < num_block_row; i++)
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
			cv::resize(dem_temp, dem_temp2, cv::Size(dem_temp.cols * interp_cell, dem_temp.rows * interp_cell),
				0, 0, cv::INTER_CUBIC);
			double upper_left_lon = lon_upperleft + col_start * lon_spacing_old;
			upper_left_lon = upper_left_lon > 180.0 ? upper_left_lon - 360.0 : upper_left_lon;
			double upper_left_lat = lat_upperleft - row_start * lat_spacing_old;

			//生成随机相位
			Mat randomAngle(dem_temp2.rows, dem_temp2.cols, CV_32F);
			cv::RNG rng(seed);
			rng.fill(randomAngle, cv::RNG::UNIFORM, 0.0, 2.0 * PI);

			//加入热噪声项
			Mat noise_real(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary(dem_temp2.rows, dem_temp2.cols, CV_32F);
			Mat noise_real2(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary2(dem_temp2.rows, dem_temp2.cols, CV_32F);
			Mat noise_real3(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary3(dem_temp2.rows, dem_temp2.cols, CV_32F);
			Mat noise_real4(dem_temp2.rows, dem_temp2.cols, CV_32F), noise_imaginary4(dem_temp2.rows, dem_temp2.cols, CV_32F);
			double noise_sigma = sqrt(pow(10.0, -SNR / 10.0) / 2.0);
			cv::RNG rng2(seed + 1 * num_block_row * num_block_col + 1);
			rng2.fill(noise_real, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng3(seed + 1 * num_block_row * num_block_col + 2);
			rng3.fill(noise_imaginary, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng4(seed + 1 * num_block_row * num_block_col + 3);
			rng4.fill(noise_real2, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng5(seed + 1 * num_block_row * num_block_col + 4);
			rng5.fill(noise_imaginary2, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng6(seed + 1 * num_block_row * num_block_col + 5);
			rng6.fill(noise_real3, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng7(seed + 1 * num_block_row * num_block_col + 6);
			rng7.fill(noise_imaginary3, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng8(seed + 1 * num_block_row * num_block_col + 7);
			rng8.fill(noise_real4, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng9(seed + 1 * num_block_row * num_block_col + 8);
			rng9.fill(noise_imaginary4, cv::RNG::NORMAL, 0, noise_sigma);
			//noise_imaginary4 = 0.0;
			//noise_imaginary3 = 0.0;
			//noise_imaginary2 = 0.0;
			//noise_imaginary = 0.0;
			//noise_real4 = 0.0;
			//noise_real3 = 0.0;
			//noise_real2 = 0.0;
			//noise_real = 0.0;
			//计算DEM点的主图像成像卫星位置
			Mat imaging_time1(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range1(dem_temp2.rows, dem_temp2.cols, CV_64F);
			Mat imaging_time2(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range2(dem_temp2.rows, dem_temp2.cols, CV_64F);
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
					int numOrbitVec = stateVectors1.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors1.newStateVectors.at<double>(iii, 1), stateVectors1.newStateVectors.at<double>(iii, 2),
							stateVectors1.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors1.newStateVectors.at<double>(iii, 4), stateVectors1.newStateVectors.at<double>(iii, 5),
							stateVectors1.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
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
						stateVectors1.getPosition(midTime, pos);
						stateVectors1.getVelocity(midTime, vel);
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
					imaging_time1.at<double>(ii, jj) = zeroDopplerTime;
					slant_range1.at<double>(ii, jj) = distance;
				}
			}

			//计算DEM点的辅图像成像卫星位置
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
					int numOrbitVec = stateVectors2.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors2.newStateVectors.at<double>(iii, 1), stateVectors2.newStateVectors.at<double>(iii, 2),
							stateVectors2.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors2.newStateVectors.at<double>(iii, 4), stateVectors2.newStateVectors.at<double>(iii, 5),
							stateVectors2.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
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
						stateVectors2.getPosition(midTime, pos);
						stateVectors2.getVelocity(midTime, vel);
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
					imaging_time2.at<double>(ii, jj) = zeroDopplerTime;
					slant_range2.at<double>(ii, jj) = distance;
				}
			}

			//逐点计算DEM像素点的实部虚部
			for (int ii = 0; ii < DEM_rows; ii++)
			{
				for (int jj = 0; jj < DEM_cols; jj++)
				{
					double real, imaginary, theta;
					double zeroDopplerTime1 = imaging_time1.at<double>(ii, jj);
					double distance1 = slant_range1.at<double>(ii, jj);
					int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
					int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);

					double zeroDopplerTime2 = imaging_time2.at<double>(ii, jj);
					double distance2 = slant_range2.at<double>(ii, jj);
					int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
					int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);

					if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1)
					{

					}
					else
					{
						theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(ii, jj);
						real = 1.0 * (cos(theta) + noise_real.at<float>(ii, jj));
						imaginary = 1.0 * (sin(theta) + noise_imaginary.at<float>(ii, jj));
						slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
						slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;


						theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(ii, jj);
						real = 1.0 * (cos(theta) + noise_real4.at<float>(ii, jj));
						imaginary = 1.0 * (sin(theta) + noise_imaginary4.at<float>(ii, jj));
						slc4.re.at<float>(azimuthIndex1, rangeIndex1) += real;
						slc4.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;
					}

					
					if (azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
					{

					}
					else
					{
						theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(ii, jj);
						real = 1.0 * (cos(theta) + noise_real2.at<float>(ii, jj));
						imaginary = 1.0 * (sin(theta) + noise_imaginary2.at<float>(ii, jj));
						slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
						slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;


						theta = -4.0 * PI * distance2 / wavelength + randomAngle.at<float>(ii, jj);
						real = 1.0 * (cos(theta) + noise_real3.at<float>(ii, jj));
						imaginary = 1.0 * (sin(theta) + noise_imaginary3.at<float>(ii, jj));
						slc3.re.at<float>(azimuthIndex2, rangeIndex2) += real;
						slc3.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;

						
					}

				}
			}
			//控制点信息
			if ((i % 3 == 0) && (j % 3 == 0) && i != 0 && j != 0)
			{
				int gcp_row = DEM_rows / 2;
				int gcp_col = DEM_cols / 2;
				double gcp_sigma = 1000.0;
				double zeroDopplerTime1 = imaging_time1.at<double>(gcp_row, gcp_col);
				double distance1 = slant_range1.at<double>(gcp_row, gcp_col);
				int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
				int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);

				double zeroDopplerTime2 = imaging_time2.at<double>(gcp_row, gcp_col);
				double distance2 = slant_range2.at<double>(gcp_row, gcp_col);
				int azimuthIndex2 = floor((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
				int rangeIndex2 = floor((distance2 - nearRange2) / rangeSpacing);

				if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1 ||
					azimuthIndex2 < 0 || azimuthIndex2 > sceneHeight2 - 1 || rangeIndex2 < 0 || rangeIndex2 > sceneWidth2 - 1)
				{

				}
				else
				{
					double lat, lon, height;
					lat = upper_left_lat - (double)gcp_row * lat_spacing;
					lon = upper_left_lon + (double)gcp_col * lon_spacing;
					lon = lon > 180.0 ? (lon - 360.0) : lon;
					height = dem_temp2.at<float>(gcp_row, gcp_col);
					GCPs.push_back(azimuthIndex1 + 1);
					GCPs.push_back(rangeIndex1 + 1);
					GCPs.push_back(lon);
					GCPs.push_back(lat);
					GCPs.push_back(height);
					GCPs.push_back(distance1);
					GCPs.push_back(distance2);

					GCPs2.push_back(azimuthIndex2 + 1);
					GCPs2.push_back(rangeIndex2 + 1);
					GCPs2.push_back(lon);
					GCPs2.push_back(lat);
					GCPs2.push_back(height);
					GCPs2.push_back(distance2);
					GCPs2.push_back(distance1);

					double theta = -4.0 * PI * distance1 / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					double real = gcp_sigma * (cos(theta) + noise_real.at<float>(gcp_row, gcp_col));
					double imaginary = gcp_sigma * (sin(theta) + noise_imaginary.at<float>(gcp_row, gcp_col));
					slc1.re.at<float>(azimuthIndex1, rangeIndex1) += real;
					slc1.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;

					theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					real = gcp_sigma * (cos(theta) + noise_real2.at<float>(gcp_row, gcp_col));
					imaginary = gcp_sigma * (sin(theta) + noise_imaginary2.at<float>(gcp_row, gcp_col));
					slc2.re.at<float>(azimuthIndex2, rangeIndex2) += real;
					slc2.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;

					theta = -4.0 * PI * distance2 / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					real = gcp_sigma * (cos(theta) + noise_real3.at<float>(gcp_row, gcp_col));
					imaginary = gcp_sigma * (sin(theta) + noise_imaginary3.at<float>(gcp_row, gcp_col));
					slc3.re.at<float>(azimuthIndex2, rangeIndex2) += real;
					slc3.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;

					theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					real = gcp_sigma * (cos(theta) + noise_real4.at<float>(gcp_row, gcp_col));
					imaginary = gcp_sigma * (sin(theta) + noise_imaginary4.at<float>(gcp_row, gcp_col));
					slc4.re.at<float>(azimuthIndex1, rangeIndex1) += real;
					slc4.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;
				}
			}
			printf("\rprocess %lf: %d / %d", (double)seed / (double)(num_block_row * num_block_col) * 100.0, seed,
				num_block_row* num_block_col);
			fflush(stdout);
		}
	}

	////二维复数卷积
	int wright = 16;
	double Drange = rangeSpacing;
	double Dazimuth = azimuthSpacing;
	double res_range = rangeSpacing * 3.0;
	double res_azimuth = azimuthSpacing * 3.0;

	Mat sinc(wright * 2 + 1, wright * 2 + 1, CV_32F);
	sinc = 0.0;
	for (int i = -wright; i <= wright; i++)
	{
		double bb;
		double r = Drange * i;
		if (i == 0) bb = 1.0;
		else
		{
			bb = sin(PI * r / res_range) / (PI * r / res_range);
		}
		for (int j = -wright; j <= wright; j++)
		{
			double aa;
			double x = Dazimuth * j;
			if (j == 0) aa = 1.0;
			else
			{
				aa = sin(PI * x / res_azimuth) / (PI * x / res_azimuth);
			}
			sinc.at<float>(i + wright, j + wright) = aa * bb;
		}
	}

	Mat x = conv2(slc1.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc1.re);
	x = conv2(slc1.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc1.im);

	x = conv2(slc2.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc2.re);
	x = conv2(slc2.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc2.im);

	x = conv2(slc3.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc3.re);
	x = conv2(slc3.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc3.im);

	x = conv2(slc4.re, sinc, CONVOLUTION_SAME);
	x.copyTo(slc4.re);
	x = conv2(slc4.im, sinc, CONVOLUTION_SAME);
	x.copyTo(slc4.im);

	int total_rows = GCPs.size() / 7;
	if (total_rows > 0)
	{
		GCP1.create(total_rows, 7, CV_64F);

		for (int i = 0; i < total_rows; i++)
		{
			GCP1.at<double>(i, 0) = GCPs[i * 7];
			GCP1.at<double>(i, 1) = GCPs[i * 7 + 1];
			GCP1.at<double>(i, 2) = GCPs[i * 7 + 2];
			GCP1.at<double>(i, 3) = GCPs[i * 7 + 3];
			GCP1.at<double>(i, 4) = GCPs[i * 7 + 4];
			GCP1.at<double>(i, 5) = GCPs[i * 7 + 5];
			GCP1.at<double>(i, 6) = GCPs[i * 7 + 6];
		}
	}

	total_rows = GCPs2.size() / 7;
	if (total_rows > 0)
	{
		GCP2.create(total_rows, 7, CV_64F);

		for (int i = 0; i < total_rows; i++)
		{
			GCP2.at<double>(i, 0) = GCPs2[i * 7];
			GCP2.at<double>(i, 1) = GCPs2[i * 7 + 1];
			GCP2.at<double>(i, 2) = GCPs2[i * 7 + 2];
			GCP2.at<double>(i, 3) = GCPs2[i * 7 + 3];
			GCP2.at<double>(i, 4) = GCPs2[i * 7 + 4];
			GCP2.at<double>(i, 5) = GCPs2[i * 7 + 5];
			GCP2.at<double>(i, 6) = GCPs2[i * 7 + 6];
		}
	}
	return 0;
}

int SLC_simulator::generateSlantrange(
	Mat& stateVec1,
	Mat& stateVec2,
	Mat& dem,
	double lon_upperleft,
	double lat_upperleft,
	int sceneHeight1,
	int sceneWidth1,
	double nearRange1,
	double prf, 
	double wavelength,
	double rangeSpacing,
	double azimuthSpacing,
	double acquisitionStartTime1,
	double acquisitionStopTime1,
	double acquisitionStartTime2, 
	double acquisitionStopTime2,
	Mat& R1, 
	Mat& R2
)
{
	if (stateVec1.cols != 7 ||
		stateVec1.rows < 7 ||
		stateVec1.type() != CV_64F ||
		stateVec2.cols != 7 ||
		stateVec2.rows < 7 ||
		stateVec2.type() != CV_64F ||
		dem.type() != CV_16S ||
		dem.empty() ||
		fabs(lon_upperleft) > 180.0 ||
		fabs(lat_upperleft) > 90.0 ||
		sceneHeight1 < 1 ||
		sceneWidth1 < 1 ||
		nearRange1 <= 0.0 ||
		prf <= 0.0 ||
		wavelength <= 0.0 ||
		rangeSpacing <= 0.0 ||
		azimuthSpacing <= 0.0 ||
		acquisitionStartTime1 <= 0.0 ||
		acquisitionStartTime1 >= acquisitionStopTime1 ||
		acquisitionStartTime2 <= 0.0 ||
		acquisitionStartTime2 >= acquisitionStopTime2
		)
	{
		fprintf(stderr, "generateSlantrange(): input check failed!\n");
		return -1;
	}

	R1.create(sceneHeight1, sceneWidth1, CV_64F); R1 = 0.0;
	R2.create(sceneHeight1, sceneWidth1, CV_64F); R2 = 0.0;
	Mat pixel_count = Mat::zeros(sceneHeight1, sceneWidth1, CV_8U); pixel_count = 0;
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
	num_block_row = num_block_row;
	num_block_col = num_block_col;
	vector<double> GCPs;//控制点信息
	//初始化轨道类
	orbitStateVectors stateVectors1(stateVec1, acquisitionStartTime1, acquisitionStopTime1);
	stateVectors1.applyOrbit();
	orbitStateVectors stateVectors2(stateVec2, acquisitionStartTime2, acquisitionStopTime2);
	stateVectors2.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	char process[512];
	for (int i = 0; i < num_block_row; i++)
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
			cv::resize(dem_temp, dem_temp2, cv::Size(dem_temp.cols * interp_cell, dem_temp.rows * interp_cell),
				0, 0, cv::INTER_CUBIC);
			double upper_left_lon = lon_upperleft + col_start * lon_spacing_old;
			upper_left_lon = upper_left_lon > 180.0 ? upper_left_lon - 360.0 : upper_left_lon;
			double upper_left_lat = lat_upperleft - row_start * lat_spacing_old;


			//计算DEM点的主图像成像卫星位置
			Mat imaging_time1(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range1(dem_temp2.rows, dem_temp2.cols, CV_64F);
			Mat imaging_time2(dem_temp2.rows, dem_temp2.cols, CV_64F), slant_range2(dem_temp2.rows, dem_temp2.cols, CV_64F);
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
					int numOrbitVec = stateVectors1.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors1.newStateVectors.at<double>(iii, 1), stateVectors1.newStateVectors.at<double>(iii, 2),
							stateVectors1.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors1.newStateVectors.at<double>(iii, 4), stateVectors1.newStateVectors.at<double>(iii, 5),
							stateVectors1.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors1.newStateVectors.at<double>(iii, 0);
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
						stateVectors1.getPosition(midTime, pos);
						stateVectors1.getVelocity(midTime, vel);
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
					imaging_time1.at<double>(ii, jj) = zeroDopplerTime;
					slant_range1.at<double>(ii, jj) = distance;
				}
			}

			//计算DEM点的辅图像成像卫星位置
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
					int numOrbitVec = stateVectors2.newStateVectors.rows;
					double firstVecTime = 0.0;
					double secondVecTime = 0.0;
					double firstVecFreq = 0.0;
					double secondVecFreq = 0.0;
					double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
					for (int iii = 0; iii < numOrbitVec; iii++) {
						Position orb_pos(stateVectors2.newStateVectors.at<double>(iii, 1), stateVectors2.newStateVectors.at<double>(iii, 2),
							stateVectors2.newStateVectors.at<double>(iii, 3));
						Velocity orb_vel(stateVectors2.newStateVectors.at<double>(iii, 4), stateVectors2.newStateVectors.at<double>(iii, 5),
							stateVectors2.newStateVectors.at<double>(iii, 6));
						currentFreq = 0;
						xdiff = groundPosition.x - orb_pos.x;
						ydiff = groundPosition.y - orb_pos.y;
						zdiff = groundPosition.z - orb_pos.z;
						distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
						currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
						if (iii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
							firstVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
							firstVecFreq = currentFreq;
						}
						else {
							secondVecTime = stateVectors2.newStateVectors.at<double>(iii, 0);
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
						stateVectors2.getPosition(midTime, pos);
						stateVectors2.getVelocity(midTime, vel);
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
					imaging_time2.at<double>(ii, jj) = zeroDopplerTime;
					slant_range2.at<double>(ii, jj) = distance;
				}
			}

			//逐点计算DEM像素点的实部虚部
			for (int ii = 0; ii < DEM_rows; ii++)
			{
				for (int jj = 0; jj < DEM_cols; jj++)
				{
					double real, imaginary, theta;
					double zeroDopplerTime1 = imaging_time1.at<double>(ii, jj);
					double distance1 = slant_range1.at<double>(ii, jj);
					int azimuthIndex1 = floor((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
					int rangeIndex1 = floor((distance1 - nearRange1) / rangeSpacing);

					double distance2 = slant_range2.at<double>(ii, jj);

					if (azimuthIndex1 < 0 || azimuthIndex1 > sceneHeight1 - 1 || rangeIndex1 < 0 || rangeIndex1 > sceneWidth1 - 1)
					{

					}
					else
					{
						pixel_count.at<uchar>(azimuthIndex1, rangeIndex1) += 1;
						R1.at<double>(azimuthIndex1, rangeIndex1) += distance1;
						R2.at<double>(azimuthIndex1, rangeIndex1) += distance2;
					}
				}
			}
			printf("\rprocess: %d / %d", seed, num_block_row * num_block_col);
			fflush(stdout);
		}
	}

	for (int i = 0; i < sceneHeight1; i++)
	{
		for (int j = 0; j < sceneWidth1; j++)
		{
			if (pixel_count.at<uchar>(i, j) > 0)
			{
				R1.at<double>(i, j) = R1.at<double>(i, j) / (double)pixel_count.at<uchar>(i, j);
				R2.at<double>(i, j) = R2.at<double>(i, j) / (double)pixel_count.at<uchar>(i, j);
			}
		}
	}

	return 0;
}

int SLC_simulator::SLC_deramp(
	Mat& mappedDEM, 
	Mat& mappedLat, 
	Mat& mappedLon, 
	const char* slcH5File1,
	const char* slcH5File2,
	const char* slcH5File3,
	const char* slcH5File4, 
	const char* slcH5File1_out,
	const char* slcH5File2_out, 
	const char* slcH5File3_out, 
	const char* slcH5File4_out
)
{
	if (mappedDEM.rows != mappedLat.rows ||
		mappedDEM.rows != mappedLon.rows ||
		mappedDEM.cols != mappedLat.cols ||
		mappedDEM.cols != mappedLon.cols ||
		mappedDEM.type() != CV_16S ||
		mappedLat.type() != CV_32F ||
		mappedLon.type() != CV_32F ||
		mappedDEM.empty() ||
		!slcH5File1 || !slcH5File1_out || !slcH5File2 || !slcH5File2_out || !slcH5File3 ||
		!slcH5File3_out || !slcH5File4 || !slcH5File4_out
		)
	{
		fprintf(stderr, "SLC_deramp(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	ComplexMat slc;
	int ret;
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, prf, start, end, start2, end2;
	int sceneHeight, sceneWidth, sceneHeight2, sceneWidth2, offset_row = 0, offset_col = 0;
	Mat lon_coef, lat_coef, statevec, statevec2;
	string start_time, end_time;
	ret = conversion.read_int_from_h5(slcH5File1, "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File1, "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "range_len", &sceneWidth2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "azimuth_len", &sceneHeight2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_double_from_h5(slcH5File1, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File1, "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(slcH5File1, "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start2);
	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end2);
	ret = conversion.read_array_from_h5(slcH5File3, "state_vec", statevec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	


	Mat sate1 = Mat::zeros(sceneHeight, 3, CV_64F);

	Mat sate2 = Mat::zeros(sceneHeight2, 3, CV_64F);

	orbitStateVectors stateVectors(statevec, start, end);
	stateVectors.applyOrbit();

	orbitStateVectors stateVectors2(statevec2, start2, end2);
	stateVectors2.applyOrbit();

	//计算主星成像位置
	double dopplerFrequency = 0.0;
	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<float>(0, 0);
	lon = mappedLon.at<float>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
	}

	//计算辅星成像位置
	numOrbitVec = stateVectors2.newStateVectors.rows;
	firstVecTime = 0.0;
	secondVecTime = 0.0;
	firstVecFreq = 0.0;
	secondVecFreq = 0.0;
	currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors2.newStateVectors.at<double>(ii, 1), stateVectors2.newStateVectors.at<double>(ii, 2),
			stateVectors2.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors2.newStateVectors.at<double>(ii, 4), stateVectors2.newStateVectors.at<double>(ii, 5),
			stateVectors2.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	lowerBoundTime = firstVecTime;
	upperBoundTime = secondVecTime;
	lowerBoundFreq = firstVecFreq;
	upperBoundFreq = secondVecFreq;
	diffTime = fabs(upperBoundTime - lowerBoundTime);
	absLineTimeInterval = 1.0 / prf;

	totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	numIterations = 0;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors2.getPosition(midTime, pos);
		stateVectors2.getVelocity(midTime, vel);
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

	for (int i = 0; i < sceneHeight2; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors2.getPosition(time, pos);
		sate2.at<double>(i, 0) = pos.x;
		sate2.at<double>(i, 1) = pos.y;
		sate2.at<double>(i, 2) = pos.z;
	}

	//主星发主星收图像去参考
	ret = conversion.read_slc_from_h5(slcH5File1, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			r = -r / wavelength * 4 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File1_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File1_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5File1_out);
	ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File1_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File1_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File1_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File1_out, "azimuth_len", sceneHeight);
	//辅星发主星收图像去参考
	ret = conversion.read_slc_from_h5(slcH5File2, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r += cv::norm(tt, cv::NORM_L2);
			r = -r / wavelength * 2.0 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File2_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File2_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File2, slcH5File2_out);
	ret = conversion.read_int_from_h5(slcH5File2, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File2_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File2, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File2_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File2_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File2_out, "azimuth_len", sceneHeight);
	//辅星发辅星收图像去参考
	ret = conversion.read_slc_from_h5(slcH5File3, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			r = -r / wavelength * 4 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File3_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File3_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File3, slcH5File3_out);
	ret = conversion.read_int_from_h5(slcH5File3, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File3_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File3, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File3_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File3_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File3_out, "azimuth_len", sceneHeight);
	//主星发辅星收图像去参考
	ret = conversion.read_slc_from_h5(slcH5File4, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r += cv::norm(tt, cv::NORM_L2);
			r = -r / wavelength * 2.0 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File4_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File4_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File4, slcH5File4_out);
	ret = conversion.read_int_from_h5(slcH5File4, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File4_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File4, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File4_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File4_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File4_out, "azimuth_len", sceneHeight);
	return 0;
}

int SLC_simulator::SLC_deramp_14(
	vector<string>& slcH5FilesList,
	vector<string>& slcH5FilesListOut,
	int master_index,
	Mat& mappedDEM,
	Mat& mappedLat,
	Mat& mappedLon, 
	int mode
)
{
	if (mappedDEM.rows != mappedLat.rows ||
		mappedDEM.rows != mappedLon.rows ||
		mappedDEM.cols != mappedLat.cols ||
		mappedDEM.cols != mappedLon.cols ||
		mappedDEM.type() != CV_16S ||
		mappedLat.type() != CV_32F ||
		mappedLon.type() != CV_32F ||
		mappedDEM.empty() ||
		slcH5FilesList.size() != slcH5FilesListOut.size() ||
		slcH5FilesList.size() < 2 || slcH5FilesList.size() > 8 ||
		mode < 1 || mode > 4
		)
	{
		fprintf(stderr, "SLC_deramp_14(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	ComplexMat slc, slc2;
	int ret;
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, wavelength2, prf, start, end, start2, end2;
	int sceneHeight, sceneWidth, sceneHeight2, sceneWidth2, offset_row = 0, offset_col = 0;
	Mat lon_coef, lat_coef, statevec, statevec2;
	string start_time, end_time;
	const char* slcH5File1 = NULL, *slcH5File3 = NULL;
	//单发单收/单发双收模式
	if (mode == 1 || mode == 2)
	{
		slcH5File1 = slcH5FilesList[0].c_str();
		slcH5File3 = slcH5FilesList[1].c_str();
	}
	else if (mode == 3)//乒乓模式
	{
		slcH5File1 = slcH5FilesList[0].c_str();
		slcH5File3 = slcH5FilesList[3].c_str();
	}
	else//双频乒乓模式
	{
		slcH5File1 = slcH5FilesList[0].c_str();
		slcH5File3 = slcH5FilesList[7].c_str();
	}

	ret = conversion.read_int_from_h5(slcH5File1, "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File1, "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "range_len", &sceneWidth2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "azimuth_len", &sceneHeight2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_double_from_h5(slcH5File1, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;

	ret = conversion.read_double_from_h5(slcH5File1, "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	ret = conversion.read_double_from_h5(slcH5File3, "carrier_frequency", &wavelength2);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength2 = VEL_C / wavelength2;

	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(slcH5File1, "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start2);
	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end2);
	ret = conversion.read_array_from_h5(slcH5File3, "state_vec", statevec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	ret = conversion.read_int_from_h5(slcH5FilesList[master_index - 1].c_str(), "offset_row", &offset_row);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5FilesList[master_index - 1].c_str(), "offset_col", &offset_col);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;


	Mat sate1 = Mat::zeros(sceneHeight, 3, CV_64F);

	Mat sate2 = Mat::zeros(sceneHeight2, 3, CV_64F);

	orbitStateVectors stateVectors(statevec, start, end);
	stateVectors.applyOrbit();

	orbitStateVectors stateVectors2(statevec2, start2, end2);
	stateVectors2.applyOrbit();

	//计算主星成像位置
	double dopplerFrequency = 0.0;
	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<float>(0, 0);
	lon = mappedLon.at<float>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
	}

	//计算辅星成像位置
	numOrbitVec = stateVectors2.newStateVectors.rows;
	firstVecTime = 0.0;
	secondVecTime = 0.0;
	firstVecFreq = 0.0;
	secondVecFreq = 0.0;
	currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors2.newStateVectors.at<double>(ii, 1), stateVectors2.newStateVectors.at<double>(ii, 2),
			stateVectors2.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors2.newStateVectors.at<double>(ii, 4), stateVectors2.newStateVectors.at<double>(ii, 5),
			stateVectors2.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	lowerBoundTime = firstVecTime;
	upperBoundTime = secondVecTime;
	lowerBoundFreq = firstVecFreq;
	upperBoundFreq = secondVecFreq;
	diffTime = fabs(upperBoundTime - lowerBoundTime);
	absLineTimeInterval = 1.0 / prf;

	totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	numIterations = 0;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors2.getPosition(midTime, pos);
		stateVectors2.getVelocity(midTime, vel);
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

	for (int i = 0; i < sceneHeight2; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors2.getPosition(time, pos);
		sate2.at<double>(i, 0) = pos.x;
		sate2.at<double>(i, 1) = pos.y;
		sate2.at<double>(i, 2) = pos.z;
	}

	if (mode == 1)//单发单收模式
	{
		//主星图像去参考
		ret = conversion.read_slc_from_h5(slcH5File1, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		Mat R = Mat::zeros(sceneHeight, sceneWidth, CV_64F);

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[0].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[0].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[0].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5FilesListOut[0].c_str());
		ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "azimuth_len", sceneHeight);

		//辅星去参考
		R = 0.00;
		ret = conversion.read_slc_from_h5(slcH5File3, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[1].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(slcH5FilesListOut[1].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[1].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File3, slcH5FilesListOut[1].c_str());
		ret = conversion.read_int_from_h5(slcH5File3, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File3, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "azimuth_len", sceneHeight);
	}
	else if (mode == 2)//单发双收模式
	{
		//主星图像去参考
		ret = conversion.read_slc_from_h5(slcH5File1, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		Mat R = Mat::zeros(sceneHeight, sceneWidth, CV_64F);

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[0].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[0].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[0].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5FilesListOut[0].c_str());
		ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "azimuth_len", sceneHeight);

		//辅星去参考
		R = 0.00;
		ret = conversion.read_slc_from_h5(slcH5File3, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r += cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[1].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(slcH5FilesListOut[1].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[1].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File3, slcH5FilesListOut[1].c_str());
		ret = conversion.read_int_from_h5(slcH5File3, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File3, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "azimuth_len", sceneHeight);
	}
	else if (mode == 3)//乒乓模式
	{
		//主星发主星收图像去参考
		ret = conversion.read_slc_from_h5(slcH5File1, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		Mat R = Mat::zeros(sceneHeight, sceneWidth, CV_64F);

	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[0].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[0].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[0].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5FilesListOut[0].c_str());
		ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "azimuth_len", sceneHeight);


		//辅星发主星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5FilesList[1].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r += cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[1].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[1].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[1].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[1].c_str(), slcH5FilesListOut[1].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[1].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[1].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "azimuth_len", sceneHeight);
		//辅星发辅星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5File3, slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[3].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[3].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[3].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5File3, slcH5FilesListOut[3].c_str());
		ret = conversion.read_int_from_h5(slcH5File3, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File3, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "azimuth_len", sceneHeight);
		//主星发辅星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5FilesList[2].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r += cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[2].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[2].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[2].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[2].c_str(), slcH5FilesListOut[2].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[2].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[2].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "azimuth_len", sceneHeight);
	}
	else //双频乒乓模式
	{
		//主星发主星收图像去参考
		ret = conversion.read_slc_from_h5(slcH5FilesList[0].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		ret = conversion.read_slc_from_h5(slcH5FilesList[4].c_str(), slc2);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		if (slc2.type() != CV_32F) slc.convertTo(slc2, CV_32F);
		Mat R = Mat::zeros(sceneHeight, sceneWidth, CV_64F);

	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;

				r = -R.at<double>(i, j) / wavelength2 * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc2.re.at<float>(i, j);
				imagine2 = slc2.im.at<float>(i, j);
				slc2.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc2.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[0].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.creat_new_h5(slcH5FilesListOut[4].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[0].c_str(), "slantRange", R);
		ret = conversion.write_array_to_h5(slcH5FilesListOut[4].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[0].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[4].c_str(), slc2);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;


		ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5FilesListOut[0].c_str());
		ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[0].c_str(), "azimuth_len", sceneHeight);

		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[4].c_str(), slcH5FilesListOut[4].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[4].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[4].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[4].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[4].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[4].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[4].c_str(), "azimuth_len", sceneHeight);

		//辅星发主星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5FilesList[1].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		ret = conversion.read_slc_from_h5(slcH5FilesList[5].c_str(), slc2);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		if (slc2.type() != CV_32F) slc.convertTo(slc2, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r += cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;

				r = -R.at<double>(i, j) / wavelength2 * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc2.re.at<float>(i, j);
				imagine2 = slc2.im.at<float>(i, j);
				slc2.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc2.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[1].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.creat_new_h5(slcH5FilesListOut[5].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[1].c_str(), "slantRange", R);
		ret = conversion.write_array_to_h5(slcH5FilesListOut[5].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[1].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[5].c_str(), slc2);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;


		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[1].c_str(), slcH5FilesListOut[1].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[1].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[1].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[1].c_str(), "azimuth_len", sceneHeight);

		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[5].c_str(), slcH5FilesListOut[5].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[5].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[5].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[5].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[5].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[5].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[5].c_str(), "azimuth_len", sceneHeight);
		//辅星发辅星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5FilesList[3].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		ret = conversion.read_slc_from_h5(slcH5FilesList[7].c_str(), slc2);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		if (slc2.type() != CV_32F) slc.convertTo(slc2, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;


				r = -R.at<double>(i, j) / wavelength2 * 4 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc2.re.at<float>(i, j);
				imagine2 = slc2.im.at<float>(i, j);
				slc2.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc2.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[3].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.creat_new_h5(slcH5FilesListOut[7].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;

		ret = conversion.write_array_to_h5(slcH5FilesListOut[3].c_str(), "slantRange", R);
		ret = conversion.write_array_to_h5(slcH5FilesListOut[7].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[3].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[7].c_str(), slc2);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;


		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[3].c_str(), slcH5FilesListOut[3].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[3].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[3].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[3].c_str(), "azimuth_len", sceneHeight);

		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[7].c_str(), slcH5FilesListOut[7].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[7].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[7].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[7].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[7].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[7].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[7].c_str(), "azimuth_len", sceneHeight);
		//主星发辅星收图像去参考
		R = 0.0;
		ret = conversion.read_slc_from_h5(slcH5FilesList[2].c_str(), slc);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
		ret = conversion.read_slc_from_h5(slcH5FilesList[6].c_str(), slc2);
		if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
		if (slc2.type() != CV_32F) slc.convertTo(slc2, CV_32F);
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < sceneHeight; i++)
		{
			for (int j = 0; j < sceneWidth; j++)
			{
				double r, real, imagine, real2, imagine2;
				Mat XYZ, LLH(1, 3, CV_64F), tt;
				LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
				LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
				LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
				util.ell2xyz(LLH, XYZ);
				tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
				r = cv::norm(tt, cv::NORM_L2);
				tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
				r += cv::norm(tt, cv::NORM_L2);
				R.at<double>(i, j) = r;
				r = -r / wavelength * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc.re.at<float>(i, j);
				imagine2 = slc.im.at<float>(i, j);
				slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;

				r = -R.at<double>(i, j) / wavelength2 * 2.0 * PI;
				real = cos(r);
				imagine = sin(r);
				real2 = slc2.re.at<float>(i, j);
				imagine2 = slc2.im.at<float>(i, j);
				slc2.re.at<float>(i, j) = real * real2 + imagine * imagine2;
				slc2.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
			}
		}
		ret = conversion.creat_new_h5(slcH5FilesListOut[2].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;
		ret = conversion.creat_new_h5(slcH5FilesListOut[6].c_str());
		if (return_check(ret, "creat_new_h5()", error_head)) return -1;


		ret = conversion.write_array_to_h5(slcH5FilesListOut[2].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[2].c_str(), slc);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
		ret = conversion.write_array_to_h5(slcH5FilesListOut[6].c_str(), "slantRange", R);
		ret = conversion.write_slc_to_h5(slcH5FilesListOut[6].c_str(), slc2);
		if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;



		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[2].c_str(), slcH5FilesListOut[2].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[2].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[2].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[2].c_str(), "azimuth_len", sceneHeight);

		ret = conversion.Copy_para_from_h5_2_h5(slcH5FilesList[6].c_str(), slcH5FilesListOut[6].c_str());
		ret = conversion.read_int_from_h5(slcH5FilesList[6].c_str(), "offset_row", &offset_row);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[6].c_str(), "offset_row", offset_row);
		ret = conversion.read_int_from_h5(slcH5FilesList[6].c_str(), "offset_col", &offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[6].c_str(), "offset_col", offset_col);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[6].c_str(), "range_len", sceneWidth);
		ret = conversion.write_int_to_h5(slcH5FilesListOut[6].c_str(), "azimuth_len", sceneHeight);
	}
	return 0;
}

int SLC_simulator::SLC_reramp(
	Mat& mappedDEM,
	Mat& mappedLat,
	Mat& mappedLon,
	const char* slcH5File1,
	const char* slcH5File2,
	const char* slcH5File3,
	const char* slcH5File4,
	const char* slcH5File1_out,
	const char* slcH5File2_out,
	const char* slcH5File3_out,
	const char* slcH5File4_out
)
{
	if (mappedDEM.rows != mappedLat.rows ||
		mappedDEM.rows != mappedLon.rows ||
		mappedDEM.cols != mappedLat.cols ||
		mappedDEM.cols != mappedLon.cols ||
		mappedDEM.type() != CV_16S ||
		mappedLat.type() != CV_32F ||
		mappedLon.type() != CV_32F ||
		mappedDEM.empty() ||
		!slcH5File1 || !slcH5File1_out || !slcH5File2 || !slcH5File2_out || !slcH5File3 ||
		!slcH5File3_out || !slcH5File4 || !slcH5File4_out
		)
	{
		fprintf(stderr, "SLC_deramp(): input check failed!\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	ComplexMat slc;
	int ret;
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft, rangeSpacing,
		nearRangeTime, wavelength, prf, start, end, start2, end2;
	int sceneHeight, sceneWidth, sceneHeight2, sceneWidth2, offset_row = 0, offset_col = 0;
	Mat lon_coef, lat_coef, statevec, statevec2;
	string start_time, end_time;
	ret = conversion.read_int_from_h5(slcH5File1, "range_len", &sceneWidth);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File1, "azimuth_len", &sceneHeight);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "range_len", &sceneWidth2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(slcH5File3, "azimuth_len", &sceneHeight2);
	if (return_check(ret, "read_int_from_h5()", error_head)) return -1;

	ret = conversion.read_double_from_h5(slcH5File1, "prf", &prf);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	ret = conversion.read_double_from_h5(slcH5File1, "carrier_frequency", &wavelength);
	if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
	wavelength = VEL_C / wavelength;
	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start);
	ret = conversion.read_str_from_h5(slcH5File1, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end);
	ret = conversion.read_array_from_h5(slcH5File1, "state_vec", statevec);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_start_time", start_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.utc2gps(start_time.c_str(), &start2);
	ret = conversion.read_str_from_h5(slcH5File3, "acquisition_stop_time", end_time);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	conversion.utc2gps(end_time.c_str(), &end2);
	ret = conversion.read_array_from_h5(slcH5File3, "state_vec", statevec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;




	Mat sate1 = Mat::zeros(sceneHeight, 3, CV_64F);

	Mat sate2 = Mat::zeros(sceneHeight2, 3, CV_64F);

	orbitStateVectors stateVectors(statevec, start, end);
	stateVectors.applyOrbit();

	orbitStateVectors stateVectors2(statevec2, start2, end2);
	stateVectors2.applyOrbit();

	//计算主星成像位置
	double dopplerFrequency = 0.0;
	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<float>(0, 0);
	lon = mappedLon.at<float>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
	}

	//计算辅星成像位置
	numOrbitVec = stateVectors2.newStateVectors.rows;
	firstVecTime = 0.0;
	secondVecTime = 0.0;
	firstVecFreq = 0.0;
	secondVecFreq = 0.0;
	currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors2.newStateVectors.at<double>(ii, 1), stateVectors2.newStateVectors.at<double>(ii, 2),
			stateVectors2.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors2.newStateVectors.at<double>(ii, 4), stateVectors2.newStateVectors.at<double>(ii, 5),
			stateVectors2.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	lowerBoundTime = firstVecTime;
	upperBoundTime = secondVecTime;
	lowerBoundFreq = firstVecFreq;
	upperBoundFreq = secondVecFreq;
	diffTime = fabs(upperBoundTime - lowerBoundTime);
	absLineTimeInterval = 1.0 / prf;

	totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	numIterations = 0;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors2.getPosition(midTime, pos);
		stateVectors2.getVelocity(midTime, vel);
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

	for (int i = 0; i < sceneHeight2; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors2.getPosition(time, pos);
		sate2.at<double>(i, 0) = pos.x;
		sate2.at<double>(i, 1) = pos.y;
		sate2.at<double>(i, 2) = pos.z;
	}

	//主星发主星收图像加参考
	ret = conversion.read_slc_from_h5(slcH5File1, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			r = r / wavelength * 4 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File1_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File1_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File1, slcH5File1_out);
	ret = conversion.read_int_from_h5(slcH5File1, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File1_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File1, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File1_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File1_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File1_out, "azimuth_len", sceneHeight);
	//辅星发主星收图像加参考
	ret = conversion.read_slc_from_h5(slcH5File2, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r += cv::norm(tt, cv::NORM_L2);
			r = r / wavelength * 2.0 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File2_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File2_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File2, slcH5File2_out);
	ret = conversion.read_int_from_h5(slcH5File2, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File2_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File2, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File2_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File2_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File2_out, "azimuth_len", sceneHeight);
	//辅星发辅星收图像加参考
	ret = conversion.read_slc_from_h5(slcH5File3, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			r = r / wavelength * 4 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File3_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File3_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File3, slcH5File3_out);
	ret = conversion.read_int_from_h5(slcH5File3, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File3_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File3, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File3_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File3_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File3_out, "azimuth_len", sceneHeight);
	//主星发辅星收图像加参考
	ret = conversion.read_slc_from_h5(slcH5File4, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slc.type() != CV_32F) slc.convertTo(slc, CV_32F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r, real, imagine, real2, imagine2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			r = cv::norm(tt, cv::NORM_L2);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			r += cv::norm(tt, cv::NORM_L2);
			r = r / wavelength * 2.0 * PI;
			real = cos(r);
			imagine = sin(r);
			real2 = slc.re.at<float>(i, j);
			imagine2 = slc.im.at<float>(i, j);
			slc.re.at<float>(i, j) = real * real2 + imagine * imagine2;
			slc.im.at<float>(i, j) = real * imagine2 - real2 * imagine;
		}
	}
	ret = conversion.creat_new_h5(slcH5File4_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File4_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.Copy_para_from_h5_2_h5(slcH5File4, slcH5File4_out);
	ret = conversion.read_int_from_h5(slcH5File4, "offset_row", &offset_row);
	ret = conversion.write_int_to_h5(slcH5File4_out, "offset_row", offset_row);
	ret = conversion.read_int_from_h5(slcH5File4, "offset_col", &offset_col);
	ret = conversion.write_int_to_h5(slcH5File4_out, "offset_col", offset_col);
	ret = conversion.write_int_to_h5(slcH5File4_out, "range_len", sceneWidth);
	ret = conversion.write_int_to_h5(slcH5File4_out, "azimuth_len", sceneHeight);
	return 0;
}

int SLC_simulator::MB_phase_estimation(
	int estimation_wndsize,
	const char* slcH5File1,
	const char* slcH5File2,
	const char* slcH5File3,
	const char* slcH5File4,
	const char* slcH5File1_out,
	const char* slcH5File2_out,
	const char* slcH5File3_out,
	const char* slcH5File4_out
)
{
	if (estimation_wndsize % 2 == 0 || estimation_wndsize < 5 ||
		!slcH5File1 || !slcH5File1_out || !slcH5File2 ||
		!slcH5File2_out || !slcH5File3 || !slcH5File3_out ||
		!slcH5File4 || !slcH5File4_out
		)
	{
		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
		return -1;
	}
	vector<string>coregis_slc_files; vector<string>coregis_slc_files_out;
	coregis_slc_files.push_back(slcH5File1);
	coregis_slc_files.push_back(slcH5File2);
	coregis_slc_files.push_back(slcH5File3);
	coregis_slc_files.push_back(slcH5File4);
	coregis_slc_files_out.push_back(slcH5File1_out);
	coregis_slc_files_out.push_back(slcH5File2_out);
	coregis_slc_files_out.push_back(slcH5File3_out);
	coregis_slc_files_out.push_back(slcH5File4_out);
	
	FormatConversion conversion; Utils util;
	double thresh_c1_to_c2 = 0.7;
	int master_indx = 1;
	int n_images = 4, nr, nc, blocksize_row = 1000, blocksize_col = 1000;
	int homotest_radius = (estimation_wndsize - 1) / 2;
	int left, right, top, bottom, block_num_row, block_num_col, left_pad, right_pad, top_pad, bottom_pad, ret;
	vector<ComplexMat> slc_series, slc_series_filter;
	ComplexMat slc;
	Mat ph, phase;

	//预先填充
	ret = conversion.read_slc_from_h5(slcH5File1, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	nr = slc.GetRows(); nc = slc.GetCols();
	ret = conversion.creat_new_h5(slcH5File1_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.creat_new_h5(slcH5File2_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.creat_new_h5(slcH5File3_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	ret = conversion.creat_new_h5(slcH5File4_out);
	if (return_check(ret, "creat_new_h5()", error_head)) return -1;


	ret = conversion.write_slc_to_h5(slcH5File1_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File2_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File3_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	ret = conversion.write_slc_to_h5(slcH5File4_out, slc);
	if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;

	if (nr % blocksize_row == 0) block_num_row = nr / blocksize_row;
	else block_num_row = int(floor((double)nr / (double)blocksize_row)) + 1;
	if (nc % blocksize_col == 0) block_num_col = nc / blocksize_col;
	else block_num_col = int(floor((double)nc / (double)blocksize_col)) + 1;
	for (int i = 0; i < block_num_row; i++)
	{
		for (int j = 0; j < block_num_col; j++)
		{
			top = i * blocksize_row;
			top_pad = top - homotest_radius; top_pad = top_pad < 0 ? 0 : top_pad;
			bottom = top + blocksize_row; bottom = bottom > nr ? nr : bottom;
			bottom_pad = bottom + homotest_radius; bottom_pad = bottom_pad > nr ? nr : bottom_pad;
			left = j * blocksize_col;
			left_pad = left - homotest_radius; left_pad = left_pad < 0 ? 0 : left_pad;
			right = left + blocksize_col; right = right > nc ? nc : right;
			right_pad = right + homotest_radius; right_pad = right_pad > nc ? nc : right_pad;

			//读取数据
			for (int k = 0; k < n_images; k++)
			{
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_re",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.re);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_im",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.im);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				if (slc.type() != CV_64F) slc.convertTo(slc, CV_64F);
				slc_series.push_back(slc);
				slc_series_filter.push_back(slc);
			}
			//phase.create(slc.GetRows(), slc.GetCols(), CV_64F); phase = 0.0;
			//计算
#pragma omp parallel for schedule(guided)
			for (int ii = (top - top_pad); ii < (bottom - top_pad); ii++)
			{
				ComplexMat coherence_matrix, eigenvector; Mat eigenvalue; int ret1;
				Mat real1, imag1, real2, imag2, real, imag;
				double value1, value2;
				for (int jj = (left - left_pad); jj < (right - left_pad); jj++)
				{
					ret1 = util.coherence_matrix_estimation(slc_series, coherence_matrix, estimation_wndsize, estimation_wndsize, ii, jj, false, true);
					if (ret1 == 0)
					{
						ret1 = util.HermitianEVD(coherence_matrix, eigenvalue, eigenvector);
						if (!eigenvalue.empty() && ret1 == 0)
						{

							/*value1 = eigenvalue.at<double>(0, 0);
							value2 = eigenvalue.at<double>(1, 0);
							eigenvector.re(cv::Range(0, n_images), cv::Range(0, 1)).copyTo(real1);
							eigenvector.im(cv::Range(0, n_images), cv::Range(0, 1)).copyTo(imag1);
							cv::transpose(real1, real2); cv::transpose(imag1, imag2);
							real = real1 * real2 + imag1 * imag2; real = real * value1;
							imag = imag1 * real2 - real1 * imag2; imag = imag * value1;*/

							/*eigenvector.re(cv::Range(0, n_images), cv::Range(1, 2)).copyTo(real1);
							eigenvector.im(cv::Range(0, n_images), cv::Range(1, 2)).copyTo(imag1);
							cv::transpose(real1, real2); cv::transpose(imag1, imag2);
							real += (real1 * real2 + imag1 * imag2) * value2;
							imag += (imag1 * real2 - real1 * imag2) * value2;*/
							/*coherence_matrix.SetRe(real);
							coherence_matrix.SetIm(imag);
							real = coherence_matrix.GetPhase();
							phase.at<double>(ii, jj) = real.at<double>(0, 2);*/

							//cout << eigenvalue << endl;
							if (eigenvalue.at<double>(1, 0) / (eigenvalue.at<double>(0, 0) + 1e-10) < thresh_c1_to_c2)
							{
								for (int kk = 0; kk < n_images; kk++)
								{
									slc_series_filter[kk].re.at<double>(ii, jj) = eigenvector.re.at<double>(kk, 0);
									slc_series_filter[kk].im.at<double>(ii, jj) = eigenvector.im.at<double>(kk, 0);
								}
							}
							//else
							//{
							//	fprintf(stdout, "not processed!\n");
							//}
						}
					}


				}
			}
			//util.savephase("G:\\tmp\\phase_test2.jpg", "jet", phase); return 0;
			//储存
			for (int kk = 0; kk < n_images; kk++)
			{
				slc_series_filter[kk].re(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ph.convertTo(ph, CV_32F);
				ret = conversion.write_subarray_to_h5(coregis_slc_files_out[kk].c_str(), "s_re", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;

				slc_series_filter[kk].im(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ph.convertTo(ph, CV_32F);
				ret = conversion.write_subarray_to_h5(coregis_slc_files_out[kk].c_str(), "s_im", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;

			}
			slc_series.clear();
			slc_series_filter.clear();
			printf("\r估计进度：%lf %", double(i * block_num_col + j + 1) / double((block_num_col) * (block_num_row)) * 100.0);
			fflush(stdout);
		}
	}

	return 0;
}

int SLC_simulator::MB_phase_estimation(
	int estimation_wndsize,
	const char* slcH5File1,
	const char* slcH5File2,
	const char* slcH5File3,
	const char* slcH5File4,
	const char* slcH5File5,
	const char* slcH5File6,
	const char* slcH5File7, 
	const char* slcH5File8,
	const char* slcH5File1_out,
	const char* slcH5File2_out,
	const char* slcH5File3_out,
	const char* slcH5File4_out,
	const char* slcH5File5_out,
	const char* slcH5File6_out,
	const char* slcH5File7_out,
	const char* slcH5File8_out
)
{
	if (estimation_wndsize % 2 == 0 || estimation_wndsize < 5 ||
		!slcH5File1 || !slcH5File1_out || !slcH5File2 ||
		!slcH5File2_out || !slcH5File3 || !slcH5File3_out ||
		!slcH5File4 || !slcH5File4_out || !slcH5File5 || !slcH5File5_out || !slcH5File6 ||
		!slcH5File6_out || !slcH5File7 || !slcH5File7_out ||
		!slcH5File8 || !slcH5File8_out
		)
	{
		fprintf(stderr, "MB_phase_estimation(): input check failed!\n");
		return -1;
	}
	vector<string>coregis_slc_files; vector<string>coregis_slc_files_out;
	coregis_slc_files.push_back(slcH5File1);
	coregis_slc_files.push_back(slcH5File2);
	coregis_slc_files.push_back(slcH5File3);
	coregis_slc_files.push_back(slcH5File4);
	coregis_slc_files.push_back(slcH5File5);
	coregis_slc_files.push_back(slcH5File6);
	coregis_slc_files.push_back(slcH5File7);
	coregis_slc_files.push_back(slcH5File8);
	coregis_slc_files_out.push_back(slcH5File1_out);
	coregis_slc_files_out.push_back(slcH5File2_out);
	coregis_slc_files_out.push_back(slcH5File3_out);
	coregis_slc_files_out.push_back(slcH5File4_out);
	coregis_slc_files_out.push_back(slcH5File5_out);
	coregis_slc_files_out.push_back(slcH5File6_out);
	coregis_slc_files_out.push_back(slcH5File7_out);
	coregis_slc_files_out.push_back(slcH5File8_out);

	FormatConversion conversion; Utils util;
	double thresh_c1_to_c2 = 0.9;
	int master_indx = 1;
	int n_images = 8, nr, nc, blocksize_row = 1000, blocksize_col = 1000;
	int homotest_radius = (estimation_wndsize - 1) / 2;
	int left, right, top, bottom, block_num_row, block_num_col, left_pad, right_pad, top_pad, bottom_pad, ret;
	vector<ComplexMat> slc_series, slc_series_filter;
	ComplexMat slc;
	Mat ph, phase;

	////预先填充
	ret = conversion.read_slc_from_h5(slcH5File1, slc);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	nr = slc.GetRows(); nc = slc.GetCols();
	//ret = conversion.creat_new_h5(slcH5File1_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File2_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File3_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File4_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File5_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File6_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File7_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.creat_new_h5(slcH5File8_out);
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;


	//ret = conversion.write_slc_to_h5(slcH5File1_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File2_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File3_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File4_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File5_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File6_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File7_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;
	//ret = conversion.write_slc_to_h5(slcH5File8_out, slc);
	//if (return_check(ret, "write_slc_to_h5()", error_head)) return -1;

	//phase.create(nr, nc, CV_64F); phase = 0.0;
	//ret = conversion.creat_new_h5("G:\\tmp\\dual_phase.h5");
	//if (return_check(ret, "creat_new_h5()", error_head)) return -1;
	//ret = conversion.write_array_to_h5("G:\\tmp\\dual_phase.h5", "phase", phase);
	//if (return_check(ret, "write_array_to_h5()", error_head)) return -1;

	if (nr % blocksize_row == 0) block_num_row = nr / blocksize_row;
	else block_num_row = int(floor((double)nr / (double)blocksize_row)) + 1;
	if (nc % blocksize_col == 0) block_num_col = nc / blocksize_col;
	else block_num_col = int(floor((double)nc / (double)blocksize_col)) + 1;
	for (int i = 0; i < block_num_row; i++)
	{
		for (int j = 0; j < block_num_col; j++)
		{
			top = i * blocksize_row;
			top_pad = top - homotest_radius; top_pad = top_pad < 0 ? 0 : top_pad;
			bottom = top + blocksize_row; bottom = bottom > nr ? nr : bottom;
			bottom_pad = bottom + homotest_radius; bottom_pad = bottom_pad > nr ? nr : bottom_pad;
			left = j * blocksize_col;
			left_pad = left - homotest_radius; left_pad = left_pad < 0 ? 0 : left_pad;
			right = left + blocksize_col; right = right > nc ? nc : right;
			right_pad = right + homotest_radius; right_pad = right_pad > nc ? nc : right_pad;

			//读取数据
			for (int k = 0; k < n_images; k++)
			{
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_re",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.re);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				ret = conversion.read_subarray_from_h5(coregis_slc_files[k].c_str(), "s_im",
					top_pad, left_pad, bottom_pad - top_pad, right_pad - left_pad, slc.im);
				if (return_check(ret, "read_subarray_from_h5()", error_head)) return -1;
				if (slc.type() != CV_64F) slc.convertTo(slc, CV_64F);
				slc_series.push_back(slc);
				//slc_series_filter.push_back(slc);
			}
			phase.create(slc.GetRows(), slc.GetCols(), CV_64F); phase = 0.0;
			//计算
#pragma omp parallel for schedule(guided)
			for (int ii = (top - top_pad); ii < (bottom - top_pad); ii++)
			{
				ComplexMat coherence_matrix, eigenvector; Mat eigenvalue; int ret1;
				Mat real1, imag1, real2, imag2, real, imag; 
				double value1, value2;
				double re, im, re1, im1, re2, im2;
				for (int jj = (left - left_pad); jj < (right - left_pad); jj++)
				{
					ret1 = util.coherence_matrix_estimation(slc_series, coherence_matrix, estimation_wndsize, estimation_wndsize, ii, jj, false, true);
					if (ret1 == 0)
					{
						ret1 = util.HermitianEVD(coherence_matrix, eigenvalue, eigenvector);
						if (!eigenvalue.empty() && ret1 == 0)
						{
							//cout << coherence_matrix.GetMod() << endl;
							//eigenvector * eigenvector;
							value1 = eigenvalue.at<double>(0, 0);
							value2 = eigenvalue.at<double>(1, 0);
							re1 = eigenvector.re.at<double>(0, 0); im1 = eigenvector.im.at<double>(0, 0);
							re2 = eigenvector.re.at<double>(2, 0); im2 = eigenvector.im.at<double>(2, 0);
							re = re1 * re2 + im1 * im2; re *= value1;
							im = im1 * re2 - re1 * im2; im *= value1;

							re1 = eigenvector.re.at<double>(0, 1); im1 = eigenvector.im.at<double>(0, 1);
							re2 = eigenvector.re.at<double>(2, 1); im2 = eigenvector.im.at<double>(2, 1);
							re += (re1 * re2 + im1 * im2) * value2;
							im += (im1 * re2 - re1 * im2) * value2;

							/*eigenvector.re(cv::Range(0, n_images), cv::Range(0, 1)).copyTo(real1);
							eigenvector.im(cv::Range(0, n_images), cv::Range(0, 1)).copyTo(imag1);
							cv::transpose(real1, real2); cv::transpose(imag1, imag2);
							real = real1 * real2 + imag1 * imag2; real = real * value1;
							imag = imag1 * real2 - real1 * imag2; imag = imag * value1;*/

							/*eigenvector.re(cv::Range(0, n_images), cv::Range(1, 2)).copyTo(real1);
							eigenvector.im(cv::Range(0, n_images), cv::Range(1, 2)).copyTo(imag1);
							cv::transpose(real1, real2); cv::transpose(imag1, imag2);
							real += (real1 * real2 + imag1 * imag2) * value2;
							imag += (imag1 * real2 - real1 * imag2) * value2;
							coherence_matrix.SetRe(real);
							coherence_matrix.SetIm(imag);
							real = coherence_matrix.GetPhase();*/
							phase.at<double>(ii, jj) = atan2(im, re);
							/*if (eigenvalue.at<double>(1, 0) / (eigenvalue.at<double>(0, 0) + 1e-10) < thresh_c1_to_c2)
							{
								for (int kk = 0; kk < n_images; kk++)
								{
									slc_series_filter[kk].re.at<double>(ii, jj) = eigenvector.re.at<double>(kk, 0);
									slc_series_filter[kk].im.at<double>(ii, jj) = eigenvector.im.at<double>(kk, 0);
								}
							}
							else
							{
								fprintf(stdout, "not processed!\n");
							}*/
						}
					}


				}
			}
			//util.cvmat2bin("G:\\tmp\\phase_test.bin", phase); return 0;
			////储存
			phase(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
			ret = conversion.write_subarray_to_h5("G:\\tmp\\dual_phase.h5", "phase", ph, top, left, bottom - top, right - left);
			if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;

			/*for (int kk = 0; kk < n_images; kk++)
			{
				slc_series_filter[kk].re(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ph.convertTo(ph, CV_32F);
				ret = conversion.write_subarray_to_h5(coregis_slc_files_out[kk].c_str(), "s_re", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;

				slc_series_filter[kk].im(cv::Range(top - top_pad, bottom - top_pad), cv::Range(left - left_pad, right - left_pad)).copyTo(ph);
				ph.convertTo(ph, CV_32F);
				ret = conversion.write_subarray_to_h5(coregis_slc_files_out[kk].c_str(), "s_im", ph, top, left, bottom - top, right - left);
				if (return_check(ret, "write_subarray_to_h5()", error_head)) return -1;
			}*/
			slc_series.clear();
			//slc_series_filter.clear();
			printf("\r估计进度：%lf %", double(i * block_num_col + j + 1) / double((block_num_col) * (block_num_row)) * 100.0);
			fflush(stdout);
		}
	}
	return 0;
}

int SLC_simulator::pingpong_MLE(
	Mat& phase_reference,
	Mat& wrapped_phase_low,
	Mat& wrapped_phase_high,
	Mat& outphase, 
	Mat& lat_coef,
	Mat& lon_coef,
	double wavelength_low,
	double wavelength_high,
	double nearRangeTime, 
	double rangeSpacing,
	double azimuthSpacing,
	int offset_row, 
	int offset_col,
	double start1,
	double end1, 
	double start2,
	double end2,
	Mat& statevec1,
	Mat& statevec2,
	double prf,
	string demPath
)
{
	if (phase_reference.size() != wrapped_phase_high.size() ||
		phase_reference.size() != wrapped_phase_low.size() ||
		phase_reference.type() != CV_64F ||
		phase_reference.type() != wrapped_phase_high.type() ||
		phase_reference.type() != wrapped_phase_low.type() ||
		wavelength_high <= 0.0 ||
		wavelength_low <= 0.0 || nearRangeTime <= 0.0 ||
		rangeSpacing <= 0.0 || azimuthSpacing <= 0.0 || offset_row < 0 || offset_col < 0 || start1 <= 0.0 ||
		start2 <= 0.0 || end2 <= start2 || end1 <= start1 ||
		statevec1.type() != CV_64F ||
		statevec2.type() != CV_64F ||
		statevec1.cols != 7 || statevec2.cols != 7 ||
		prf <= 0.0 ||
		lon_coef.rows != 1 ||
		lon_coef.cols != 32 ||
		lon_coef.type() != CV_64F ||
		lat_coef.rows != 1 ||
		lat_coef.cols != 32 ||
		lat_coef.type() != CV_64F || demPath.empty()
		)
	{
		fprintf(stderr, "pingpong_MLE(): input check failed!\n");
		return -1;
	}
	Deflat flat; Utils util;
	//去参考平面
	double lonMax, lonMin, latMax, latMin, lon_upperleft, lat_upperleft;
	Mat dem, mappedDEM, mappedLat, mappedLon;
	int ret;
	int sceneHeight = wrapped_phase_low.rows;
	int sceneWidth = wrapped_phase_low.cols;
	ret = Utils::computeImageGeoBoundry(lat_coef, lon_coef, sceneHeight, sceneWidth, offset_row, offset_col,
		&lonMax, &latMax, &lonMin, &latMin);
	ret = Utils::getSRTMDEM(demPath.c_str(), dem, &lon_upperleft, &lat_upperleft, lonMin, lonMax, latMin, latMax);
	ret = flat.demMapping(dem, mappedDEM, mappedLat, mappedLon, lon_upperleft, lat_upperleft, offset_row, offset_col, sceneHeight, sceneWidth,
		prf, rangeSpacing, wavelength_high, nearRangeTime, start1, end1, statevec1, 20);
	mappedDEM = 0;
	//计算卫星位置
	Mat sate1 = Mat::zeros(sceneHeight, 3, CV_64F);

	Mat sate2 = Mat::zeros(sceneHeight, 3, CV_64F);

	orbitStateVectors stateVectors(statevec1, start1, end1);
	stateVectors.applyOrbit();

	orbitStateVectors stateVectors2(statevec2, start2, end2);
	stateVectors2.applyOrbit();

	//计算主星成像位置
	double dopplerFrequency = 0.0;
	Position groundPosition;
	double lat, lon, height;
	lat = mappedLat.at<float>(0, 0);
	lon = mappedLon.at<float>(0, 0);
	lon = lon > 180.0 ? (lon - 360.0) : lon;
	height = mappedDEM.at<short>(0, 0);
	Utils::ell2xyz(lon, lat, height, groundPosition);
	int numOrbitVec = stateVectors.newStateVectors.rows;
	double firstVecTime = 0.0;
	double secondVecTime = 0.0;
	double firstVecFreq = 0.0;
	double secondVecFreq = 0.0;
	double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			stateVectors.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			stateVectors.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength_high * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	double lowerBoundTime = firstVecTime;
	double upperBoundTime = secondVecTime;
	double lowerBoundFreq = firstVecFreq;
	double upperBoundFreq = secondVecFreq;
	double midTime, midFreq;
	double diffTime = fabs(upperBoundTime - lowerBoundTime);
	double absLineTimeInterval = 1.0 / prf;

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
		midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength_high * distance);
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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors.getPosition(time, pos);
		sate1.at<double>(i, 0) = pos.x;
		sate1.at<double>(i, 1) = pos.y;
		sate1.at<double>(i, 2) = pos.z;
	}

	//计算辅星成像位置
	numOrbitVec = stateVectors2.newStateVectors.rows;
	firstVecTime = 0.0;
	secondVecTime = 0.0;
	firstVecFreq = 0.0;
	secondVecFreq = 0.0;
	currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
	for (int ii = 0; ii < numOrbitVec; ii++) {
		Position orb_pos(stateVectors2.newStateVectors.at<double>(ii, 1), stateVectors2.newStateVectors.at<double>(ii, 2),
			stateVectors2.newStateVectors.at<double>(ii, 3));
		Velocity orb_vel(stateVectors2.newStateVectors.at<double>(ii, 4), stateVectors2.newStateVectors.at<double>(ii, 5),
			stateVectors2.newStateVectors.at<double>(ii, 6));
		currentFreq = 0;
		xdiff = groundPosition.x - orb_pos.x;
		ydiff = groundPosition.y - orb_pos.y;
		zdiff = groundPosition.z - orb_pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength_high * distance);
		if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			firstVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			firstVecFreq = currentFreq;
		}
		else {
			secondVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			secondVecFreq = currentFreq;
			break;
		}
	}

	if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
		fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
		return -1;
	}

	lowerBoundTime = firstVecTime;
	upperBoundTime = secondVecTime;
	lowerBoundFreq = firstVecFreq;
	upperBoundFreq = secondVecFreq;
	diffTime = fabs(upperBoundTime - lowerBoundTime);
	absLineTimeInterval = 1.0 / prf;

	totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
	numIterations = 0;
	while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

		midTime = (upperBoundTime + lowerBoundTime) / 2.0;
		stateVectors2.getPosition(midTime, pos);
		stateVectors2.getVelocity(midTime, vel);
		xdiff = groundPosition.x - pos.x;
		ydiff = groundPosition.y - pos.y;
		zdiff = groundPosition.z - pos.z;
		distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
		midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength_high * distance);
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

	for (int i = 0; i < sceneHeight; i++)
	{
		double time = zeroDopplerTime + (double)i * (1.0 / prf);
		stateVectors2.getPosition(time, pos);
		sate2.at<double>(i, 0) = pos.x;
		sate2.at<double>(i, 1) = pos.y;
		sate2.at<double>(i, 2) = pos.z;
	}

	Mat R1(sceneHeight, sceneWidth, CV_64F), R2(sceneHeight, sceneWidth, CV_64F);
	Mat phase_reference_deflat(sceneHeight, sceneWidth, CV_64F);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double r1, r2;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = mappedLat.at<float>(i, j);
			LLH.at<double>(0, 1) = mappedLon.at<float>(i, j);
			LLH.at<double>(0, 2) = mappedDEM.at<short>(i, j);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1(cv::Range(i, i + 1), cv::Range(0, 3));
			R1.at<double>(i, j) = cv::norm(tt, cv::NORM_L2);
			r1 = R1.at<double>(i, j);
			tt = XYZ - sate2(cv::Range(i, i + 1), cv::Range(0, 3));
			R2.at<double>(i, j) = cv::norm(tt, cv::NORM_L2);
			r2 = R2.at<double>(i, j);
			phase_reference_deflat.at<double>(i, j) = phase_reference.at<double>(i, j) + 4 * PI * (r1 - r2) / wavelength_high;
			wrapped_phase_high.at<double>(i, j) = wrapped_phase_high.at<double>(i, j) + 4 * PI * (r1 - r2) / wavelength_high;
			wrapped_phase_low.at<double>(i, j) = wrapped_phase_low.at<double>(i, j) + 4 * PI * (r1 - r2) / wavelength_low;
		}
	}
	util.wrap(wrapped_phase_high, wrapped_phase_high);
	util.wrap(wrapped_phase_low, wrapped_phase_low);
	wrapped_phase_high.copyTo(outphase);
	Mat coherence;
	util.phase_coherence(wrapped_phase_high, 7, 7, coherence);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < sceneHeight; i++)
	{
		for (int j = 0; j < sceneWidth; j++)
		{
			double gamma = coherence.at<double>(i, j)/*0.98*/;
			double gamma_square = gamma * gamma;
			double sigma = sqrt(1.0 - gamma * gamma) / (gamma * sqrt(2.0 * 1.0));
			int n_bin = 10000;
			double delta = 6 * sigma / n_bin;
			double phi_high, phi_low, ref;
			phi_high = wrapped_phase_high.at<double>(i, j);
			phi_low = wrapped_phase_low.at<double>(i, j);
			ref = phase_reference_deflat.at<double>(i, j);
			Mat psd(n_bin, 1, CV_64F);
			 
			for (int k = 0; k < n_bin; k++)
			{
				double reference = ref - 3 * sigma + k * delta;
				double cos_high = cos(phi_high - reference);
				double cos_low = cos(phi_low - reference * wavelength_high / wavelength_low);
				double f_high = 1.0 / (2.0 * PI) * (1.0 - gamma_square) / (1.0 - gamma_square * cos_high * cos_high) *
					(1.0 + gamma * cos_high * acos(-1.0 * gamma * cos_high) / sqrt(1.0 - gamma_square * cos_high * cos_high));
				double f_low = 1.0 / (2.0 * PI) * (1.0 - gamma_square) / (1.0 - gamma_square * cos_low * cos_low) *
					(1.0 + gamma * cos_low * acos(-1.0 * gamma * cos_low) / sqrt(1.0 - gamma_square * cos_low * cos_low));
				psd.at<double>(k, 0) = f_low * f_high;
			}
			cv::Point peak;
			cv::minMaxLoc(psd, NULL, NULL, NULL, &peak);
			//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\psd.bin", psd); return 0;
			outphase.at<double>(i, j) = ((peak.x + 1) * (peak.y + 1) - 1) * delta + ref - 3 * sigma;
		}
	}
	outphase = outphase - 4 * PI * (R1 - R2) / wavelength_high;
	return 0;
}


