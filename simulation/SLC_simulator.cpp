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
	//Mat pixel_count(sceneHeight, sceneWidth, CV_8U); pixel_count = 0;
	//Mat R(sceneHeight, sceneWidth, CV_32F); R = 0.0;
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
	num_block_row = 2;
	num_block_col = num_block_col;
	vector<double> GCPs;//控制点信息
	//初始化轨道类
	orbitStateVectors stateVectors(stateVec, acquisitionStartTime, acquisitionStopTime);
	stateVectors.applyOrbit();
	int ret;
	double time_interval = 1.0 / prf;
	double dopplerFrequency = 0.0;
	uint64 seed = 0;
	for (int i = 15; i < 17; i++)
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

			for (int ii = 0; ii < DEM_rows; ii++)
			{
				for (int jj = 0; jj < DEM_cols; jj++)
				{
					double zeroDopplerTime = imaging_time.at<double>(ii, jj);
					double distance = slant_range.at<double>(ii, jj);
					double real, theta, imaginary;
					int azimuthIndex = floor((zeroDopplerTime - acquisitionStartTime) / time_interval);
					int rangeIndex = floor((distance - nearRange) / rangeSpacing);
					if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
					{

					}
					else
					{
						//pixel_count.at<uchar>(azimuthIndex, rangeIndex) += 1;
						//R.at<float>(azimuthIndex, rangeIndex) += distance;
						theta = -4.0 * PI * distance / wavelength + randomAngle.at<float>(ii, jj);
						real = 1.0 * (cos(theta) + noise_real.at<float>(ii, jj));
						imaginary = 1.0 * (sin(theta) + noise_imaginary.at<float>(ii, jj));
						slc.re.at<float>(azimuthIndex, rangeIndex) += real;
						slc.im.at<float>(azimuthIndex, rangeIndex) += imaginary;
						
					}
				}
			}
			
			//控制点信息
			if ((i % 5 == 0) && (j % 5 == 0) && i != 0 && j != 0)
			{
				int gcp_row = DEM_rows / 2;
				int gcp_col = DEM_cols / 2;
				double gcp_sigma = 1000.0;
				double zeroDopplerTime = imaging_time.at<double>(gcp_row, gcp_col);
				double distance = slant_range.at<double>(gcp_row, gcp_col);
				int azimuthIndex = round((zeroDopplerTime - acquisitionStartTime) / time_interval);
				int rangeIndex = round((distance - nearRange) / rangeSpacing);
				if (azimuthIndex < 0 || azimuthIndex > sceneHeight - 1 || rangeIndex < 0 || rangeIndex > sceneWidth - 1)
				{

				}
				else
				{
					double lat, lon, height;
					lat = upper_left_lat - (double)gcp_row * lat_spacing;
					lon = upper_left_lon + (double)gcp_col * lon_spacing;
					lon = lon > 180.0 ? (lon - 360.0) : lon;
					height = dem_temp2.at<float>(gcp_row, gcp_col);
					GCPs.push_back(azimuthIndex + 1);
					GCPs.push_back(rangeIndex + 1);
					GCPs.push_back(lon);
					GCPs.push_back(lat);
					GCPs.push_back(height);
					GCPs.push_back(distance);
					GCPs.push_back(distance);
					double theta = -4.0 * PI * distance / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					double real = gcp_sigma * (cos(theta) + noise_real.at<float>(gcp_row, gcp_col));
					double imaginary = gcp_sigma * (sin(theta) + noise_imaginary.at<float>(gcp_row, gcp_col));
					slc.re.at<float>(azimuthIndex, rangeIndex) += real;
					slc.im.at<float>(azimuthIndex, rangeIndex) += imaginary;
				}
			}
			
			fprintf(stdout, "process %lf: %d / %d\n", (double)seed / (double)(num_block_row * num_block_col) * 100.0, seed,
				num_block_row * num_block_col);
		}
	}
	//pixel_count.convertTo(pixel_count, CV_64F);
	//util.cvmat2bin("G:\\tmp\\pixel_count.bin", pixel_count);

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
	double res_range = rangeSpacing * 2.0;
	double res_azimuth = azimuthSpacing * 2.0;

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
	Mat& GCP
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
			cv::RNG rng2(seed + 1);
			rng2.fill(noise_real, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng3(seed + 2);
			rng3.fill(noise_imaginary, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng4(seed + 3);
			rng4.fill(noise_real2, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng5(seed + 4);
			rng5.fill(noise_imaginary2, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng6(seed + 5);
			rng6.fill(noise_real3, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng7(seed + 6);
			rng7.fill(noise_imaginary3, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng8(seed + 7);
			rng8.fill(noise_real4, cv::RNG::NORMAL, 0, noise_sigma);
			cv::RNG rng9(seed + 8);
			rng9.fill(noise_imaginary4, cv::RNG::NORMAL, 0, noise_sigma);

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
					int azimuthIndex1 = round((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
					int rangeIndex1 = round((distance1 - nearRange1) / rangeSpacing);

					double zeroDopplerTime2 = imaging_time2.at<double>(ii, jj);
					double distance2 = slant_range2.at<double>(ii, jj);
					int azimuthIndex2 = round((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
					int rangeIndex2 = round((distance2 - nearRange2) / rangeSpacing);

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


						theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(ii, jj);
						real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real4.at<float>(ii, jj));
						imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary4.at<float>(ii, jj));
						slc4.re.at<float>(azimuthIndex1, rangeIndex1) += real;
						slc4.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;
					}

					
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


						theta = -2.0 * PI * (distance2 + distance2) / wavelength + randomAngle.at<float>(ii, jj);
						real = /*sigma.at<double>(ii, jj)*/ 1.0 * (cos(theta) + noise_real3.at<float>(ii, jj));
						imaginary = /*sigma.at<double>(ii, jj)*/ 1.0 * (sin(theta) + noise_imaginary3.at<float>(ii, jj));
						slc3.re.at<float>(azimuthIndex2, rangeIndex2) += real;
						slc3.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;

						
					}

				}
			}
			//控制点信息
			if ((i % 5 == 0) && (j % 5 == 0) && i != 0 && j != 0)
			{
				int gcp_row = DEM_rows / 2;
				int gcp_col = DEM_cols / 2;
				double gcp_sigma = 1000.0;
				double zeroDopplerTime1 = imaging_time1.at<double>(gcp_row, gcp_col);
				double distance1 = slant_range1.at<double>(gcp_row, gcp_col);
				int azimuthIndex1 = round((zeroDopplerTime1 - acquisitionStartTime1) / time_interval);
				int rangeIndex1 = round((distance1 - nearRange1) / rangeSpacing);

				double zeroDopplerTime2 = imaging_time2.at<double>(gcp_row, gcp_col);
				double distance2 = slant_range2.at<double>(gcp_row, gcp_col);
				int azimuthIndex2 = round((zeroDopplerTime2 - acquisitionStartTime2) / time_interval);
				int rangeIndex2 = round((distance2 - nearRange2) / rangeSpacing);

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
					slc2.re.at<float>(azimuthIndex1, rangeIndex1) += real;
					slc2.im.at<float>(azimuthIndex1, rangeIndex1) += imaginary;

					theta = -2.0 * PI * (distance2 + distance2) / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					real = gcp_sigma * (cos(theta) + noise_real3.at<float>(gcp_row, gcp_col));
					imaginary = gcp_sigma * (sin(theta) + noise_imaginary3.at<float>(gcp_row, gcp_col));
					slc3.re.at<float>(azimuthIndex2, rangeIndex2) += real;
					slc3.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;

					theta = -2.0 * PI * (distance1 + distance2) / wavelength + randomAngle.at<float>(gcp_row, gcp_col);
					real = gcp_sigma * (cos(theta) + noise_real4.at<float>(gcp_row, gcp_col));
					imaginary = gcp_sigma * (sin(theta) + noise_imaginary4.at<float>(gcp_row, gcp_col));
					slc4.re.at<float>(azimuthIndex2, rangeIndex2) += real;
					slc4.im.at<float>(azimuthIndex2, rangeIndex2) += imaginary;
				}
			}
			//memset(process, 0, 512);
			printf("\rprocess %lf: %d / %d", (double)seed / (double)(num_block_row * num_block_col) * 100.0, seed,
				num_block_row* num_block_col);
			fflush(stdout);
			//std::cout << "\r" << process;
			//std::cout.flush();
		}
	}

	////二维复数卷积
	int wright = 16;
	double Drange = rangeSpacing;
	double Dazimuth = azimuthSpacing;
	double res_range = rangeSpacing * 2.0;
	double res_azimuth = azimuthSpacing * 2.0;

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


