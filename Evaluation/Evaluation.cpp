#include<Evaluation.h>
#include<FormatConversion.h>
#include<Utils.h>
#include<ComplexMat.h>
#include<Deflat.h>

#ifdef _DEBUG
#pragma comment(lib, "FormatConversion_d.lib")
#pragma comment(lib, "Utils_d.lib")
#pragma comment(lib, "ComplexMat_d.lib")
#pragma comment(lib, "Deflat_d.lib")
#else
#pragma comment(lib, "FormatConversion.lib")
#pragma comment(lib, "Utils.lib")
#pragma comment(lib, "ComplexMat.lib")
#pragma comment(lib, "Deflat.lib")
#endif

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
inline bool parallel_check(volatile bool parallel_flag, const char* detail_info, const char* parallel_error_head)
{
	if (!parallel_flag)
	{
		fprintf(stderr, "%s %s\n\n", parallel_error_head, detail_info);
		return true;
	}
	else
	{
		return false;
	}
}

Evaluation::Evaluation()
{
	memset(this->error_head, 0, 256);
	memset(this->parallel_error_head, 0, 256);
	strcpy(this->error_head, "Evaluation_DLL_ERROR: error happens when using ");
	strcpy(this->parallel_error_head, "Evaluation_DLL_ERROR: error happens when using parallel computing in function: ");
}

Evaluation::~Evaluation()
{

}

int Evaluation::PhasePreserve(const char* master_h5,
	const char* slave_h5,
	double* Output)
{
	if (master_h5 == NULL ||
		slave_h5 == NULL ||
		Output == NULL)
	{
		fprintf(stderr, "PhasePreserve(): input check failed!\n\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;
	
	int ret = 0;
	double PhaseError = 0;
	ComplexMat master, slave;
	//读取SLC数据
	ret = conversion.read_slc_from_h5(master_h5, master);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (master.type() != CV_64F) master.convertTo(master, CV_64F);
	ret = conversion.read_slc_from_h5(slave_h5, slave);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
	int rows = master.GetRows();
	int cols = master.GetCols();
	Mat GCPS, slave_Gcps;
	//读取标志点数据
	ret = conversion.read_array_from_h5(master_h5, "GCP", GCPS);
	if (return_check(ret, "read_Gcps_from_h5()", error_head)) return -1;
	int Gcps_number = GCPS.rows;
	bool RealPhaseIsExisted = false;
	double lonMax, lonMin, latMax, latMin, wavelength, wavelength2, prf, prf2, start, start2, end, end2;
	int sceneHeight = rows, sceneWidth = cols;
	Mat lon_coef, lat_coef, statevec, statevec2;
	string start_time, start_time2, end_time, end_time2;
	//if (GCPS.cols == 6) RealPhaseIsExisted = true;
		//主星参数
		ret = conversion.read_double_from_h5(master_h5, "prf", &prf);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		ret = conversion.read_double_from_h5(master_h5, "carrier_frequency", &wavelength);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		wavelength = VEL_C / wavelength;
		ret = conversion.read_str_from_h5(master_h5, "acquisition_start_time", start_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time.c_str(), &start);
		ret = conversion.read_str_from_h5(master_h5, "acquisition_stop_time", end_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		conversion.utc2gps(end_time.c_str(), &end);
		ret = conversion.read_array_from_h5(master_h5, "state_vec", statevec);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		//辅星参数
		ret = conversion.read_double_from_h5(slave_h5, "prf", &prf2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		ret = conversion.read_double_from_h5(slave_h5, "carrier_frequency", &wavelength2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		wavelength2 = VEL_C / wavelength2;
		ret = conversion.read_str_from_h5(slave_h5, "acquisition_start_time", start_time2);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time2.c_str(), &start2);
		ret = conversion.read_str_from_h5(slave_h5, "acquisition_stop_time", end_time2);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		conversion.utc2gps(end_time2.c_str(), &end2);
		ret = conversion.read_array_from_h5(slave_h5, "state_vec", statevec2);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		
	int interp_times = 32;
	int win_size = 16;
	int count = 0; //符合要求的标志点个数

	Mat Inphase_pre = Mat::zeros(Size(Gcps_number, 1), CV_64FC1);
	Mat Inphase_pro = Mat::zeros(Size(Gcps_number, 1), CV_64FC1);
	Mat Error = Mat::zeros(Size(Gcps_number, 1), CV_64FC1);
	
	

	for (int i = 0; i < Gcps_number; i++)
	{
		/*后验相位-插值*/
		ComplexMat master_win;
		ComplexMat slave_win;
		ComplexMat master_interp, slave_interp;
		int row = GCPS.at<double>(i, 0);
		int col = GCPS.at<double>(i, 1);
		if (row >= win_size && row < rows - win_size && col >= win_size && col < cols - win_size)
		{
			master_win = master(Range(row - win_size, row + win_size), Range(col - win_size, col + win_size));
			slave_win = slave(Range(row - win_size, row + win_size), Range(col - win_size, col + win_size));
		}
		else
			continue;
		resize(master_win.re, master_interp.re, Size(2 * win_size * interp_times, 2 * win_size * interp_times), 0, 0, INTER_CUBIC);
		resize(master_win.im, master_interp.im, Size(2 * win_size * interp_times, 2 * win_size * interp_times), 0, 0, INTER_CUBIC);
		resize(slave_win.re, slave_interp.re, Size(2 * win_size * interp_times, 2 * win_size * interp_times), 0, 0, INTER_CUBIC);
		resize(slave_win.im, slave_interp.im, Size(2 * win_size * interp_times, 2 * win_size * interp_times), 0, 0, INTER_CUBIC);
		Mat master_mod = master_interp.GetMod();
		Mat slave_mod = slave_interp.GetMod();
		Point master_max, slave_max;
		minMaxLoc(master_mod, NULL, NULL, NULL, &master_max);
		minMaxLoc(slave_mod, NULL, NULL, NULL, &slave_max);
		double master_phase = atan2(master_interp.im.at<double>(master_max.y, master_max.x),
			master_interp.re.at<double>(master_max.y, master_max.x));
		double slave_phase = atan2(slave_interp.im.at<double>(slave_max.y, slave_max.x),
			slave_interp.re.at<double>(slave_max.y, slave_max.x));
		double Inphase_post = atan2(sin(master_phase - slave_phase), cos(master_phase - slave_phase));
		Inphase_pre.at<double>(i) = Inphase_post;
		/*先验相位-斜距差*/
		double Inphase_prior = 0;
			////主卫星斜距
			//Mat sate1 = Mat::zeros(1, 3, CV_64F);
			//orbitStateVectors stateVectors(statevec, start, end);
			//stateVectors.applyOrbit();

			//double dopplerFrequency = 0.0;

			//Position groundPosition;
			//double lat, lon, height;
			//lat = GCPS.at<double>(i, 3);
			//lon = GCPS.at<double>(i, 2);
			//lon = lon > 180.0 ? (lon - 360.0) : lon;
			//height = GCPS.at<double>(i, 4);
			//Utils::ell2xyz(lon, lat, height, groundPosition);
			//int numOrbitVec = stateVectors.newStateVectors.rows;
			//double firstVecTime = 0.0;
			//double secondVecTime = 0.0;
			//double firstVecFreq = 0.0;
			//double secondVecFreq = 0.0;
			//double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
			////检测标志点位于哪两个轨道点之间
			//for (int ii = 0; ii < numOrbitVec; ii++) {
			//	Position orb_pos(stateVectors.newStateVectors.at<double>(ii, 1), stateVectors.newStateVectors.at<double>(ii, 2),
			//		stateVectors.newStateVectors.at<double>(ii, 3));
			//	Velocity orb_vel(stateVectors.newStateVectors.at<double>(ii, 4), stateVectors.newStateVectors.at<double>(ii, 5),
			//		stateVectors.newStateVectors.at<double>(ii, 6));
			//	currentFreq = 0;
			//	xdiff = groundPosition.x - orb_pos.x;
			//	ydiff = groundPosition.y - orb_pos.y;
			//	zdiff = groundPosition.z - orb_pos.z;
			//	distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
			//	currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength * distance);
			//	if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			//		firstVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			//		firstVecFreq = currentFreq;
			//	}
			//	else {
			//		secondVecTime = stateVectors.newStateVectors.at<double>(ii, 0);
			//		secondVecFreq = currentFreq;
			//		break;
			//	}
			//}

			//if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
			//	fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
			//	return -1;
			//}

			//double lowerBoundTime = firstVecTime;
			//double upperBoundTime = secondVecTime;
			//double lowerBoundFreq = firstVecFreq;
			//double upperBoundFreq = secondVecFreq;
			//double midTime, midFreq;
			//double diffTime = fabs(upperBoundTime - lowerBoundTime);
			//double absLineTimeInterval = 1.0 / prf;

			//int totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
			//int numIterations = 0; Position pos; Velocity vel;
			////对两个点之间（相差10s）进行进一步插值检测找到标志点对应的具体卫星位置
			//while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

			//	midTime = (upperBoundTime + lowerBoundTime) / 2.0;
			//	stateVectors.getPosition(midTime, pos);
			//	stateVectors.getVelocity(midTime, vel);
			//	xdiff = groundPosition.x - pos.x;
			//	ydiff = groundPosition.y - pos.y;
			//	zdiff = groundPosition.z - pos.z;
			//	distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
			//	midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength * distance);
			//	if ((midFreq - dopplerFrequency) * (lowerBoundFreq - dopplerFrequency) > 0.0) {
			//		lowerBoundTime = midTime;
			//		lowerBoundFreq = midFreq;
			//	}
			//	else if ((midFreq - dopplerFrequency) * (upperBoundFreq - dopplerFrequency) > 0.0) {
			//		upperBoundTime = midTime;
			//		upperBoundFreq = midFreq;
			//	}
			//	else if (fabs(midFreq - dopplerFrequency) < 0.01) {
			//		zeroDopplerTime = midTime;
			//		break;
			//	}

			//	diffTime = fabs(upperBoundTime - lowerBoundTime);
			//	numIterations++;
			//}
			//zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
			//stateVectors.getPosition(zeroDopplerTime, pos);
			//sate1.at<double>(0) = pos.x;
			//sate1.at<double>(1) = pos.y;
			//sate1.at<double>(2) = pos.z;
			//double r;
			//Mat XYZ, LLH(1, 3, CV_64F), tt;
			//LLH.at<double>(0, 0) = GCPS.at<double>(i, 3);
			//LLH.at<double>(0, 1) = GCPS.at<double>(i, 2);
			//LLH.at<double>(0, 2) = GCPS.at<double>(i, 4);
			//util.ell2xyz(LLH, XYZ);
			//tt = XYZ - sate1;
			//r = cv::norm(tt, cv::NORM_L2);
			//master_phase = -r / wavelength * 4 * PI;

			////辅卫星斜距
			//Mat sate2 = Mat::zeros(1, 3, CV_64F);
			//orbitStateVectors stateVectors2(statevec2, start2, end2);
			//stateVectors2.applyOrbit();

			//dopplerFrequency = 0.0;

			//numOrbitVec = stateVectors2.newStateVectors.rows;
			//firstVecTime = 0.0;
			//secondVecTime = 0.0;
			//firstVecFreq = 0.0;
			//secondVecFreq = 0.0;
			//distance = 1.0;
			//for (int ii = 0; ii < numOrbitVec; ii++) {
			//	Position orb_pos(stateVectors2.newStateVectors.at<double>(ii, 1), stateVectors2.newStateVectors.at<double>(ii, 2),
			//		stateVectors2.newStateVectors.at<double>(ii, 3));
			//	Velocity orb_vel(stateVectors2.newStateVectors.at<double>(ii, 4), stateVectors2.newStateVectors.at<double>(ii, 5),
			//		stateVectors2.newStateVectors.at<double>(ii, 6));
			//	currentFreq = 0;
			//	xdiff = groundPosition.x - orb_pos.x;
			//	ydiff = groundPosition.y - orb_pos.y;
			//	zdiff = groundPosition.z - orb_pos.z;
			//	distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
			//	currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength2 * distance);
			//	if (ii == 0 || (firstVecFreq - dopplerFrequency) * (currentFreq - dopplerFrequency) > 0) {
			//		firstVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			//		firstVecFreq = currentFreq;
			//	}
			//	else {
			//		secondVecTime = stateVectors2.newStateVectors.at<double>(ii, 0);
			//		secondVecFreq = currentFreq;
			//		break;
			//	}
			//}

			//if ((firstVecFreq - dopplerFrequency) * (secondVecFreq - dopplerFrequency) >= 0.0) {
			//	fprintf(stderr, "SLC_deramp(): orbit mismatch!\n");
			//	return -1;
			//}

			//lowerBoundTime = firstVecTime;
			//upperBoundTime = secondVecTime;
			//lowerBoundFreq = firstVecFreq;
			//upperBoundFreq = secondVecFreq;

			//diffTime = fabs(upperBoundTime - lowerBoundTime);
			//absLineTimeInterval = 1.0 / prf2;

			//totalIterations = (int)(diffTime / absLineTimeInterval) + 1;
			//numIterations = 0;
			//while (diffTime > absLineTimeInterval * 0.1 && numIterations <= totalIterations) {

			//	midTime = (upperBoundTime + lowerBoundTime) / 2.0;
			//	stateVectors2.getPosition(midTime, pos);
			//	stateVectors2.getVelocity(midTime, vel);
			//	xdiff = groundPosition.x - pos.x;
			//	ydiff = groundPosition.y - pos.y;
			//	zdiff = groundPosition.z - pos.z;
			//	distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
			//	midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength2 * distance);
			//	if ((midFreq - dopplerFrequency) * (lowerBoundFreq - dopplerFrequency) > 0.0) {
			//		lowerBoundTime = midTime;
			//		lowerBoundFreq = midFreq;
			//	}
			//	else if ((midFreq - dopplerFrequency) * (upperBoundFreq - dopplerFrequency) > 0.0) {
			//		upperBoundTime = midTime;
			//		upperBoundFreq = midFreq;
			//	}
			//	else if (fabs(midFreq - dopplerFrequency) < 0.01) {
			//		zeroDopplerTime = midTime;
			//		break;
			//	}

			//	diffTime = fabs(upperBoundTime - lowerBoundTime);
			//	numIterations++;
			//}
			//zeroDopplerTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);

			//double r2;
			//stateVectors2.getPosition(zeroDopplerTime, pos);
			//sate2.at<double>(0) = pos.x;
			//sate2.at<double>(1) = pos.y;
			//sate2.at<double>(2) = pos.z;
			//LLH.at<double>(0, 0) = GCPS.at<double>(i, 3);
			//LLH.at<double>(0, 1) = GCPS.at<double>(i, 2);
			//LLH.at<double>(0, 2) = GCPS.at<double>(i, 4);
			//util.ell2xyz(LLH, XYZ);
			//tt = XYZ - sate2;
			//r2 = cv::norm(tt, cv::NORM_L2);
			//slave_phase = -r2 / wavelength2 * 4 * PI;

			//先验相位
		double r = GCPS.at<double>(i, 5);
		double r2 = GCPS.at<double>(i, 6);
		master_phase = -r * 4 * PI / wavelength;
		slave_phase = -r2 * 4 * PI / wavelength2;
		Inphase_prior = atan2(sin(master_phase - slave_phase), cos(master_phase - slave_phase));
		Inphase_pro.at<double>(i) = Inphase_prior;
		Error.at<double>(i) = atan2(sin(Inphase_post - Inphase_prior), cos(Inphase_post - Inphase_prior));
		PhaseError += Error.at<double>(i);
		count++;
	}
	double relevant_error = 0;
	double count2 = 0;
	for (int i = 0; i < Gcps_number; i++)
		for (int j = i + 1; j < Gcps_number; j++)
		{
				relevant_error += pow(atan2(sin(Inphase_pre.at<double>(i) - Inphase_pre.at<double>(j) - (Inphase_pro.at<double>(i) - Inphase_pro.at<double>(j))),
					cos(Inphase_pre.at<double>(i) - Inphase_pre.at<double>(j) - (Inphase_pro.at<double>(i) - Inphase_pro.at<double>(j)))),2);
				count2++;
		}
	relevant_error = sqrt(relevant_error / count2);
	util.cvmat2bin("D:\\Test\\Error.bin", Error);
	*Output = sqrt(PhaseError / count);
}

int Evaluation::Regis(const char* master_h5,  const char* slave_regis_h5, Mat& coherence, double* Output)
{
	if(master_h5 == NULL ||
		slave_regis_h5 == NULL ||
		Output == NULL)
	{
		fprintf(stderr, "PhasePreserve(): input check failed!\n\n");
		return -1;
	}
	FormatConversion conversion; Utils util;

	int ret = 0;
	double PhaseError = 0;
	ComplexMat master, slave;
	//读取SLC数据
	ret = conversion.read_slc_from_h5(master_h5, master);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (master.type() != CV_64F) master.convertTo(master, CV_64F);
	ret = conversion.read_slc_from_h5(slave_regis_h5, slave);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
	int rows = master.GetRows();
	int cols = master.GetCols();
	Mat GCPS, slave_Gcps;
	//读取标志点数据
	ret = conversion.read_array_from_h5(master_h5, "GCP", GCPS);
	if (return_check(ret, "read_Gcps_from_h5()", error_head)) return -1;
	int offset_row, offset_col;
	ret = conversion.read_int_from_h5(master_h5, "offset_row", &offset_row);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_int_from_h5(master_h5, "offset_col", &offset_col);

	int Gcps_number = GCPS.rows;
	int interp_times = 32;
	int win_size = 32;
	int interp_size = interp_times * win_size;
	int count = 0; //符合要求的标志点个数
	double regis_row_error = 0, regis_col_error = 0, regis_row_error2 = 0, regis_col_error2 = 0;
	FILE* fp3 = NULL;
	fp3 = fopen("D:\\Test_File\\Error.txt", "w+");
	fprintf(fp3, "误差\n");
	for (int i = 0; i < Gcps_number; i++)
	{
		ComplexMat master_win, master_fft;
		ComplexMat slave_win, slave_fft;
		ComplexMat master_interp(interp_size, interp_size),
			slave_interp(interp_size, interp_size),
			master_interp_fft(interp_size, interp_size),
			slave_interp_fft(interp_size, interp_size);
		int row = GCPS.at<double>(i, 0);
		int col = GCPS.at<double>(i, 1);
		if ((row - offset_row) >= win_size / 2 && (row - offset_row) < rows - win_size / 2 
			&& (col - offset_col) >= win_size / 2 && (col - offset_col) < cols - win_size / 2)
		{
			master_win = master(Range(row - offset_row - win_size/2, row - offset_row + win_size/2), Range(col - offset_col - win_size/2, col - offset_col + win_size/2));
			slave_win = slave(Range(row - offset_row - win_size / 2, row - offset_row + win_size / 2), Range(col - offset_col - win_size / 2, col - offset_col + win_size / 2));
		}
		else
			continue;
		//resize(master_win.re, master_interp.re, Size( win_size * interp_times,  win_size * interp_times), 0, 0, INTER_CUBIC);
		//resize(master_win.im, master_interp.im, Size( win_size * interp_times, win_size * interp_times), 0, 0, INTER_CUBIC);
		//resize(slave_win.re, slave_interp.re, Size( win_size * interp_times, win_size * interp_times), 0, 0, INTER_CUBIC);
		//resize(slave_win.im, slave_interp.im, Size( win_size * interp_times, win_size * interp_times), 0, 0, INTER_CUBIC);
		
		FFT2(master_win, master_interp, win_size, interp_times);
		FFT2(slave_win, slave_interp, win_size, interp_times);;
		Mat master_mod = master_interp.GetMod();
		Mat slave_mod = slave_interp.GetMod();
		Point master_max, slave_max;
		double max1, max2;
		minMaxLoc(master_mod, NULL, &max1, NULL, &master_max);
		minMaxLoc(slave_mod, NULL, &max2, NULL, &slave_max);
		util.cvmat2bin("D:/Test_File/master.bin", master_mod);
		util.cvmat2bin("D:/Test_File/slave.bin", slave_mod);
		util.saveSLC("D:/Test_File/master.bmp", 60, master_interp);
		util.saveSLC("D:/Test_File/slave.bmp", 60, slave_interp);
		regis_col_error += powf(double(master_max.x - slave_max.x) / interp_times, 2);
		regis_col_error2 += (double(master_max.x - slave_max.x) / interp_times);
		regis_row_error	+= powf(double(master_max.y - slave_max.y) / interp_times, 2);
		regis_row_error2 += (double(master_max.y - slave_max.y) / interp_times);
		count++;
	}
	
	regis_row_error = sqrt(regis_row_error / count);
	regis_col_error = sqrt(regis_col_error / count);
	fprintf(fp3, "%f\t", regis_row_error);
	fprintf(fp3, "%f\t", regis_col_error);
	fclose(fp3);
	*Output = regis_row_error;

	//util.real_coherence(master, slave, coherence);
	return 0;
}

int Evaluation::Unwrap(const char* master_h5, const char* slave_regis_h5, const char* phase_unwrapped_h5, double* Output)
{
	if (master_h5 == NULL ||
		slave_regis_h5 == NULL ||
		phase_unwrapped_h5 == NULL ||
		Output == NULL)
	{
		fprintf(stderr, "PhasePreserve(): input check failed!\n\n");
		return -1;
	}
	FormatConversion conversion; Deflat flat; Utils util;

	int ret = 0;
	double PhaseError = 0;
	ComplexMat master, slave;
	//读取SLC数据
	ret = conversion.read_slc_from_h5(master_h5, master);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (master.type() != CV_64F) master.convertTo(master, CV_64F);
	ret = conversion.read_slc_from_h5(slave_regis_h5, slave);
	if (return_check(ret, "read_slc_from_h5()", error_head)) return -1;
	if (slave.type() != CV_64F) slave.convertTo(slave, CV_64F);
	int rows = master.GetRows();
	int cols = master.GetCols();
	Mat GCPS, slave_Gcps;
	//读取标志点数据
	ret = conversion.read_array_from_h5(master_h5, "GCP", GCPS);
	if (return_check(ret, "read_Gcps_from_h5()", error_head)) return -1;
	int Gcps_number = GCPS.rows;
	Mat GCPS_New = Mat::zeros(Size(6, Gcps_number), CV_64FC1);
	bool RealPhaseIsExisted = false;
	double lonMax, lonMin, latMax, latMin, wavelength, wavelength2, prf, prf2, start, start2, end, end2;
	int sceneHeight = rows, sceneWidth = cols;
	Mat lon_coef, lat_coef, statevec, statevec2;
	string start_time, start_time2, end_time, end_time2;
	if (GCPS.cols == 6) RealPhaseIsExisted = true;
	else
	{
		GCPS.copyTo(GCPS_New(Range(0, Gcps_number), Range(0, 5)));
		//主星参数
		ret = conversion.read_double_from_h5(master_h5, "prf", &prf);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		ret = conversion.read_double_from_h5(master_h5, "carrier_frequency", &wavelength);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		wavelength = VEL_C / wavelength;
		ret = conversion.read_str_from_h5(master_h5, "acquisition_start_time", start_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time.c_str(), &start);
		ret = conversion.read_str_from_h5(master_h5, "acquisition_stop_time", end_time);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		conversion.utc2gps(end_time.c_str(), &end);
		ret = conversion.read_array_from_h5(master_h5, "state_vec", statevec);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
		//辅星参数
		ret = conversion.read_double_from_h5(slave_regis_h5, "prf", &prf2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		ret = conversion.read_double_from_h5(slave_regis_h5, "carrier_frequency", &wavelength2);
		if (return_check(ret, "read_double_from_h5()", error_head)) return -1;
		wavelength2 = VEL_C / wavelength2;
		ret = conversion.read_str_from_h5(slave_regis_h5, "acquisition_start_time", start_time2);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		ret = conversion.utc2gps(start_time2.c_str(), &start2);
		ret = conversion.read_str_from_h5(slave_regis_h5, "acquisition_stop_time", end_time2);
		if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
		conversion.utc2gps(end_time2.c_str(), &end2);
		ret = conversion.read_array_from_h5(slave_regis_h5, "state_vec", statevec2);
		if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	}

	for (int i = 0; i < Gcps_number; i++)
	{
		/*解缠相位*/
		int row = GCPS.at<double>(i, 0);
		int col = GCPS.at<double>(i, 1);
		Mat phase_unwrapped;
		ret = conversion.read_subarray_from_h5(phase_unwrapped_h5, "phase", row - 1, col - 1, 1, 1, phase_unwrapped);
		double Inphase_unwrapped = phase_unwrapped.at<double>(0, 0);
		/*真实相位*/
		double Inphase_prior = 0;
		if (RealPhaseIsExisted)
			Inphase_prior = GCPS.at<double>(i, 5);
		else
		{
			//主卫星斜距
			Mat sate1 = Mat::zeros(1, 3, CV_64F);
			orbitStateVectors stateVectors(statevec, start, end);
			stateVectors.applyOrbit();

			double dopplerFrequency = 0.0;

			Position groundPosition;
			double lat, lon, height;
			lat = GCPS.at<double>(i, 3);
			lon = GCPS.at<double>(i, 2);
			lon = lon > 180.0 ? (lon - 360.0) : lon;
			height = GCPS.at<double>(i, 4);
			Utils::ell2xyz(lon, lat, height, groundPosition);
			int numOrbitVec = stateVectors.newStateVectors.rows;
			double firstVecTime = 0.0;
			double secondVecTime = 0.0;
			double firstVecFreq = 0.0;
			double secondVecFreq = 0.0;
			double currentFreq, xdiff, ydiff, zdiff, distance = 1.0, zeroDopplerTime;
			//检测标志点位于哪两个轨道点之间
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
			//对两个点之间（相差10s）进行进一步插值检测找到标志点对应的具体卫星位置
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
			zeroDopplerTime = lowerBoundTime + lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);
			stateVectors.getPosition(zeroDopplerTime, pos);
			sate1.at<double>(0) = pos.x;
			sate1.at<double>(1) = pos.y;
			sate1.at<double>(2) = pos.z;
			double r;
			Mat XYZ, LLH(1, 3, CV_64F), tt;
			LLH.at<double>(0, 0) = GCPS.at<double>(i, 3);
			LLH.at<double>(0, 1) = GCPS.at<double>(i, 2);
			LLH.at<double>(0, 2) = GCPS.at<double>(i, 4);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate1;
			r = cv::norm(tt, cv::NORM_L2);
			double master_phase = -r / wavelength * 4 * PI;

			//辅卫星斜距
			Mat sate2 = Mat::zeros(1, 3, CV_64F);
			orbitStateVectors stateVectors2(statevec2, start2, end2);
			stateVectors2.applyOrbit();

			dopplerFrequency = 0.0;

			numOrbitVec = stateVectors2.newStateVectors.rows;
			firstVecTime = 0.0;
			secondVecTime = 0.0;
			firstVecFreq = 0.0;
			secondVecFreq = 0.0;
			distance = 1.0;
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
				currentFreq = 2.0 * (xdiff * orb_vel.vx + ydiff * orb_vel.vy + zdiff * orb_vel.vz) / (wavelength2 * distance);
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
			absLineTimeInterval = 1.0 / prf2;

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
				midFreq = 2.0 * (xdiff * vel.vx + ydiff * vel.vy + zdiff * vel.vz) / (wavelength2 * distance);
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
			zeroDopplerTime = lowerBoundTime + lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);

			double r2;
			stateVectors2.getPosition(zeroDopplerTime, pos);
			sate2.at<double>(0) = pos.x;
			sate2.at<double>(1) = pos.y;
			sate2.at<double>(2) = pos.z;
			LLH.at<double>(0, 0) = GCPS.at<double>(i, 3);
			LLH.at<double>(0, 1) = GCPS.at<double>(i, 2);
			LLH.at<double>(0, 2) = GCPS.at<double>(i, 4);
			util.ell2xyz(LLH, XYZ);
			tt = XYZ - sate2;
			r2 = cv::norm(tt, cv::NORM_L2);
			double slave_phase = -r2 / wavelength2 * 4 * PI;

			//先验相位
			Inphase_prior = atan2(sin(master_phase - slave_phase), cos(master_phase - slave_phase));
			GCPS_New.at<double>(i, 5) = Inphase_prior;
		}
		PhaseError += pow(atan2(sin(Inphase_unwrapped - Inphase_prior), cos(Inphase_unwrapped - Inphase_prior)), 2);
	}
	*Output = sqrt(PhaseError / Gcps_number);
	return 0;
}

int Evaluation::FFT2(ComplexMat src, ComplexMat& dst, int win_size, int interp_times)
{
	Point master_max, slave_max;
	double max1 = 0, max2;
	Utils util;
	int interp_size = win_size * interp_times;
	ComplexMat tmp, tmp_fft, tmp_fft_half;
	ComplexMat tmp_interp_fft(1, interp_size), tmp_interp(1, interp_size),
		tmp_interp_fft2(interp_size, 1), tmp_interp2(interp_size,1);
	ComplexMat tmp_fft1(win_size, interp_size), tmp_fft2(interp_size, interp_size);
	for (int i = 0; i < win_size; i++)
	{
		tmp = src(Range(i, i + 1), Range(0, win_size));
		util.fft2(tmp, tmp_fft);
		tmp_fft = tmp_fft * interp_times;
		tmp_interp_fft = ComplexMat(1, interp_size);
		tmp_fft_half = tmp_fft(Range(0, 1), Range(0, win_size / 2));
		tmp_interp_fft.SetValue(Range(0, 1), Range(0, win_size / 2), tmp_fft_half);
		tmp_fft_half = tmp_fft(Range(0, 1), Range(win_size/2, win_size));
		tmp_interp_fft.SetValue(Range(0, 1), Range(interp_size-win_size/2, interp_size), tmp_fft_half);
		util.ifft2(tmp_interp_fft, tmp_interp);
		tmp_interp = tmp_interp * ((double)1 / interp_times / interp_times);
		tmp_fft1.SetValue(Range(i, i + 1), Range(0, interp_size),tmp_interp);
		minMaxLoc(tmp_interp.GetMod(), NULL, &max1, NULL, &master_max);
	}

	for (int i = 0; i < interp_size; i++)
	{
		tmp = tmp_fft1(Range(0, win_size), Range(i, i+1));
		util.fft2(tmp, tmp_fft);
		tmp_fft = tmp_fft * interp_times;
		tmp_interp_fft2 = ComplexMat(interp_size,1);
		tmp_fft_half = tmp_fft(Range(0, win_size / 2), Range(0, 1));
		tmp_interp_fft2.SetValue(Range(0, win_size / 2), Range(0, 1), tmp_fft_half);
		tmp_fft_half = tmp_fft( Range(win_size / 2, win_size), Range(0, 1));
		tmp_interp_fft2.SetValue(Range(interp_size - win_size / 2, interp_size), Range(0, 1), tmp_fft_half);
		util.ifft2(tmp_interp_fft2, tmp_interp2);
		tmp_interp2 = tmp_interp2 * ((double)1 / interp_times / interp_times);
		tmp_fft2.SetValue(Range(0, interp_size), Range(i, i + 1), tmp_interp2);
	}
	
	minMaxLoc(tmp_fft2.GetMod(), NULL, &max1, NULL, &master_max);
	dst = tmp_fft2;
	return 0;
}

int Evaluation::Pos(const char* unwrapped_phase_file, const char* project_path, const char* GCP_path, double* lat_abs, double* lat_rel, double* lon_abs, double* lon_rel, double* height_abs, double* height_rel)
{
	FormatConversion conversion; Utils util;
	int nr, nc, ret, offset_row, offset_col;
	double time_interval1, time_interval2;
	string source_1, source_2, tmp;
	string project(project_path);
	Mat unwrapped_phase, flat_phase_coefficient, gcps, temp, range_spacing,
		stateVec1, stateVec2, lat_coefficient, lon_coefficient, prf1, prf2, carrier_frequency;
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "phase", unwrapped_phase);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	nr = unwrapped_phase.rows; nc = unwrapped_phase.cols;
	if (nr < 1 || nc < 1)
	{
		fprintf(stderr, "dem_newton_iter(): invalide unwrapped_phase !\n");
		return -1;
	}
	ret = conversion.read_array_from_h5(unwrapped_phase_file, "flat_phase_coefficient", flat_phase_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_1", source_1);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	ret = conversion.read_str_from_h5(unwrapped_phase_file, "source_2", source_2);
	if (return_check(ret, "read_str_from_h5()", error_head)) return -1;
	tmp = project;
	source_1 = tmp + source_1;
	source_2 = tmp + source_2;
	ret = conversion.read_array_from_h5(GCP_path, "GCP", gcps);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_row", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_row = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "offset_col", temp);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	offset_col = temp.at<int>(0, 0);
	ret = conversion.read_array_from_h5(source_1.c_str(), "range_spacing", range_spacing);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "state_vec", stateVec1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_2.c_str(), "state_vec", stateVec2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "prf", prf1);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval1 = 1.0 / (prf1.at<double>(0, 0) + 1e-10);
	ret = conversion.read_array_from_h5(source_2.c_str(), "prf", prf2);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	time_interval2 = 1.0 / (prf2.at<double>(0, 0) + 1e-10);
	ret = conversion.read_array_from_h5(source_1.c_str(), "lat_coefficient", lat_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "lon_coefficient", lon_coefficient);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;
	ret = conversion.read_array_from_h5(source_1.c_str(), "carrier_frequency", carrier_frequency);
	if (return_check(ret, "read_array_from_h5()", error_head)) return -1;

	//寻找图像范围内的控制点信息

	int num_gcps = gcps.rows;
	int row, col, i = 0;
	bool b_gcp = false;
	for (i = 0; i < num_gcps; i++)
	{
		row = (int)gcps.at<double>(i, 0);
		col = (int)gcps.at<double>(i, 1);
		if ((row - offset_row) >= 0 && (row - offset_row) < nr && (col - offset_col) >= 0 && (col - offset_col) < nc)
		{
			b_gcp = true;
			break;
		}
	}

	Mat llh(1, 3, CV_64F), xyz_ground(1, 3, CV_64F);
	Mat row_coord(1, 1, CV_64F), col_coord(1, 1, CV_64F), lat, lon;
	if (b_gcp)
	{
		llh.at<double>(0, 0) = gcps.at<double>(i, 3);
		llh.at<double>(0, 1) = gcps.at<double>(i, 2);
		llh.at<double>(0, 2) = gcps.at<double>(i, 4);
		row = row - offset_row; col = col - offset_col;
	}
	else
	{
		row_coord.at<double>(0, 0) = offset_row + double(nr) / 2.0;
		col_coord.at<double>(0, 0) = offset_col + double(nc) / 2.0;
		ret = util.coord_conversion(lat_coefficient, row_coord, col_coord, lat);
		if (return_check(ret, "coord_conversion", error_head)) return -1;
		ret = util.coord_conversion(lon_coefficient, row_coord, col_coord, lon);
		if (return_check(ret, "coord_conversion", error_head)) return -1;
		llh.at<double>(0, 0) = lat.at<double>(0, 0);
		llh.at<double>(0, 1) = lon.at<double>(0, 0);
		llh.at<double>(0, 2) = 0.0;
		row = (int)nr / 2; col = (int)nc / 2;
	}
	ret = util.ell2xyz(llh, xyz_ground);
	if (return_check(ret, "ell2xyz", error_head)) return -1;


	/*
	* 轨道插值
	*/
	ret = util.stateVec_interp(stateVec1, time_interval1, stateVec1);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	ret = util.stateVec_interp(stateVec2, time_interval2, stateVec2);
	if (return_check(ret, "stateVec_interp()", error_head)) return -1;
	/*
	* 寻找图像左上角成像卫星位置
	*/
	Mat sate1_xyz, sate2_xyz, sate1_v, sate2_v;
	Mat sate1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat satev1 = Mat::zeros(nr, 3, CV_64F);//主星位置
	Mat sate2 = Mat::zeros(nr, 3, CV_64F);//辅星位置
	stateVec1(cv::Range(0, stateVec1.rows), cv::Range(1, 4)).copyTo(sate1_xyz);
	stateVec1(cv::Range(0, stateVec1.rows), cv::Range(4, 7)).copyTo(sate1_v);
	stateVec2(cv::Range(0, stateVec2.rows), cv::Range(1, 4)).copyTo(sate2_xyz);
	stateVec2(cv::Range(0, stateVec2.rows), cv::Range(4, 7)).copyTo(sate2_v);
	//找到零多普勒位置
	Mat xyz, llh_upperleft(1, 3, CV_64F);
	row_coord.at<double>(0, 0) = offset_row;
	col_coord.at<double>(0, 0) = offset_col + double(nc) / 2.0;
	ret = util.coord_conversion(lat_coefficient, row_coord, col_coord, lat);
	if (return_check(ret, "coord_conversion", error_head)) return -1;
	ret = util.coord_conversion(lon_coefficient, row_coord, col_coord, lon);
	if (return_check(ret, "coord_conversion", error_head)) return -1;
	llh_upperleft.at<double>(0, 0) = lat.at<double>(0, 0);
	llh_upperleft.at<double>(0, 1) = lon.at<double>(0, 0);
	llh_upperleft.at<double>(0, 2) = 0.0;
	ret = util.ell2xyz(llh_upperleft, xyz);
	if (return_check(ret, "ell2xyz", error_head)) return -1;
	Mat dop = Mat::zeros(sate1_xyz.rows, 1, CV_64F);
	Mat r;
	for (int j = 0; j < sate1_xyz.rows; j++)
	{
		r = xyz - sate1_xyz(Range(j, j + 1), Range(0, 3));
		dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate1_v(Range(j, j + 1), Range(0, 3))))[0]);
	}
	Point peak_loc;
	cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
	int xxxx;
	for (int j = 0; j < nr; j++)
	{
		xxxx = (peak_loc.y + j) > (sate1_xyz.rows - 1) ? (sate1_xyz.rows - 1) : (peak_loc.y + j);
		sate1_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate1(Range(j, j + 1), Range(0, 3)));
		sate1_v(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(satev1(Range(j, j + 1), Range(0, 3)));
	}
	//卫星2
	dop = Mat::zeros(sate2_xyz.rows, 1, CV_64F);
	for (int j = 0; j < sate2_xyz.rows; j++)
	{
		r = xyz - sate2_xyz(Range(j, j + 1), Range(0, 3));
		dop.at<double>(j, 0) = fabs(cv::sum(r.mul(sate2_v(Range(j, j + 1), Range(0, 3))))[0]);
	}
	cv::minMaxLoc(dop, NULL, NULL, &peak_loc, NULL);
	for (int j = 0; j < nr; j++)
	{
		xxxx = (peak_loc.y + j) > (sate2_xyz.rows - 1) ? (sate2_xyz.rows - 1) : (peak_loc.y + j);
		sate2_xyz(Range(xxxx, xxxx + 1), Range(0, 3)).copyTo(sate2(Range(j, j + 1), Range(0, 3)));
	}
	//加回平地相位
	//Mat coef;
	//cv::transpose(flat_phase_coefficient, coef);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nr; i++)
	{
		Mat temp(1, 6, CV_64F);
		for (int j = 0; j < nc; j++)
		{
			temp.at<double>(0, 0) = 1.0;
			temp.at<double>(0, 1) = i;
			temp.at<double>(0, 2) = j;
			temp.at<double>(0, 3) = i * j;
			temp.at<double>(0, 4) = i * i;
			temp.at<double>(0, 5) = j * j;
			unwrapped_phase.at<double>(i, j) = unwrapped_phase.at<double>(i, j) + sum(temp.mul(flat_phase_coefficient))[0];
		}
	}
	//控制点绝对相位计算
	double r_main = sqrt(sum((sate1(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground).mul(sate1(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground))[0]);
	double r_slave = sqrt(sum((sate2(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground).mul(sate2(cv::Range(row, row + 1), cv::Range(0, 3)) - xyz_ground))[0]);
	int mode = 1;
	double C = 4*PI; //自发自收
	double lambda = 3e8 / (carrier_frequency.at<double>(0, 0) + 1e-10);
	double phase_real = (r_slave - r_main) / lambda * C;
	double K = round((phase_real - unwrapped_phase.at<double>(row, col)) / (2 * PI));
	unwrapped_phase = unwrapped_phase + K * 2 * PI;//相位校正

	//util.cvmat2bin("E:\\working_dir\\projects\\software\\InSAR\\bin\\unwrapped_phase_abs.bin", unwrapped_phase);

	/*
	* 反演高程
	*/

	Mat R_M(1, nc, CV_64F);
	for (int i = 0; i < nc; i++)
	{
		R_M.at<double>(0, i) = r_main + range_spacing.at<double>(0, 0) * (double(i) - col);
	}
	Mat ones = Mat::ones(nr, 1, CV_64F);
	R_M = ones * R_M;
	Mat R_F = R_M * 2.0 + lambda * unwrapped_phase / (2 * PI);
	Mat Satellite_M_T_Position = sate1;//主星发射位置
	Mat Satellite_S_T_Position;
	if (mode == 1)
	{
		Satellite_S_T_Position = sate2;//辅星发射位置
	}
	else
	{
		Satellite_S_T_Position = sate1;//辅星发射位置
	}
	Mat Satellite_M_R_Position = sate1;//主星接收位置
	Mat Satellite_S_R_Position = sate2;//辅星接收位置
	Mat Satellite_M = (Satellite_M_R_Position + Satellite_M_T_Position) / 2;
	Mat Vs = satev1;
	ones = Mat::ones(nr, nc, CV_64F);
	Mat P1 = ones * xyz_ground.at<double>(0, 0);
	Mat P2 = ones * xyz_ground.at<double>(0, 1);
	Mat P3 = ones * xyz_ground.at<double>(0, 2);
	Mat M_T, M_R, S_T, S_R;
	Mat f1, f2, f3;
	Mat det_Df, Df_ni11, Df_ni12, Df_ni13, Df_ni21, Df_ni22, Df_ni23;
	Mat Df11, Df12, Df13, Df21, Df22, Df23, Df31, Df32, Df33, Df_ni31, Df_ni32, Df_ni33;
	Mat delta_Rt1, delta_Rt2, delta_Rt3;
	ones = Mat::ones(1, nc, CV_64F);
	Mat temp_var, temp_var1;
	Mat fd = Mat::zeros(1, nc, CV_64F);
	for (int i = 0; i < 15; i++)
	{
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(M_T);
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		M_T = M_T + temp_var;
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		M_T = M_T + temp_var;
		cv::sqrt(M_T, f1);
		f1 = f1 * 2 - 2 * R_M;

		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(S_T);
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		S_T = S_T + temp_var;
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		S_T = S_T + temp_var;

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(0, 1)) * ones - P1;
		temp_var = temp_var.mul(temp_var);
		temp_var.copyTo(S_R);
		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(1, 2)) * ones - P2;
		temp_var = temp_var.mul(temp_var);
		S_R = S_R + temp_var;
		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(2, 3)) * ones - P3;
		temp_var = temp_var.mul(temp_var);
		S_R = S_R + temp_var;

		cv::sqrt(S_T, f2);
		cv::sqrt(S_R, temp_var);
		f2 = f2 + temp_var - R_F;


		temp_var = Vs(Range(0, Vs.rows), Range(0, 1)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(0, 1)) * ones - P1;
		f3 = temp_var.mul(temp_var1);
		temp_var = Vs(Range(0, Vs.rows), Range(1, 2)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(1, 2)) * ones - P2;
		f3 = f3 + temp_var.mul(temp_var1);
		temp_var = Vs(Range(0, Vs.rows), Range(2, 3)) * ones;
		temp_var1 = Satellite_M(Range(0, Satellite_M.rows), Range(2, 3)) * ones - P3;
		f3 = f3 + temp_var.mul(temp_var1);
		ones = Mat::ones(nr, 1, CV_64F);
		temp_var = ones * fd;
		temp_var1 = R_M * lambda / 2.0;
		f3 = f3 + temp_var.mul(temp_var1);

		//Dff
		//第一行：f(1)的x，y，z的导数
		ones = Mat::ones(1, nc, CV_64F);
		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df11 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df11 = Df11 + temp_var.mul(temp_var1);

		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df12 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df12 = Df12 + temp_var.mul(temp_var1);

		temp_var = Satellite_M_T_Position(Range(0, Satellite_M_T_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df13 = temp_var.mul(temp_var1);

		temp_var = Satellite_M_R_Position(Range(0, Satellite_M_R_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(M_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df13 = Df13 + temp_var.mul(temp_var1);


		//第二行：f(2)的x，y，z的导数
		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df21 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(0, 1)) * ones - P1;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df21 = Df21 + temp_var.mul(temp_var1);


		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df22 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(1, 2)) * ones - P2;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df22 = Df22 + temp_var.mul(temp_var1);


		temp_var = Satellite_S_T_Position(Range(0, Satellite_S_T_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(S_T, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df23 = temp_var.mul(temp_var1);

		temp_var = Satellite_S_R_Position(Range(0, Satellite_S_R_Position.rows), Range(2, 3)) * ones - P3;
		cv::sqrt(S_R, temp_var1);
		temp_var1 = 1 / temp_var1;
		temp_var1 = -temp_var1;
		Df23 = Df23 + temp_var.mul(temp_var1);

		//第三行：f(3)的x，y，z的导数
		Df31 = -Vs(Range(0, Vs.rows), Range(0, 1)) * ones;
		Df32 = -Vs(Range(0, Vs.rows), Range(1, 2)) * ones;
		Df33 = -Vs(Range(0, Vs.rows), Range(2, 3)) * ones;

		temp_var = Df11.mul(Df22);
		temp_var = temp_var.mul(Df33);
		temp_var.copyTo(det_Df);

		temp_var = Df12.mul(Df23);
		temp_var = temp_var.mul(Df31);
		det_Df = det_Df + temp_var;

		temp_var = Df13.mul(Df21);
		temp_var = temp_var.mul(Df32);
		det_Df = det_Df + temp_var;

		temp_var = Df31.mul(Df22);
		temp_var = temp_var.mul(Df13);
		det_Df = det_Df - temp_var;

		temp_var = Df32.mul(Df23);
		temp_var = temp_var.mul(Df11);
		det_Df = det_Df - temp_var;

		temp_var = Df33.mul(Df21);
		temp_var = temp_var.mul(Df12);
		det_Df = det_Df - temp_var;

		Df_ni11 = (Df22.mul(Df33) - Df32.mul(Df23)) / det_Df;
		Df_ni12 = -(Df12.mul(Df33) - Df32.mul(Df13)) / det_Df;
		Df_ni13 = (Df12.mul(Df23) - Df22.mul(Df13)) / det_Df;
		delta_Rt1 = Df_ni11.mul(f1) + Df_ni12.mul(f2) + Df_ni13.mul(f3);


		Df_ni21 = -(Df21.mul(Df33) - Df31.mul(Df23)) / det_Df;
		Df_ni22 = (Df11.mul(Df33) - Df31.mul(Df13)) / det_Df;
		Df_ni23 = -(Df11.mul(Df23) - Df21.mul(Df13)) / det_Df;
		delta_Rt2 = Df_ni21.mul(f1) + Df_ni22.mul(f2) + Df_ni23.mul(f3);

		Df_ni31 = (Df21.mul(Df32) - Df31.mul(Df22)) / det_Df;
		Df_ni32 = -(Df11.mul(Df32) - Df31.mul(Df12)) / det_Df;
		Df_ni33 = (Df22.mul(Df11) - Df21.mul(Df12)) / det_Df;
		delta_Rt3 = Df_ni31.mul(f1) + Df_ni32.mul(f2) + Df_ni33.mul(f3);

		P1 = P1 - delta_Rt1;
		P2 = P2 - delta_Rt2;
		P3 = P3 - delta_Rt3;
	}
	delta_Rt1.release(); delta_Rt2.release(); delta_Rt3.release(); Df_ni31.release(); Df_ni32.release();
	Df_ni33.release(); Df_ni21.release(); Df_ni22.release(); Df_ni23.release(); Df_ni11.release(); Df_ni12.release();
	Df_ni13.release();
	Mat GCPs;
	ret = conversion.read_array_from_h5(GCP_path, "GCP", GCPs);
	int GCP_count = GCPs.rows;
	int count = 0, count2 = 0;
	double lat_1 = 0, lat_2 = 0, lon_1 = 0, lon_2 = 0, height_1 = 0, height_2 = 0;
	double height_offset = 0;
	for (int i = 0; i < GCP_count; i++)
	{
		int row1 = GCPs.at<double>(i, 0);
		int col1 = GCPs.at<double>(i, 1);
		if (row1 >= offset_row && row1 < offset_row + nr && col1 >= offset_col && col1 < offset_col + nc)
		{
			Mat xyz, llh, llh2;
			xyz = Mat::zeros(1, 3, CV_64F);
			xyz.at<double>(0, 0) = P1.at<double>(row1 - offset_row, col1 - offset_col);
			xyz.at<double>(0, 1) = P2.at<double>(row1 - offset_row, col1 - offset_col);
			xyz.at<double>(0, 2) = P3.at<double>(row1 - offset_row, col1 - offset_col);

			ret = util.xyz2ell(xyz, llh);
			height_offset = llh.at<double>(0, 2);
			lat_1 += GCPs.at<double>(i, 3) - llh.at<double>(0, 0);
			lon_1 += GCPs.at<double>(i, 2) - llh.at<double>(0, 1);
			height_1 += GCPs.at<double>(i, 4) - llh.at<double>(0, 2);
			count++;
			for (int j = i; j < GCP_count; j++)
			{
				int row2 = GCPs.at<double>(j, 0);
				int col2 = GCPs.at<double>(j, 1);
				if (row2 >= offset_row && row2 < offset_row + nr && col2 >= offset_col && col2 < offset_col + nc)
				{
					xyz.at<double>(0, 0) = P1.at<double>(row2 - offset_row, col2 - offset_col);
					xyz.at<double>(0, 1) = P2.at<double>(row2 - offset_row, col2 - offset_col);
					xyz.at<double>(0, 2) = P3.at<double>(row2 - offset_row, col2 - offset_col);

					ret = util.xyz2ell(xyz, llh2);
					lat_2 += GCPs.at<double>(i, 3) - llh.at<double>(0, 0) - (GCPs.at<double>(j, 3) - llh2.at<double>(0, 0));
					lon_2 += GCPs.at<double>(i, 2) - llh.at<double>(0, 1) - (GCPs.at<double>(j, 2) - llh2.at<double>(0, 1));
					height_2 += GCPs.at<double>(i, 4) - llh.at<double>(0, 2) - (GCPs.at<double>(j, 4) - llh2.at<double>(0, 2));
					count2++;
				}
				else continue;
			}

		}
		else continue;
	}
	height_1 = height_1 + count * (llh.at<double>(0, 2) - height_offset);
		lon_1 /= count; lat_1 /= count;
		lon_2 /= count2; lat_2 /= count2;
		height_1 /= count; height_2 /= count2;
	*lon_abs = lon_1; *lon_rel = lon_2; *lat_abs = lat_1; *lat_rel = lat_2; *height_abs = height_1; *height_rel = height_2;
	return 0;
}
