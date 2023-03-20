#include"stdafx.h"
#include"Eigen/Dense"
#include"..\include\ComplexMat.h"

ComplexMat::ComplexMat()
{
	Mat tmp = Mat::zeros(1, 1, CV_64F);
	tmp.copyTo(this->re);
	tmp.copyTo(this->im);
}

ComplexMat::ComplexMat(Mat& re, Mat& im)
{
	//this->re = re;
	//this->im = im;
	re.copyTo(this->re);
	im.copyTo(this->im);
}

ComplexMat::ComplexMat(int rows, int cols)
{
	if (rows > 0 && cols > 0)
	{
		Mat tmp = Mat::zeros(rows, cols, CV_64F);
		tmp.copyTo(this->im);
		tmp.copyTo(this->re);
	}
}

ComplexMat::ComplexMat(const ComplexMat& b)
{
	b.re.copyTo(this->re);
	b.im.copyTo(this->im);
}

ComplexMat::~ComplexMat()
{

}

void ComplexMat::release()
{
	this->re.release();
	this->im.release();
}

int ComplexMat::type() const
{
	if (re.type() != im.type()) return -1;
	return this->re.type();
}

int ComplexMat::Mul(const ComplexMat& Src, ComplexMat& Dst, bool bConj) const
{
	ComplexMat result;
	if (this->GetRows() < 1 ||
		this->GetCols() < 1 ||
		Src.GetRows() < 1 ||
		Src.GetCols() < 1||
		this->GetCols() != Src.GetCols()||
		this->GetRows() != Src.GetRows()||
		this->type() != Src.type())
	{
		fprintf(stderr, "ComplexMat::Mul(): input check failed!\n\n");
		return -1;
	}

	if (bConj)
	{
		result.re = this->re.mul(Src.re) + this->im.mul(Src.im);
		result.im = Src.re.mul(this->im) - this->re.mul(Src.im);
	}
	else
	{
		result.re = this->re.mul(Src.re) - this->im.mul(Src.im);
		result.im = Src.re.mul(this->im) + this->re.mul(Src.im);
	}
	Dst.SetRe(result.re);
	Dst.SetIm(result.im);
	return 0;
}

ComplexMat ComplexMat::operator*(const ComplexMat& b) const
{
	ComplexMat result;
	if (this->GetRows() < 1 ||
		this->GetCols() < 1 ||
		b.GetRows() < 1 ||
		b.GetCols() < 1 ||
		(this->GetCols() != b.GetCols()) && b.GetCols() != 1 ||
		(this->GetRows() != b.GetRows()) && b.GetRows() != 1 ||
		this->type() != b.type())
	{
		fprintf(stderr, "ComplexMat::Mul(): input check failed!\n\n");
		return ComplexMat();
	}
	if (b.GetCols() == 1 && b.GetRows() == 1)
	{
		result.re = this->re * b.re.at<double>(0, 0) - this->im * b.im.at<double>(0, 0);
		result.im = this->im * b.re.at<double>(0, 0) + this->re * b.im.at<double>(0, 0);
	}
	else
	{
		result.re = this->re.mul(b.re) - this->im.mul(b.im);
		result.im = b.re.mul(this->im) + this->re.mul(b.im);
	}
	return result;
}

ComplexMat ComplexMat::operator*(const Mat& a) const
{
	if (a.cols != this->GetCols() ||
		a.rows != this->GetRows()||
		a.type() != this->type()||
		a.channels() != 1)
	{
		fprintf(stderr, "ComplexMat::operator*(const Mat& a): input check failed!\n\n");
		return ComplexMat();
	}
	Mat out_re, out_im;
	ComplexMat out;
	out_re = this->re.mul(a);
	out_im = this->im.mul(a);
	out.SetRe(out_re);
	out.SetIm(out_im);
	return out;
}

ComplexMat ComplexMat::operator*(const double& a) const
{
	Mat out_re, out_im;
	ComplexMat out;
	out_re = this->re * a;
	out_im = this->im * a;
	out.SetRe(out_re);
	out.SetIm(out_im);
	return out;
}

ComplexMat ComplexMat::operator()(cv::Range _rowRange, cv::Range _colRange) const
{
	if (_rowRange.start < 0 ||
		_rowRange.end > this->GetRows() ||
		_colRange.start < 0 ||
		_colRange.end > this->GetCols())
	{
		fprintf(stderr, "ComplexMat::operator()(cv::Range _rowRange, cv::Range _colRange): \n Range exceeds legal value!\n\n");
		return ComplexMat();
	}
	ComplexMat out;
	this->re(cv::Range(_rowRange.start, _rowRange.end), cv::Range(_colRange.start, _colRange.end)).copyTo(out.re);
	this->im(cv::Range(_rowRange.start, _rowRange.end), cv::Range(_colRange.start, _colRange.end)).copyTo(out.im);
	return out;
}

int ComplexMat::SetValue(cv::Range _rowRange, cv::Range _colRange, ComplexMat& src)
{
	if ((_rowRange.end - _rowRange.start) != src.GetRows() ||
		(_colRange.end - _colRange.start) != src.GetCols() ||
		_rowRange.start < 0 ||
		_rowRange.end > this->GetRows() ||
		_colRange.start < 0 ||
		_colRange.end > this->GetCols()||
		src.type() != this->type()||
		(src.type() != CV_64F && src.type() != CV_16S)
		)
	{
		fprintf(stderr, "ComplexMat::SetValue(): input check failed!\n\n");
		return -1;
	}
	if (src.type() == CV_64F)
	{
		for (int i = _rowRange.start; i < _rowRange.end; i++)
		{
			for (int j = _colRange.start; j < _colRange.end; j++)
			{

				this->re.at<double>(i, j) = src.re.at<double>(i - _rowRange.start, j - _colRange.start);
				this->im.at<double>(i, j) = src.im.at<double>(i - _rowRange.start, j - _colRange.start);
			}
		}
	}
	if (src.type() == CV_16S)
	{
		for (int i = _rowRange.start; i < _rowRange.end; i++)
		{
			for (int j = _colRange.start; j < _colRange.end; j++)
			{

				this->re.at<short>(i, j) = src.re.at<short>(i - _rowRange.start, j - _colRange.start);
				this->im.at<short>(i, j) = src.im.at<short>(i - _rowRange.start, j - _colRange.start);
			}
		}
	}
	
	return 0;
}

Mat ComplexMat::GetIm() const
{
	return this->im;
}

Mat ComplexMat::GetMod() const
{
	Mat tmp;
	magnitude(this->re, this->im, tmp);
	return tmp;
}

Mat ComplexMat::GetPhase()
{
	int nr = this->GetRows();
	int nc = this->GetCols();
	if (nr < 1 || nc < 1)
	{
		return Mat::zeros(1, 1, CV_64F);
	}
	Mat phase(nr, nc, CV_64F);
	if (this->type() == CV_64F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				phase.at <double>(i, j) = atan2(this->im.at<double>(i, j), this->re.at<double>(i, j));
			}
		}
	}
	else if (this->type() == CV_32F)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				phase.at <double>(i, j) = atan2(this->im.at<float>(i, j), this->re.at<float>(i, j));
			}
		}
	}
	else if (this->type() == CV_32S)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				phase.at <double>(i, j) = atan2((double)this->im.at<int>(i, j), (double)this->re.at<int>(i, j));
			}
		}
	}
	else if (this->type() == CV_16S)
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				phase.at <double>(i, j) = atan2((double)this->im.at<short>(i, j), (double)this->re.at<short>(i, j));
			}
		}
	}

	return phase;
}

Mat ComplexMat::GetRe() const
{
	return this->re;
}

void ComplexMat::SetIm(Mat& im)
{
	im.copyTo(this->im);

}

void ComplexMat::SetRe(Mat& re)
{
	re.copyTo(this->re);
}

int ComplexMat::GetCols() const
{
	if (re.cols != im.cols) return -1;
	return this->re.cols;
}

int ComplexMat::mul(const ComplexMat& Src, ComplexMat& Dst, bool bConj)
{
	if (Src.isempty() ||
		this->isempty() ||
		Src.type() != this->type() ||
		this->GetCols() != Src.GetRows()
		)
	{
		fprintf(stderr, "ComplexMat::mul(): input check failed!\n");
		return -1;
	}
	if (!bConj)
	{
		Dst.re = this->re * Src.re - this->im * Src.im;
		Dst.im = this->re * Src.im + this->im * Src.re;
	}
	else
	{
		Dst.re = this->re * Src.re + this->im * Src.im;
		Dst.im = this->im * Src.re - this->re * Src.im;
	}

	return 0;
}

int ComplexMat::GetRows() const
{
	if (re.rows != im.rows) return -1;
	return this->re.rows;
}

ComplexMat ComplexMat::operator+(const ComplexMat& b) const
{
	if (this->GetCols() != b.GetCols() ||
		b.GetRows() != b.GetRows() ||
		this->type() != b.type()||
		this->GetCols() < 1||
		this->GetRows() < 1
		)
	{
		fprintf(stderr, "ComplexMat operator+: input check failed!\n\n");
		return *this;
	}
	ComplexMat out;
	Mat out_re, out_im;
	out_re = this->re + b.re;
	out_im = this->im + b.im;
	out.SetIm(out_im);
	out.SetRe(out_re);
	return out;
}

ComplexMat ComplexMat::operator=(const ComplexMat& b)
{
	b.re.copyTo(this->re);
	b.im.copyTo(this->im);
	return *this;
}

ComplexMat ComplexMat::sum(int dim) const
{
	ComplexMat out;
	Mat re, im;
	int nr = this->GetRows();
	int nc = this->GetCols();
	if (nr < 0 || nc < 0)
	{
		fprintf(stderr, "ComplexMat::sum(): rows or cols < 0!\n\n");
		return ComplexMat();
	}
	if (dim == 1)
	{
		re = Mat::zeros(nr, 1, CV_64F);
		im = Mat::zeros(nr, 1, CV_64F);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nr; i++)
		{
			double tmp_re, tmp_im;
			tmp_re = 0.0; tmp_im = 0.0;
			for (int j = 0; j < nc; j++)
			{
				tmp_re += this->re.at<double>(i, j);
				tmp_im += this->im.at<double>(i, j);
			}
			re.at<double>(i, 0) = tmp_re;
			im.at<double>(i, 0) = tmp_im;
		}
	}
	else
	{
		re = Mat::zeros(1, nc, CV_64F);
		im = Mat::zeros(1, nc, CV_64F);
		
#pragma omp parallel for schedule(guided)
		for (int j = 0; j < nc; j++)
		{
			double tmp_re, tmp_im;
			tmp_re = 0.0; tmp_im = 0.0;
			for (int i = 0; i < nr; i++)
			{
				tmp_re += this->re.at<double>(i, j);
				tmp_im += this->im.at<double>(i, j);
			}
			re.at<double>(0, j) = tmp_re;
			im.at<double>(0, j) = tmp_im;
		}
	}
	out.SetRe(re);
	out.SetIm(im);
	return out;
}

complex<double> ComplexMat::determinant() const
{
	int nr, nc;
	nr = this->GetRows();
	nc = this->GetCols();
	Eigen::MatrixXcd x(nr, nc);
	complex<double> d;
	if (this->type() == CV_32F)
	{
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				d.real(this->re.at<float>(i, j));
				d.imag(this->im.at<float>(i, j));
				x(i, j) = d;
			}
		}
	}
	else
	{
		for (int i = 0; i < nr; i++)
		{
			for (int j = 0; j < nc; j++)
			{
				d.real(this->re.at<double>(i, j));
				d.imag(this->im.at<double>(i, j));
				x(i, j) = d;
			}
		}
	}
	
	return x.determinant();
}

ComplexMat ComplexMat::conj() const
{
	ComplexMat out;
	Mat im, re;
	this->re.copyTo(re);
	this->im.copyTo(im);
	im = -im;
	out.SetRe(re);
	out.SetIm(im);
	return out;
}

ComplexMat ComplexMat::transpose(bool conj) const
{
	ComplexMat out;
	Mat im, re;
	this->re.copyTo(re);
	this->im.copyTo(im);
	if (conj) im = -im;
	cv::transpose(re, re);
	cv::transpose(im, im);
	out.SetRe(re);
	out.SetIm(im);
	return out;
}

int ComplexMat::reshape(int rows, int cols, ComplexMat& dst)
{
	if (rows * cols != this->GetCols() * this->GetRows())
	{
		fprintf(stderr, "ComplexMat::reshape(): input check failed!\n");
		return -1;
	}
	Mat real, imagine;
	cv::transpose(this->re, real);
	cv::transpose(this->im, imagine);
	real = real.reshape(1, rows);
	imagine = imagine.reshape(1, rows);
	real.copyTo(dst.re);
	imagine.copyTo(dst.im);
	return 0;
}

int ComplexMat::countNonzero() const
{
	if (this->re.rows == 0 ||
		this->im.rows == 0 ||
		this->re.cols == 0 ||
		this->im.cols == 0
		)
	{
		return 0;
	}
	int count = 0;
	int nr = GetRows();
	int nc = GetCols();
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (fabs(this->re.at<double>(i, j)) > DBL_EPSILON || fabs(this->im.at<double>(i, j)) > DBL_EPSILON)
			{
				count++;
			}
		}
	}
	return count;
}

bool ComplexMat::isempty() const
{
	if (this->GetRows() < 1 || this->GetCols() < 1) return true;
	return false;
}

void ComplexMat::convertTo(ComplexMat& out, int type) const 
{
	if (this->type() == type)
	{
		out = *this;
	}
	else
	{
		this->re.convertTo(out.re, type);
		this->im.convertTo(out.im, type);
	}
}
