#pragma once
#ifndef __PACKAGE__H__
#define __PACKAGE__H__

#include"opencv2\core\core.hpp"
#include"opencv2\highgui\highgui.hpp"
#include"opencv2\imgproc\imgproc.hpp"
#include"opencv2\opencv.hpp"
#include <omp.h>  /*多线程计算库*/


#define InSAR_API __declspec(dllexport)
#ifdef _DEBUG

#pragma comment(lib, "opencv_world450d.lib")


#else

#pragma comment(lib, "opencv_world450.lib")


#endif // DEBUG


#endif // !__PACKAGE__H__

