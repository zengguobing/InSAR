#pragma once
#ifndef __PACKAGE__H__
#define __PACKAGE__H__
#define PI 3.141592653589793238
#define VEL_C 299792458.0
#define INPUTMAXSIZE 1024
#include"opencv2\core\core.hpp"
#include"opencv2\highgui\highgui.hpp"
#include"opencv2\imgproc\imgproc.hpp"
#include"opencv2\opencv.hpp"
#include <omp.h>  /*多线程计算库*/

/*-------------------------------------------------------*/
/*                    三维位置矢量                       */
/*-------------------------------------------------------*/
struct Position
{
	double x;
	double y;
	double z;
	Position()
	{
		this->x = 0.0;
		this->y = 0.0;
		this->z = 0.0;
	}
	/*传值构造函数*/
	Position(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	/*拷贝构造函数*/
	Position(const Position& cp)
	{
		this->x = cp.x;
		this->y = cp.y;
		this->z = cp.z;
	}
	/*赋值函数(深拷贝)*/
	Position operator=(const Position& cp)
	{
		this->x = cp.x;
		this->y = cp.y;
		this->z = cp.z;
		return *this;
	}

};

/*-------------------------------------------------------*/
/*                    三维速度矢量                       */
/*-------------------------------------------------------*/
struct Velocity
{
	double vx;
	double vy;
	double vz;
	Velocity()
	{
		this->vx = 0.0;
		this->vy = 0.0;
		this->vz = 0.0;
	}
	/*传值构造函数*/
	Velocity(double vx, double vy, double vz)
	{
		this->vx = vx;
		this->vy = vy;
		this->vz = vz;
	}
	/*拷贝构造函数*/
	Velocity(const Velocity& cp)
	{
		this->vx = cp.vx;
		this->vy = cp.vy;
		this->vz = cp.vz;
	}
	/*赋值函数(深拷贝)*/
	Velocity operator=(const Velocity& cp)
	{
		this->vx = cp.vx;
		this->vy = cp.vy;
		this->vz = cp.vz;
		return *this;
	}

};

/*-------------------------------------------------------*/
/*                   卫星轨道信息                        */
/*-------------------------------------------------------*/
struct OSV
{
	double time;
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	OSV()
	{
		time = x = y = z = vx = vy = vz = 0.0;
	}
	OSV(double time, double x, double y, double z, double vx, double vy, double vz)
	{
		this->time = time;
		this->x = x;
		this->y = y;
		this->z = z;
		this->vx = vx;
		this->vy = vy;
		this->vz = vz;
	}
	/*拷贝构造函数*/
	OSV(const OSV& osv)
	{
		this->time = osv.time;
		this->x = osv.x;
		this->y = osv.y;
		this->z = osv.z;
		this->vx = osv.vx;
		this->vy = osv.vy;
		this->vz = osv.vz;
	}
	/*赋值函数*/
	OSV operator=(const OSV& osv)
	{
		this->time = osv.time;
		this->x = osv.x;
		this->y = osv.y;
		this->z = osv.z;
		this->vx = osv.vx;
		this->vy = osv.vy;
		this->vz = osv.vz;
		return *this;
	}

};

#define InSAR_API __declspec(dllexport)
#ifdef _DEBUG

#pragma comment(lib, "opencv_world450d.lib")


#else

#pragma comment(lib, "opencv_world450.lib")


#endif // DEBUG


#endif // !__PACKAGE__H__

