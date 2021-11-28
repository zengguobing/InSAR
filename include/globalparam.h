#ifndef GDALBRIDGE_H_INCLUDED
#define GDALBRIDGE_H_INCLUDED

/* -------------------------------------------------------------------- */
/*      Start C context.                                                */
/* -------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
    
/* ==================================================================== */
/*      Standard types and defines normally supplied by cpl_port.h.     */
/* ==================================================================== */
#if UINT_MAX == 65535
typedef long            GInt32;
typedef unsigned long   GUInt32;
#else
typedef int             GInt32;
typedef unsigned int    GUInt32;
#endif

typedef short           GInt16;
typedef unsigned short  GUInt16;
typedef unsigned char   GByte;
typedef int             GBool;
typedef long long        GIntBig;
typedef unsigned long long GUIntBig;

#ifndef FALSE
#define FALSE		0
#define TRUE		1
#endif

#ifndef NULL
#  define NULL		0
#endif

#ifndef GDAL_ENTRY
#  define GDAL_ENTRY extern
#  define GDAL_NULL
#endif

typedef unsigned long long uint64_t; 

typedef struct {
	long tv_sec;
	long tv_usec;
} ModTime;
/* -------------------------------------------------------------------- */
/*      Significant constants.                                          */
/* -------------------------------------------------------------------- */

/*! Pixel data types */
typedef enum {
    /*! Unknown or unspecified type */ 		    GDT_Unknown = 0,
    /*! Eight bit unsigned integer */ 		    GDT_Byte = 1,
    /*! Sixteen bit unsigned integer */         GDT_UInt16 = 2,
    /*! Sixteen bit signed integer */           GDT_Int16 = 3,
    /*! Thirty two bit unsigned integer */      GDT_UInt32 = 4,
    /*! Thirty two bit signed integer */        GDT_Int32 = 5,
    /*! Thirty two bit floating point */        GDT_Float32 = 6,
    /*! Sixty four bit floating point */        GDT_Float64 = 7,
    /*! Complex Int16 */                        GDT_CInt16 = 8,
    /*! Complex Int32 */                        GDT_CInt32 = 9,
    /*! Complex Float32 */                      GDT_CFloat32 = 10,
    /*! Complex Float64 */                      GDT_CFloat64 = 11,
    GDT_TypeCount = 12		/* maximum type # + 1 */
} GDALDataType;

typedef enum {
    GA_ReadOnly = 0,
    GA_Update = 1
} GDALAccess;

typedef enum {
    GF_Read = 0,
    GF_Write = 1
} GDALRWFlag;

/*! Types of color interpretation for raster bands. */
typedef enum
{
    GCI_Undefined=0,
    /*! Greyscale */                                      GCI_GrayIndex=1,
    /*! Paletted (see associated color table) */          GCI_PaletteIndex=2,
    /*! Red band of RGBA image */                         GCI_RedBand=3,
    /*! Green band of RGBA image */                       GCI_GreenBand=4,
    /*! Blue band of RGBA image */                        GCI_BlueBand=5,
    /*! Alpha (0=transparent, 255=opaque) */              GCI_AlphaBand=6,
    /*! Hue band of HLS image */                          GCI_HueBand=7,
    /*! Saturation band of HLS image */                   GCI_SaturationBand=8,
    /*! Lightness band of HLS image */                    GCI_LightnessBand=9,
    /*! Cyan band of CMYK image */                        GCI_CyanBand=10,
    /*! Magenta band of CMYK image */                     GCI_MagentaBand=11,
    /*! Yellow band of CMYK image */                      GCI_YellowBand=12,
    /*! Black band of CMLY image */                       GCI_BlackBand=13
} GDALColorInterp;

/*! Types of color interpretations for a GDALColorTable. */
typedef enum 
{
  /*! Grayscale (in GDALColorEntry.c1) */                      GPI_Gray=0,
  /*! Red, Green, Blue and Alpha in (in c1, c2, c3 and c4) */  GPI_RGB=1,
  /*! Cyan, Magenta, Yellow and Black (in c1, c2, c3 and c4)*/ GPI_CMYK=2,
  /*! Hue, Lightness and Saturation (in c1, c2, and c3) */     GPI_HLS=3
} GDALPaletteInterp;

/* -------------------------------------------------------------------- */
/*      GDAL Specific error codes.                                      */
/*                                                                      */
/*      error codes 100 to 299 reserved for GDAL.                       */
/* -------------------------------------------------------------------- */
typedef enum
{
    CE_None = 0,
    CE_Log = 1,
    CE_Warning = 2,
    CE_Failure = 3,
    CE_Fatal = 4
  
} CPLErr;

#define CPLE_AppDefined			1
#define CPLE_OutOfMemory		2
#define CPLE_FileIO			3
#define CPLE_OpenFailed			4
#define CPLE_IllegalArg			5
#define CPLE_NotSupported		6
#define CPLE_AssertionFailed		7
#define CPLE_NoWriteAccess		8

#define CPLE_WrongFormat	200

typedef int OGRErr;




/* ==================================================================== */
/*      GDAL_GCP                                                        */
/* ==================================================================== */

/** Ground Control Point */
typedef struct
{
    /** Unique identifier, often numeric */
    char	*pszId; 

    /** Informational message or "" */
    char	*pszInfo;

    /** Pixel (x) location of GCP on raster */
    double 	dfGCPPixel;
    /** Line (y) location of GCP on raster */
    double	dfGCPLine;

    /** X position of GCP in georeferenced space */
    double	dfGCPX;

    /** Y position of GCP in georeferenced space */
    double	dfGCPY;

    /** Elevation of GCP, or zero if not known */
    double	dfGCPZ;
} GDAL_GCP;


/* ==================================================================== */
/*      Color tables.                                                   */
/* ==================================================================== */
/** Color tuple */
typedef struct
{
    /*! gray, red, cyan or hue */
    short      c1;      

    /*! green, magenta, or lightness */    
    short      c2;      

    /*! blue, yellow, or saturation */
    short      c3;      

    /*! alpha or blackband */
    short      c4;      
} GDALColorEntry;

struct Heap //创建小顶堆以储存解缠邻接点序列
{
	int size;
	int* x = (int*)malloc(sizeof(int) * 1000000);
	int* y = (int*)malloc(sizeof(int) * 1000000);
	double* queue = (double*)malloc(sizeof(double) * 1000000);
public:
	Heap()         //初始化 
	{
		size = 0;
		/*for (int i = 0; i < 1000000; i++)
		{
			queue[i] = 0;
			x[i] = 0;
			y[i] = 0;
		}*/
	}
	~Heap()
	{
	}

	void shift_up(int i)  //上浮 
	{
		while (i > 1)
		{
			if (queue[i] < queue[i >> 1])
			{
				double temp = queue[i];
				queue[i] = queue[i >> 1];
				queue[i >> 1] = temp;
				int temp2 = x[i];
				x[i] = x[i >> 1];
				x[i >> 1] = temp2;
				temp2 = y[i];
				y[i] = y[i >> 1];
				y[i >> 1] = temp2;
			}
			i >>= 1;
		}
	}
	void shift_down(int i)   //下沉 
	{
		while ((i << 1) <= size)
		{
			int next = i << 1;
			if (next < size && queue[next + 1] < queue[next])
				next++;
			if (queue[i] > queue[next])
			{
				double temp = queue[i];
				queue[i] = queue[next];
				queue[next] = temp;
				int temp2 = x[i];
				x[i] = x[next];
				x[next] = temp2;
				temp2 = y[i];
				y[i] = y[next];
				y[next] = temp2;
				i = next;
			}
			else return;
		}
	}
	int push(double v, int map_x, int map_y)   //加入元素 
	{
		if (v < 0 ||
			map_x < 0 ||
			map_y < 0)
		{
			fprintf(stderr, "Heap::push(): input check failed!\n\n");
			return -1;
		}
		queue[++size] = v;
		x[size] = map_x;
		y[size] = map_y;
		shift_up(size);
		return 0;
	}
	int pop()         //弹出操作 
	{
		if (size <= 0)
		{
			fprintf(stderr, "Heap::pop(): no data!\n\n");
			return -1;
		}
		double temp = queue[1];
		queue[1] = queue[size];
		queue[size] = temp;
		int temp2 = x[1];
		x[1] = x[size];
		x[size] = temp2;
		temp2 = y[1];
		y[1] = y[size];
		y[size] = temp2;
		size--;
		shift_down(1);
		return 0;
	}
	int top(int* x0, int* y0)			//输出根节点坐标
	{
		if (size <= 0)
		{
			fprintf(stderr, "Heap::top(): no data!\n\n");
			return -1;
		}
		*x0 = x[1];
		*y0 = y[1];
		return 1;
	}
	bool empty()
	{
		return size;
	}
};




/* -------------------------------------------------------------------- */
/*      Terminate C context.                                            */
/* -------------------------------------------------------------------- */
#ifdef __cplusplus
}
#endif



#endif /* ndef GDALPARAM_H_INCLUDED */