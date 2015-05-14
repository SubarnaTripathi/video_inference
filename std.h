#ifndef __std
#define __std

#include <conio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _WIN64
#pragma comment(lib,"lib64\\DevIL.LIB")
#pragma comment(lib,"lib64\\ILU.LIB")
#pragma comment(lib,"lib64\\ILUT.LIB")
#else
#pragma comment(lib,"lib32\\DevIL.LIB")
#pragma comment(lib,"lib32\\ILU.LIB")
#pragma comment(lib,"lib32\\ILUT.LIB")
#endif

//#define MULTITHREAD
#define MAXTHREAD 16

namespace LMath
{
	static const double pi = (double)3.1415926535897932384626433832795;
	static const double positiveInfinity = (double)1e50;
	static const double negativeInfinity = (double)-1e50;
	static const double almostZero = (double)1e-12;
};

template <class T>
class LList
{
	private :
		int count, capacity;
		T *items;

		void Resize(int size);
	public :
		LList();
		~LList();
		
		T &operator[](int index);
		T *GetArray();
};

void _error(char *str);


//#include <conio.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//
//#ifdef _WIN64
//#pragma comment(lib,"lib64\\DevIL.LIB")
//#pragma comment(lib,"lib64\\ILU.LIB")
//#pragma comment(lib,"lib64\\ILUT.LIB")
//#else
//#pragma comment(lib,"lib32\\DevIL.LIB")
//#pragma comment(lib,"lib32\\ILU.LIB")
//#pragma comment(lib,"lib32\\ILUT.LIB")
//#endif
//
////#define MULTITHREAD
//#define MAXTHREAD 16
//
//namespace LMath
//{
//	static const double pi = (double)3.1415926535897932384626433832795;
//	static const double positiveInfinity = (double)1e50;
//	static const double negativeInfinity = (double)-1e50;
//	static const double almostZero = (double)1e-12;
//
//	void SetSeed(unsigned int seed);
//	unsigned int RandomInt();
//	unsigned int RandomInt(unsigned int maxval);
//	unsigned int RandomInt(unsigned int minval, unsigned int maxval);
//	double RandomReal();
//	double RandomGaussian(double mi, double var);
//
//	double SquareEuclidianDistance(double *v1, double *v2, int size);
//	double KLDivergence(double *histogram, double *referenceHistogram, int size, double threshold);
//	double GetAngle(double x, double y);
//};
//
//template <class T>
//class LList
//{
//	private :
//		int count, capacity;
//		T *items;
//
//		void Resize(int size);
//		void QuickSort(int from, int to, int (*sort)(T, T));
//	public :
//		LList();
//		~LList();
//		
//		T &operator[](int index);
//		T &Add(T value);
//		T &Insert(T value, int index);
//		void Delete(int index);
//		void Swap(int index1, int index2);
//		void Sort(int (*sort)(T, T));
//		T *GetArray();
//		int GetCount();
//		void Clear();
//};
//
//struct LBoundingBox
//{
//	unsigned char type, pose;
//	int x1, y1, x2, y2;
//};
//
//struct LResponse
//{
//	double val;
//	unsigned char correct;
//};
//
//void _error(char *str);
//
//char *GetFileName(const char *folder, const char *name, const char *extension);
//void ForceDirectory(const char *dir);
//void ForceDirectory(const char *dir, const char *subdir);

#endif