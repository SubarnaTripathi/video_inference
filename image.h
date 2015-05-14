#ifndef __image
#define __image

#include <cstring>
#include "std.h"

class LRgbImage;
class LGreyImage;

template <class T>
class LImage
{
	protected :
		T *data;
		int width, height, bands;
	public :
		LImage();
		LImage(int setBands);
		LImage(int setWidth, int setHeight, int setBands);
		LImage(LImage<T> &image);
		LImage(char *fileName);
		~LImage();

		void SetResolution(int setWidth, int setHeight);
		void SetResolution(int setWidth, int setHeight, int setBands);
		void CopyDataFrom(LImage<T> &image);
		void CopyDataFrom(T *fromData);
		void CopyDataTo(T *toData);

		T *GetData();
		T *operator()(int x, int y);
		T &operator()(int x, int y, int band);
		T &GetValue(int index);

		int GetWidth();
		int GetHeight();
		int GetPoints();
		int GetBands();
		int GetSize();

		virtual void Save(char *fileName);
		virtual void Load(char *fileName);
		static int Exist(char *fileName);
};

template <class T>
class LColourImage : public LImage<T>
{
	public :
		LColourImage();
		LColourImage(int setBands);
		LColourImage(int setWidth, int setHeight, int setBands);
		LColourImage(LColourImage<T> &image);

};

class LRgbImage : public LColourImage<unsigned char>
{
	public :
		
		LRgbImage();
		LRgbImage(int setWidth, int setHeight);
		LRgbImage(char *fileName);
		LRgbImage(LRgbImage &rgbImage);
		LRgbImage(LGreyImage &greyImage);
		void Load(char *fileName);
		void Load(LRgbImage &rgbImage);
		void Load(LGreyImage &greyImage);
		void Save(char *fileName);
		void Save(LRgbImage &rgbImage);
		void Save(LGreyImage &greyImage);
};



//#include "filter.h"
//#include "std.h"
//
//class LRgbImage;
//class LGreyImage;
//class LLuvImage;
//class LLabImage;
//class LLabelImage;
//class LCostImage;
//class LSegmentImage;
//class LDataset;
//class LCrfDomain;
//template <class T> class LFilter2D;
//
//template <class T>
//class LImage
//{
//	protected :
//		T *data;
//		int width, height, bands;
//	public :
//		LImage();
//		LImage(int setBands);
//		LImage(int setWidth, int setHeight, int setBands);
//		LImage(LImage<T> &image);
//		LImage(char *fileName);
//		~LImage();
//
//		void SetResolution(int setWidth, int setHeight);
//		void SetResolution(int setWidth, int setHeight, int setBands);
//		void CopyDataFrom(LImage<T> &image);
//		void CopyDataFrom(T *fromData);
//		void CopyDataTo(T *toData);
//
//		T *GetData();
//		T *operator()(int x, int y);
//		T &operator()(int x, int y, int band);
//		T &GetValue(int index);
//
//		int GetWidth();
//		int GetHeight();
//		int GetPoints();
//		int GetBands();
//		int GetSize();
//
//		void WindowTo(LImage &imageTo, int centreX, int centreY, int halfWidth, int halfHeight);
//		void WindowFrom(LImage &imageTo, int centreX, int centreY, int halfWidth, int halfHeight);
//
//		virtual void Save(char *fileName);
//		virtual void Load(char *fileName);
//		static int Exist(char *fileName);
//};
//
//template <class T>
//class LColourImage : public LImage<T>
//{
//	public :
//		LColourImage();
//		LColourImage(int setBands);
//		LColourImage(int setWidth, int setHeight, int setBands);
//		LColourImage(LColourImage<T> &image);
//
//		void FilterFrom(LColourImage<T> &fromImage, LFilter2D<T> &filter);
//		void FilterTo(LColourImage<T> &toImage, LFilter2D<T> &filter);
//		void ScaleTo(LColourImage &imageTo, double scale);
//		void ScaleFrom(LColourImage &imageFrom, double scale);
//		void ScaleTo(LColourImage &imageTo, double centreX, double centreY, double scale, int halfWidth, int halfHeight);
//		void ScaleFrom(LColourImage &imageFrom, double centreX, double centreY, double scale, int halfWidth, int halfHeight);
//		void RotateTo(LColourImage &imageTo, double centreX, double centreY, double angle, int halfWidth, int halfHeight);
//		void RotateFrom(LColourImage &imageFrom, double centreX, double centreY, double angle, int halfWidth, int halfHeight);
//		void ScaleRotateTo(LColourImage &imageTo, double centreX, double centreY, double scale, double angle, int halfWidth, int halfHeight);
//		void ScaleRotateFrom(LColourImage &imageFrom, double centreX, double centreY, double scale, double angle, int halfWidth, int halfHeight);
//		void AffineTo(LColourImage &imageTo, double centreX, double centreY, double scale, double *u, int halfWidth, int halfHeight);
//		void AffineFrom(LColourImage &imageFrom, double centreX, double centreY, double scale, double *u, int halfWidth, int halfHeight);
//		void AffineRotateTo(LColourImage &imageTo, double centreX, double centreY, double scale, double *u, double angle, int halfWidth, int halfHeight);
//		void AffineRotateFrom(LColourImage &imageFrom, double centreX, double centreY, double scale, double *u, double angle, int halfWidth, int halfHeight);
//};
//
//class LRgbImage : public LColourImage<unsigned char>
//{
//	public :
//		static void RgbToGrey(unsigned char *rgb, double *grey);
//		static void RgbToLuv(unsigned char *rgb, double *luv);
//		static void RgbToLab(unsigned char *rgb, double *lab);
//
//		LRgbImage();
//		LRgbImage(int setWidth, int setHeight);
//		LRgbImage(char *fileName);
//		LRgbImage(LRgbImage &rgbImage);
//		LRgbImage(LGreyImage &greyImage);
//		LRgbImage(LLuvImage &luvImage);
//		LRgbImage(LLabImage &labImage);
//		LRgbImage(LLabelImage &labelImage, LCrfDomain *domain);
//		LRgbImage(LLabelImage &labelImage, LDataset *dataset, void (LDataset::*newLabelToRgb)(unsigned char *, unsigned char *));
//
//		LRgbImage(LSegmentImage &segmentImage, LRgbImage &rgbImage, int showBoundaries);
//		LRgbImage(LCostImage &costImage, LCrfDomain *domain, int showMaximum);
//
//		void Load(char *fileName);
//		void Load(LRgbImage &rgbImage);
//		void Load(LGreyImage &greyImage);
//		void Load(LLuvImage &luvImage);
//		void Load(LLabImage &labImage);
//		void Load(LLabelImage &labelImage, LCrfDomain *domain);
//		void Load(LLabelImage &labelImage, LDataset *dataset, void (LDataset::*newLabelToRgb)(unsigned char *, unsigned char *));
//		void Load(LSegmentImage &segmentImage, LRgbImage &rgbImage, int showBoundaries);
//		void Load(LCostImage &costImage, LCrfDomain *domain, int showMaximum);
//
//		void Save(char *fileName);
//		void Save(LRgbImage &rgbImage);
//		void Save(LGreyImage &greyImage);
//		void Save(LLuvImage &luvImage);
//		void Save(LLabImage &labImage);
//		void Save(LLabelImage &labelImage, LCrfDomain *domain);
//		void Save(LLabelImage &labelImage, LDataset *dataset, void (LDataset::*newRgbToLabel)(unsigned char *, unsigned char *));
//};
//
//class LGreyImage : public LColourImage<double>
//{
//	public :
//		static void GreyToRgb(double *grey, unsigned char *rgb);
//		static void GreyToLuv(double *grey, double *luv);
//		static void GreyToLab(double *grey, double *lab);
//
//		LGreyImage();
//		LGreyImage(int setWidth, int setHeight);
//		LGreyImage(char *fileName);
//		LGreyImage(LRgbImage &rgbImage);
//		LGreyImage(LGreyImage &greyImage);
//		LGreyImage(LLuvImage &luvImage);
//		LGreyImage(LLabImage &labImage);
//
//		void Load(char *fileName);
//		void Load(LRgbImage &rgbImage);
//		void Load(LGreyImage &greyImage);
//		void Load(LLuvImage &luvImage);
//		void Load(LLabImage &labImage);
//
//		void Save(char *fileName);
//		void Save(LRgbImage &rgbImage);
//		void Save(LGreyImage &greyImage);
//		void Save(LLuvImage &luvImage);
//		void Save(LLabImage &labImage);
//};

//class LLuvImage : public LColourImage<double>
//{
//	public :
//		static void LuvToRgb(double *luv, unsigned char *rgb);
//		static void LuvToGrey(double *luv, double *grey);
//
//		LLuvImage();
//		LLuvImage(int setWidth, int setHeight);
//		LLuvImage(char *fileName);
//		LLuvImage(LRgbImage &rgbImage);
//		LLuvImage(LGreyImage &greyImage);
//		LLuvImage(LLuvImage &luvImage);
//		LLuvImage(LLabImage &labImage);
//
//		void Load(char *fileName);
//		void Load(LRgbImage &rgbImage);
//		void Load(LGreyImage &greyImage);
//		void Load(LLuvImage &luvImage);
//		void Load(LLabImage &labImage);
//
//		void Save(char *fileName);
//		void Save(LRgbImage &rgbImage);
//		void Save(LGreyImage &greyImage);
//		void Save(LLuvImage &luvImage);
//		void Save(LLabImage &labImage);
//};
//
//class LLabImage : public LColourImage<double>
//{
//	public :
//		static void LabToGrey(double *luv, double *grey);
//		static void LabToRgb(double *lab, unsigned char *rgb);
//
//		LLabImage();
//		LLabImage(int setWidth, int setHeight);
//		LLabImage(char *fileName);
//		LLabImage(LRgbImage &rgbImage);
//		LLabImage(LGreyImage &greyImage);
//		LLabImage(LLabImage &labImage);
//		LLabImage(LLuvImage &luvImage);
//
//		void Load(char *fileName);
//		void Load(LRgbImage &rgbImage);
//		void Load(LGreyImage &greyImage);
//		void Load(LLuvImage &luvImage);
//		void Load(LLabImage &labImage);
//
//		void Save(char *fileName);
//		void Save(LRgbImage &rgbImage);
//		void Save(LGreyImage &greyImage);
//		void Save(LLuvImage &luvImage);
//		void Save(LLabImage &labImage);
//};
//
//class LLabelImage : public LImage<unsigned char>
//{
//	public :
//		LLabelImage();
//		LLabelImage(int setWidth, int setHeight);
//		LLabelImage(char *fileName, LCrfDomain *domain);
//		LLabelImage(char *fileName, LDataset *dataset, void (LDataset::*newRgbToLabel)(unsigned char *, unsigned char *));
//		LLabelImage(LRgbImage &rgbImage, LCrfDomain *domain);
//		LLabelImage(LRgbImage &rgbImage, LDataset *dataset, void (LDataset::*newRgbToLabel)(unsigned char *, unsigned char *));
//		LLabelImage(LLabelImage &labelImage);
//		LLabelImage(LCostImage &costImage, int showMaximum);
//
//		void Load(char *fileName, LCrfDomain *domain);
//		void Load(char *fileName, LDataset *dataset, void (LDataset::*newRgbToLabel)(unsigned char *, unsigned char *));
//		void Load(LRgbImage &rgbImage, LCrfDomain *domain);
//		void Load(LRgbImage &rgbImage, LDataset *dataset, void (LDataset::*newRgbToLabel)(unsigned char *, unsigned char *));
//		void Load(LLabelImage &labelImage);
//		void Load(LCostImage &costImage, int showMaximum);
//
//		void Save(char *fileName, LCrfDomain *domain);
//		void Save(char *fileName, LDataset *dataset, void (LDataset::*newLabelToRgb)(unsigned char *, unsigned char *));
//		void Save(char *fileName, LRgbImage &rgbImage, LCrfDomain *domain);
//		void Save(char *fileName, LRgbImage &rgbImage, LDataset *dataset, void (LDataset::*newLabelToRgb)(unsigned char *, unsigned char *));
//		void Save(LRgbImage &rgbImage, LCrfDomain *domain);
//		void Save(LRgbImage &rgbImage, LDataset *dataset, void (LDataset::*newLabelToRgb)(unsigned char *, unsigned char *));
//
//		void Save(LLabelImage &labelImage);
//		void Save8bit(char *fileName);
//};
//
//class LCostImage : public LImage<double>
//{
//	private :
//		void CostToLabel(double *cost, unsigned char *label, int showMaximum);
//	public :
//		LCostImage();
//		LCostImage(LCostImage &costImage);
//		LCostImage(int setWidth, int setHeight, int setBands);
//
//		void Save(char *fileName, LCrfDomain *domain, int showMaximum);
//		void Save(LRgbImage &rgbImage, LCrfDomain *domain, int showMaximum);
//		void Save(LLabelImage &labelImage, int showMaximum);
//};
//
//class LSegmentImage : public LImage<int>
//{
//	public :
//		LSegmentImage();
//		LSegmentImage(LSegmentImage &segmentImage);
//		LSegmentImage(int setWidth, int setHeight);
//
//		void Save(char *fileName, LRgbImage &rgbImage, int showBoundaries);
//		void Save(LRgbImage &segmentRgbImage, LRgbImage &rgbImage, int showBoundaries);
//};
//
//class LIntegralImage
//{
//	protected :
//	public :
//		int width, height;
//		virtual ~LIntegralImage() {};
//	public :
//		LIntegralImage();
//		virtual int Response(int x1, int y1, int x2, int y2) { return(0); };
//		virtual double DResponse(double x1, double y1, double x2, double y2) { return(0); };
//		virtual void Load(LImage<unsigned short> &dataImage, int subSample, int index) {};
//		virtual void Copy(LImage<double> &dataImage, int subSample, int index, double scale) {};
//};
//
//class LIntegralImage4B : public LIntegralImage
//{
//	private :
//		LImage <unsigned int> image;
//	public :
//		int Response(int x1, int y1, int x2, int y2);
//		double DResponse(double x1, double y1, double x2, double y2);
//		void Load(LImage<unsigned short> &dataImage, int subSample, int index);
//		void Copy(LImage<double> &dataImage, int subSample, int index, double scale);
//};
//
//class LIntegralImage2B : public LIntegralImage
//{
//	private :
//		LImage <unsigned short> image;
//	public :
//		int Response(int x1, int y1, int x2, int y2);
//		double DResponse(double x1, double y1, double x2, double y2);
//		void Load(LImage<unsigned short> &dataImage, int subSample, int index);
//		void Copy(LImage<double> &dataImage, int subSample, int index, double scale);
//};
//
//class LIntegralImage1B : public LIntegralImage
//{
//	private :
//		LImage <unsigned char> image;
//	public :
//		int Response(int x1, int y1, int x2, int y2);
//		double DResponse(double x1, double y1, double x2, double y2);
//		void Load(LImage<unsigned short> &dataImage, int subSample, int index);
//		void Copy(LImage<double> &dataImage, int subSample, int index, double scale);
//};
//
//class LIntegralImageHB : public LIntegralImage
//{
//	private :
//		LImage <unsigned char> image;
//	public :
//		int Response(int x1, int x2, int y1, int y2);
//		double DResponse(double x1, double y1, double x2, double y2);
//		void Load(LImage<unsigned short> &dataImage, int subSample, int index);
//};

#endif