#include <math.h>
#include <string.h>
#include "image.h"
#include "IL/il.h"
#include "IL/ilu.h"
#include "IL/ilut.h"


template class LImage<unsigned char>;
template class LImage<unsigned short>;
template class LImage<unsigned int>;
template class LImage<double>;
template class LImage<int>;
template class LColourImage<unsigned char>;
template class LColourImage<double>;

template <class T>
LImage<T>::LImage()
{
	data = NULL;
	width = height = bands = 0;
}

template <class T>
LImage<T>::LImage(int setBands)
{
	data = NULL;
	width = height = 0;
	bands = setBands;
}

template <class T>
LImage<T>::LImage(int setWidth, int setHeight, int setBands)
{
	data = new T[setWidth * setHeight * setBands];
	width = setWidth, height = setHeight, bands = setBands;
}

template <class T>
LImage<T>::LImage(LImage<T> &image)
{
	width = image.GetWidth();
	height = image.GetHeight();
	bands = image.GetBands();
	data = new T[width * height * bands];
	memcpy(data, image.GetData(), width * height * bands * sizeof(T));
}

template <class T>
LImage<T>::LImage(char *fileName)
{
	data = NULL;
	Load(fileName);
}

template <class T>
LImage<T>::~LImage()
{
	if(data != NULL) delete[] data;
}

template <class T>
void LImage<T>::SetResolution(int setWidth, int setHeight)
{
	if(data != NULL) delete[] data;
	data = new T[setWidth * setHeight * bands];
	width = setWidth, height = setHeight;
}

template <class T>
void LImage<T>::SetResolution(int setWidth, int setHeight, int setBands)
{
	if(data != NULL) delete[] data;
	data = new T[setWidth * setHeight * setBands];
	width = setWidth, height = setHeight, bands = setBands;
}

template <class T>
void LImage<T>::CopyDataFrom(T *fromData)
{
	memcpy(data, fromData, width * height * bands * sizeof(T));
}

template <class T>
void LImage<T>::CopyDataTo(T *toData)
{
	memcpy(toData, data, width * height * bands * sizeof(T));
}

template <class T>
void LImage<T>::CopyDataFrom(LImage<T> &image)
{
	if(data != NULL) delete[] data;
	width = image.GetWidth();
	height = image.GetHeight();
	bands = image.GetBands();
	data = new T[width * height * bands];
	memcpy(data, image.GetData(), width * height * bands * sizeof(T));
}

template <class T>
T *LImage<T>::GetData()
{
	return(data);
}

template <class T>
T &LImage<T>::GetValue(int index)
{
	return(data[index]);
}

template <class T>
T *LImage<T>::operator()(int x, int y)
{
	return(&data[(y * width + x) * bands]);
}

template <class T>
T &LImage<T>::operator()(int x, int y, int band)
{
	return(data[(y * width + x) * bands + band]);
}

template <class T>
int LImage<T>::GetWidth()
{
	return(width);
}

template <class T>
int LImage<T>::GetHeight()
{
	return(height);
}

template <class T>
int LImage<T>::GetBands()
{
	return(bands);
}

template <class T>
int LImage<T>::GetPoints()
{
	return(width * height);
}

template <class T>
int LImage<T>::GetSize()
{
	return(width * height * bands);
}

template <class T>
void LImage<T>::Save(char *fileName)
{
	FILE *f;
	f = fopen(fileName, "wb");
	fwrite(&width, 4, 1, f);
	fwrite(&height, 4, 1, f);
	fwrite(&bands, 4, 1, f);
	fwrite(data, sizeof(T), width * height * bands, f);

	fclose(f);
}

template <class T>
void LImage<T>::Load(char *fileName)
{
	if(data != NULL) delete[] data;

	FILE *f;
	f = fopen(fileName, "rb");
	if(f == NULL) _error(fileName);

	fread(&width, 4, 1, f);
	fread(&height, 4, 1, f);
	fread(&bands, 4, 1, f);

	data = new T[width * height * bands];
	fread(data, sizeof(T), width * height * bands, f);
	fclose(f);
}

template <class T>
int LImage<T>::Exist(char *fileName)
{
	FILE *f;
	f = fopen(fileName, "rb");
	if(f == NULL) return(0);
	else
	{
   		fclose(f);
		return(1);
	}
}

template <class T>
LColourImage<T>::LColourImage() : LImage<T>()
{
}

template <class T>
LColourImage<T>::LColourImage(int setBands) : LImage<T>(setBands)
{
}

template <class T>
LColourImage<T>::LColourImage(int setWidth, int setHeight, int setBands) : LImage<T>(setWidth, setHeight, setBands)
{
}

template <class T>
LColourImage<T>::LColourImage(LColourImage<T> &image) : LImage<T>(image)
{
}


LRgbImage::LRgbImage() : LColourImage<unsigned char>(3)
{
}

LRgbImage::LRgbImage(int setWidth, int setHeight) : LColourImage<unsigned char>(setWidth, setHeight, 3)
{
}

LRgbImage::LRgbImage(char *fileName) : LColourImage<unsigned char>(3)
{
	Load(fileName);
}

LRgbImage::LRgbImage(LRgbImage &rgbImage) : LColourImage<unsigned char>(rgbImage)
{
}


void LRgbImage::Load(char *fileName)
{
	
	ILuint ilImage;

	ilGenImages(1, &ilImage);
	ilBindImage(ilImage);

    ilEnable(IL_ORIGIN_SET); 
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT); 

	ilLoadImage(fileName);
	ilConvertImage(IL_RGB, IL_UNSIGNED_BYTE);

	SetResolution(ilGetInteger(IL_IMAGE_WIDTH), ilGetInteger(IL_IMAGE_HEIGHT));
	CopyDataFrom(ilGetData());

	ilDeleteImages(1, &ilImage);
}

void LRgbImage::Load(LRgbImage &rgbImage)
{
	rgbImage.Save(*this);
}

void LRgbImage::Save(LRgbImage &rgbImage)
{
	rgbImage.SetResolution(width, height);
	CopyDataTo(rgbImage.GetData());
}

void LRgbImage::Save(char *fileName)
{
	ILuint ilImage;

	ilGenImages(1, &ilImage);
	ilBindImage(ilImage);

	ilTexImage(width, height, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);
	CopyDataTo(ilGetData());

	ilEnable(IL_FILE_OVERWRITE);
	ilSaveImage(fileName);
	ilDeleteImages(1, &ilImage);
}


