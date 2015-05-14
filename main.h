
//#pragma warning(disable : 4305)

#include <fstream>
#include <vector>
#include <io.h>
#include <time.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <Windows.h>
#include "image.h"
#include "IL/il.h"
#include "probimage.h"

//subarna
# define ALE_UNARY_CAMVID
//# define SUPERVOXEL
//# define MATLAB_DEMO

# define MIMAGE
# define INTER

# ifndef MIMAGE
# undef INTER
# endif

# ifdef MIMAGE
 void mem_initM(char **TrainFileName, int num_images);
 void loadfilesM(char **TrainFileName, int num_images);
# endif

void mem_init(char *TrainFileName);
void del_meminit();
void loadfiles(char *);
float *dataCost, *cooc_unary, *cooc_pairwise; 
int imagewidth, imageheight, K; 
unsigned char *final_labels; 
char *fileNameL;

# ifndef ALE_UNARY_CAMVID  //subarna
int num_of_labels = 21; 
# else
int num_of_labels = 11; 
# endif

int num_files; 
void labeltorgb(short *, char *); 

//subarna
void labeltorgbM(int, short *, char **); 

# ifndef MATLAB_DEMO
char **TrainFileNames; 
void read_files(); 
# else
char *FileName;
# endif


unsigned char *im_orig;

// ho_terms
int num_of_layers; 
char **layers_dir; 
char *ho_stats_pot; 
char *ho_seg_ext, *ho_sts_ext; 

// det_terms
char *det_seg_dir; 
char *det_bb_dir; 
char *det_seg_ext; 
char *det_bb_ext; 



void del_meminit()
{
	delete []dataCost; 

	delete []final_labels; 
	delete []cooc_unary; 
	delete []cooc_pairwise; 
}

# ifdef MIMAGE //subarna
void mem_initM(char **trainfileName, int num_images)
{
	fileNameL = new char[1024]; 
	sprintf(fileNameL, "Pascal\\Result\\Image\\%s.jpg", trainfileName[0]);

# ifdef ALE_UNARY_CAMVID
	sprintf(fileNameL, "CamVid\\Result\\CImage\\%s.jpg", trainfileName[0]);
# endif

	LRgbImage rgbImage(fileNameL); 
	imageheight = rgbImage.GetHeight(); imagewidth = rgbImage.GetWidth(); 

	dataCost = new float[num_images*imagewidth*imageheight*num_of_labels]; 
	final_labels = new unsigned char[imagewidth*imageheight];
	cooc_unary = new float[num_of_labels]; 
	cooc_pairwise = new float[num_of_labels*num_of_labels];
	memset(dataCost, 0, num_images*imagewidth*imageheight*num_of_labels*sizeof(float));
	memset(cooc_unary, 0, num_of_labels*sizeof(float));
	memset(cooc_pairwise, 0, num_of_labels*num_of_labels*sizeof(float));
	memset(final_labels, 0, imagewidth*imageheight*sizeof(unsigned char));
	
	loadfilesM(trainfileName, num_images); 

	im_orig = new unsigned char[num_images*3*imagewidth*imageheight];

	for (int n = 0; n < num_images; n++ )
	{
		sprintf(fileNameL, "CamVid\\Result\\CImage\\%s.jpg", trainfileName[n]);
		LRgbImage rgbImage(fileNameL); 
		unsigned char *im = rgbImage.GetData();
		long offset = n*3*imageheight*imagewidth;

		for(int i = 0; i < imageheight; i++)
		{
			for(int j = 0; j < imagewidth; j++)
			{
				im_orig[offset + 3*((imageheight-(i+1))*imagewidth+j)+0] = im[3*(i*imagewidth+j)+0];
				im_orig[offset + 3*((imageheight-(i+1))*imagewidth+j)+1] = im[3*(i*imagewidth+j)+1];
				im_orig[offset + 3*((imageheight-(i+1))*imagewidth+j)+2] = im[3*(i*imagewidth+j)+2];
			}
		}
	}

	delete []fileNameL;
}

void loadfilesM(char **trainfileName, int num_images)
{
	char *filename = new char[1024]; 
	// read unary potential
	ProbImage im_unary;

	// start ALE boosted per pixel unary potential	
	// Uncomment following lines if you wants to use unary potentials evaluated using ALE
    /************************************************************************************/	
	double *classifier_val_reverse = new double[imagewidth*imageheight*num_of_labels];
	double *classifier_val = new double[imagewidth*imageheight*num_of_labels];

	for(int i = 0; i < imagewidth*imageheight*num_of_labels; i++)
	{
		classifier_val_reverse[i] = 0.0; 
		classifier_val[i] = 0.0; 
	}

	for (int n = 0; n < num_images; n++ )
	{
		long offset = n*imageheight*imagewidth;
		//char *filename = new char[200]; 
		sprintf(filename, "CamVid\\Result\\CDenseALE\\%s.dns", trainfileName[n]);
		FILE *fp = fopen(filename, "rb");
		int temp[3]; 
		fread(temp, sizeof(int), 3, fp);
		fread(classifier_val_reverse, sizeof(double), imagewidth*imageheight*num_of_labels, fp);
		fclose(fp);
	
		// start
		for(int k = 0; k < num_of_labels; k++)
		{
			for(int i = 0; i < imageheight; i++)
			{
				for(int j = 0; j < imagewidth; j++)
				{
					classifier_val[num_of_labels*((imageheight-(1+i))*imagewidth+j)+k] = classifier_val_reverse[num_of_labels*(i*imagewidth+j)+k]; 	
				}
			}
		}

		// calculating the per pixel dataCost
		double sum = 0.0; 
		for(int i = 0; i < imageheight; i++)
		{
			for(int j = 0; j < imagewidth; j++)
			{
				sum = 0.0; double data_cost_val = 0.0; 
				for(int l = 0; l < num_of_labels; l++) 
					sum += exp(classifier_val[num_of_labels*(i*imagewidth+j)+l]);

				for(int l = 0; l < num_of_labels; l++)  
				{
					data_cost_val = -1.0 * log(exp(classifier_val[num_of_labels*(i*imagewidth+j)+l]) /sum);
					dataCost[offset*num_of_labels + num_of_labels*(i*imagewidth+j)+l] = data_cost_val;	
				}
			}
		}
	}

	delete []classifier_val; 
	delete []classifier_val_reverse; 
    /********************************************************************************/
	// end ALE boosted per pixel unary potential

	delete []filename;
}

# endif

void mem_init(char *fileName)
{
	fileNameL = new char[200]; 
	sprintf(fileNameL, "Pascal\\Result\\Image\\%s.jpg", fileName);

# ifdef ALE_UNARY_CAMVID
	sprintf(fileNameL, "CamVid\\Result\\CImage\\%s.jpg", fileName);
# endif

	LRgbImage rgbImage(fileNameL); 
	imageheight = rgbImage.GetHeight(); imagewidth = rgbImage.GetWidth(); 

	dataCost = new float[imagewidth*imageheight*num_of_labels]; 
	final_labels = new unsigned char[imagewidth*imageheight];
	cooc_unary = new float[num_of_labels]; 
	cooc_pairwise = new float[num_of_labels*num_of_labels];
	memset(dataCost, 0, imagewidth*imageheight*num_of_labels*sizeof(float));
	memset(cooc_unary, 0, num_of_labels*sizeof(float));
	memset(cooc_pairwise, 0, num_of_labels*num_of_labels*sizeof(float));
	memset(final_labels, 0, imagewidth*imageheight*sizeof(unsigned char));
	
	loadfiles(fileName); 

	im_orig = new unsigned char[3*imagewidth*imageheight];
	unsigned char *im = rgbImage.GetData();

	for(int i = 0; i < imageheight; i++)
		for(int j = 0; j < imagewidth; j++)
		{
			im_orig[3*((imageheight-(i+1))*imagewidth+j)+0] = im[3*(i*imagewidth+j)+0];
			im_orig[3*((imageheight-(i+1))*imagewidth+j)+1] = im[3*(i*imagewidth+j)+1];
			im_orig[3*((imageheight-(i+1))*imagewidth+j)+2] = im[3*(i*imagewidth+j)+2];
		}

}

void loadfiles(char *fileName)
{
	char *filename = new char[200]; 

# ifndef ALE_UNARY_CAMVID
	sprintf(filename, "Pascal\\Result\\Dense\\%s.c_unary", fileName);
# endif
	// read unary potential
	ProbImage im_unary;

# ifndef ALE_UNARY_CAMVID
	im_unary.decompress(filename);

	float *temp_data = new float[imagewidth*imageheight*num_of_labels];
	for(int i = 0; i < imagewidth*imageheight*num_of_labels; i++)
		temp_data[i] = (im_unary.data())[i]; 

	for(int i = 0; i < imagewidth*imageheight*num_of_labels; i++)
	{
		if(temp_data[i] == 0.0)	temp_data[i] = 0.00001; 
		if(temp_data[i] == 1.0) temp_data[i] = 0.9998; 
		dataCost[i] = -3.0 * log(temp_data[i]);	
	}
	
	delete []temp_data; 

# else  //ALE_UNARY_CAMVID

	// start ALE boosted per pixel unary potential
	
	// Uncomment following lines if you wants to use unary potentials evaluated using ALE

/************************************************************************************/	
	double *classifier_val_reverse = new double[imagewidth*imageheight*num_of_labels];
	double *classifier_val = new double[imagewidth*imageheight*num_of_labels];

	for(int i = 0; i < imagewidth*imageheight*num_of_labels; i++)
	{
		classifier_val_reverse[i] = 0.0; 
		classifier_val[i] = 0.0; 
	}

	//char *filename = new char[200]; 
	sprintf(filename, "CamVid\\Result\\CDenseALE\\%s.dns", fileName);
	FILE *fp = fopen(filename, "rb");
	int temp[3]; 
	fread(temp, sizeof(int), 3, fp);
	fread(classifier_val_reverse, sizeof(double), imagewidth*imageheight*num_of_labels, fp);
	fclose(fp);
	
	// start
	for(int k = 0; k < num_of_labels; k++)
	{
		for(int i = 0; i < imageheight; i++)
		{
			for(int j = 0; j < imagewidth; j++)
			{
				classifier_val[num_of_labels*((imageheight-(1+i))*imagewidth+j)+k] = classifier_val_reverse[num_of_labels*(i*imagewidth+j)+k]; 	
			}
		}
	}

    // calculating the per pixel dataCost
	double sum = 0.0; 
	for(int i = 0; i < imageheight; i++)
		for(int j = 0; j < imagewidth; j++)
		{
			sum = 0.0; double data_cost_val = 0.0; 
			for(int l = 0; l < num_of_labels; l++) 
				sum += exp(classifier_val[num_of_labels*(i*imagewidth+j)+l]);

			for(int l = 0; l < num_of_labels; l++)  
			{
				data_cost_val = -1.0 * log(exp(classifier_val[num_of_labels*(i*imagewidth+j)+l]) /sum);
				dataCost[num_of_labels*(i*imagewidth+j)+l] = data_cost_val;	
			}
		}

	delete []classifier_val; 
	delete []classifier_val_reverse; 
/********************************************************************************/
# endif
	

	// end ALE boosted per pixel unary potential



# ifndef ALE_UNARY_CAMVID //subarna
	// reading cooc unary and pairwise terms
	//char *filename2 = new char[200]; 
	sprintf(filename, "Pascal\\Result\\Cooccurrence\\cooccurence.dat");
	int cooc_total = 0; 
	FILE *fp_in = fopen(filename, "r");
	fscanf(fp_in, "%d", &cooc_total);
	
	for(int i = 0; i < num_of_labels; i++)
	{
		float temp_val = 0; 
		fscanf(fp_in, "%f", &temp_val);
		cooc_unary[i] = temp_val; 
	}

	for(int i = 0; i < num_of_labels; i++)
		for(int j = 0; j < num_of_labels; j++)
		{
			float temp_val = 0; 
			fscanf(fp_in, "%f", &temp_val);
			cooc_pairwise[i*num_of_labels+j] = temp_val; 
		}
# endif

	delete []filename;
}

# ifdef ALE_UNARY_CAMVID
unsigned char my_label(unsigned char label)
{
	unsigned char lab = label;
	switch(lab)
	{
		case 1:
			lab = 1; break;
		case 2:
			lab = 3; break;
		case 3:
			lab = 7; break;
		case 4:
			lab = 12; break;
		case 5:
			lab = 15; break;
		case 6:
			lab = 21; break;
		case 7:
			lab = 24; break;
		case 8:
			lab = 28; break;
		case 9:
			lab = 31; break;
		case 10:
			lab = 36; break;
		case 11:
			lab = 38; break;
		default:
			lab = 0; 
	}
	return lab;
}
# endif


void labeltorgb(short *map, char *filename)
{
	for(int i1 = 0; i1 < imageheight; i1++) 
		for(int j1 = 0; j1 < imagewidth; j1++) 
		{
			final_labels[(imageheight-(1+i1))*imagewidth+j1] = (unsigned char)(map[i1*imagewidth+j1])+1; 
		}

	LRgbImage rgbImage1(imagewidth, imageheight);
	unsigned char *rgbData = rgbImage1.GetData(); 
	for(int i1 = 0; i1 < imagewidth*imageheight; i1++, rgbData += 3)
	{
# ifndef ALE_UNARY_CAMVID
		unsigned char lab = final_labels[i1]; 
# else
		unsigned char lab = my_label(final_labels[i1]); 
# endif
		if(lab == 0) 
			rgbData[0] = rgbData[1] = rgbData[2] = 0; 
		else
		{
# ifndef ALE_UNARY_CAMVID
			lab--;
# endif
			rgbData[0] = rgbData[1] = rgbData[2] = 0;
			for(int i = 0; lab > 0; i++, lab >>= 3)
			{
				rgbData[0] |= (unsigned char) (((lab >> 0) & 1) << (7 - i));
				rgbData[1] |= (unsigned char) (((lab >> 1) & 1) << (7 - i));
				rgbData[2] |= (unsigned char) (((lab >> 2) & 1) << (7 - i));
			}
		}
	}
	char *final_filename = new char[200];

# ifndef ALE_UNARY_CAMVID
	sprintf(final_filename, "Pascal\\Result\\Crf\\%s.png", filename); 
# else
	sprintf(final_filename, "CamVid\\Result\\Crf\\%s.png", filename); 
# endif
	rgbImage1.Save(final_filename);
}


//subarna
void labeltorgbM(int num_images, short *map, char **filename)
{
	int offset;
	for (int n=0; n < num_images; n++ )
	{
		offset = n*imageheight*imagewidth;
		for(int i1 = 0; i1 < imageheight; i1++)
		{
			for(int j1 = 0; j1 < imagewidth; j1++) 
			{
				final_labels[(imageheight-(1+i1))*imagewidth+j1] = (unsigned char)(map[offset + i1*imagewidth+j1])+1; 
			}
		}

		LRgbImage rgbImage1(imagewidth, imageheight);
		unsigned char *rgbData = rgbImage1.GetData(); 
		for(int i1 = 0; i1 < imagewidth*imageheight; i1++, rgbData += 3)
		{
	# ifndef ALE_UNARY_CAMVID
			unsigned char lab = final_labels[i1]; 
	# else
			unsigned char lab = my_label(final_labels[i1]); 
	# endif
			if(lab == 0) 
				rgbData[0] = rgbData[1] = rgbData[2] = 0; 
			else
			{
	# ifndef ALE_UNARY_CAMVID
				lab--;
	# endif
				rgbData[0] = rgbData[1] = rgbData[2] = 0;
				for(int i = 0; lab > 0; i++, lab >>= 3)
				{
					rgbData[0] |= (unsigned char) (((lab >> 0) & 1) << (7 - i));
					rgbData[1] |= (unsigned char) (((lab >> 1) & 1) << (7 - i));
					rgbData[2] |= (unsigned char) (((lab >> 2) & 1) << (7 - i));
				}
			}
		}
		char *final_filename = new char[1024];

# ifndef ALE_UNARY_CAMVID
		sprintf(final_filename, "Pascal\\Result\\Crf\\%s.png", filename[n]); 
# else
		sprintf(final_filename, "CamVid\\Result\\Crf\\%s.png", filename[n]); 
# endif
		rgbImage1.Save(final_filename);

	}
}

# ifndef MATLAB_DEMO
void read_files()
{
	char *trainfile = new char[200];

# ifndef ALE_UNARY_CAMVID
	sprintf(trainfile,  "Pascal\\Result\\Test.txt");
# else
	sprintf(trainfile,  "CamVid\\Result\\CTest.txt");
# endif

	TrainFileNames = new char*[num_files];
	for(int i = 0; i < num_files; i++)
		TrainFileNames[i] = new char[100];

	int part1, part3; char part2; 

# ifdef ALE_UNARY_CAMVID
	char part4, part5, part6;
# endif


	char *fileid = new char[100]; 
	FILE *fp = fopen(trainfile, "r"); //string fileid; 
	for(int i = 0; i < num_files; i++)
	{

# ifndef ALE_UNARY_CAMVID
		//subarna: VOC
		fscanf(fp, "%d %c %d", &part1, &part2, &part3);
	//	printf("%d %d\n", part1, part3); //getch(); 
		sprintf(TrainFileNames[i], "%d_%06d", part1, part3);

# else
		//subarna: CAMVID
		fscanf(fp, "%d %c %c %c %d %c", &part1, &part4, &part5, &part2, &part3, &part6);
	//	printf("%d %d\n", part1, part3); //getch(); 
		sprintf(TrainFileNames[i], "%03dTP_%03d", part1, part3);
# endif

	}
	fclose(fp);

}
# endif

# if defined(ALE_UNARY_CAMVID) & defined(SUPERVOXEL)
void set_ho_layers()
{
	num_of_layers = 4; //13; //subarna
	layers_dir = new char*[num_of_layers]; 
	
	layers_dir[0] = "CamVid\\Result\\Segmentation\\camvid_1TP\\00\\";
	layers_dir[1] = "CamVid\\Result\\Segmentation\\camvid_1TP\\01\\";
	layers_dir[2] = "CamVid\\Result\\Segmentation\\camvid_1TP\\02\\";
	layers_dir[3] = "CamVid\\Result\\Segmentation\\camvid_1TP\\03\\";

	//layers_dir[4] = "CamVid\\Result\\Segmentation\\camvid_1TP\\04\\";
	//layers_dir[5] = "CamVid\\Result\\Segmentation\\camvid_1TP\\05\\";
	//layers_dir[6] = "CamVid\\Result\\Segmentation\\camvid_1TP\\06\\";
	//layers_dir[7] = "CamVid\\Result\\Segmentation\\camvid_1TP\\07\\";
	
	//layers_dir[8] = "CamVid\\Result\\Segmentation\\camvid_1TP\\08\\";
	//layers_dir[9] = "CamVid\\Result\\Segmentation\\camvid_1TP\\09\\";
	//layers_dir[10] = "CamVid\\Result\\Segmentation\\camvid_1TP\\10\\";
	//layers_dir[11] = "CamVid\\Result\\Segmentation\\camvid_1TP\\11\\";
	//layers_dir[12] = "CamVid\\Result\\Segmentation\\camvid_1TP\\12\\";
		
	ho_stats_pot = "CamVid\\Result\\Stats\\"; 
	ho_seg_ext = "msh";
	ho_sts_ext = "sts"; 
}

# elif defined(ALE_UNARY_CAMVID)
void set_ho_layers()
{
	num_of_layers = 3; 
	layers_dir = new char*[num_of_layers]; 
	
	layers_dir[0] = "CamVid\\Result\\Segmentation\\MeanShift30x03\\";
	layers_dir[1] = "CamVid\\Result\\Segmentation\\MeanShift30x06\\";
	layers_dir[2] = "CamVid\\Result\\Segmentation\\MeanShift30x09\\";
		
	ho_stats_pot = "CamVid\\Result\\Stats\\"; 
	ho_seg_ext = "msh";
	ho_sts_ext = "sts"; 
}
# else
void set_ho_layers()
{
	num_of_layers = 10; 
	layers_dir = new char*[num_of_layers]; 
	
	layers_dir[0] = "Pascal\\Result\\Segmentation\\KMeans30\\";
	layers_dir[1] = "Pascal\\Result\\Segmentation\\KMeans40\\";
	layers_dir[2] = "Pascal\\Result\\Segmentation\\KMeans50\\";
	layers_dir[3] = "Pascal\\Result\\Segmentation\\KMeans60\\";
	layers_dir[4] = "Pascal\\Result\\Segmentation\\KMeans80\\";
	layers_dir[5] = "Pascal\\Result\\Segmentation\\KMeans100\\";
	layers_dir[6] = "Pascal\\Result\\Segmentation\\MeanShift70x70\\";
	layers_dir[7] = "Pascal\\Result\\Segmentation\\MeanShift70x100\\";
	layers_dir[8] = "Pascal\\Result\\Segmentation\\MeanShift100x70\\";
	layers_dir[9] = "Pascal\\Result\\Segmentation\\MeanShift100x100\\";
		
	ho_stats_pot = "Pascal\\Result\\Stats\\"; 
	ho_seg_ext = "seg";
	ho_sts_ext = "sts"; 
}
# endif


void set_det_layers()
{
	det_seg_dir = "Pascal\\Result\\Detectors\\seg\\"; 
	det_bb_dir = "Pascal\\Result\\Detectors\\bb\\"; 
	det_seg_ext = "seg"; 
	det_bb_ext = "bb"; 
}