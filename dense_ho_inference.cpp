
/*
This software is an environment for pixel-wise labelling problems, designed mainly for object-
class segmentation problem and described in detail in

Vibhav Vineet, Jonathan Warrell, Philip H.S. Torr
Filter-based Mean-Field Inference for Random Fields with Higher Order Terms and Product Label-Spaces
Proceeding of the twelfth European Conference on Computer Vision, 2012.

This software is free ONLY for research purposes. If you want to use any part of the code you
should cite this paper. 

 THIS SOFTWARE IS PROVIDED BY Vibhav Vineet ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL Vibhav Vineet BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma warning(disable : 4305)

#include "densecrf.h"
#include <cstdio>
#include <cmath>
#include "util.h"
#include <time.h> 
#include <iostream>
#include <conio.h>
#include "main.h"


int main( int argc, char* argv[])
{
	ilInit();

	num_files = 38; 

# ifndef MATLAB_DEMO
	read_files(); 
# else
   num_files = 1; 
   FileName = new char[100];
   strcpy(FileName, argv[1]);
# endif

   printf("code running");

# ifdef MIMAGE
    int num_images = num_files;
    mem_initM(TrainFileNames, num_images);

	short *map = new short[num_images*imagewidth*imageheight];

	DenseCRF2D crf_planeM(num_images, imagewidth, imageheight, num_of_labels);

	// unary
	crf_planeM.setUnaryEnergyM(dataCost);

# ifdef INTER
	int inter = 1;
# else
	int inter = 0;
# endif

	// pairwise
	crf_planeM.addPairwiseGaussianM( num_images, inter, 3, 3, 3 );

	crf_planeM.addPairwiseBilateralM( num_images, inter, 50, 50, 15, 15, 15, im_orig, 5);
	
	int ho_on = 0, ho_det = 0, ho_cooc = 0;
	//// set PN potts ho_order
	ho_on = 1; 
	
	crf_planeM.set_ho(ho_on);
	if(ho_on)
	{
		
		set_ho_layers(); 
		crf_planeM.ho_mem_initM(num_images, imagewidth, imageheight, layers_dir, num_of_layers, ho_stats_pot, ho_seg_ext, ho_sts_ext, 0.0006, 1.0); 
		
		# ifndef MATLAB_DEMO
			crf_planeM.readSegmentsM(num_images, TrainFileNames); 
		# else
			crf_planeM.readSegmentsM(num_images, TrainFileNames); 
		# endif  // MATLAB_DEMO
	}
	// start inference 
	clock_t start=clock(); 
	crf_planeM.mapM(num_images, 5, map); //subarna: mapM ?
	clock_t end=clock(); 
	printf("time taken %f\n", (end-start)/(float)CLOCKS_PER_SEC);
	
	crf_planeM.del_mem_higherorderM(num_images); //subarna

	// save the output
# ifndef MATLAB_DEMO
	labeltorgbM(num_images, map, TrainFileNames); 
# else
	labeltorgbM(num_images, map, FileName); 
# endif
	del_meminit(); 
	delete[] map; 
	delete[] im_orig; 

# else		
	for(int files = 0; files < num_files; files++)
	{
		printf("solving image id %d\n", files);

# ifndef MATLAB_DEMO
		mem_init(TrainFileNames[files]); 
# else
		mem_init(FileName);
# endif
		
		short *map = new short[imagewidth*imageheight];
		DenseCRF2D crf_plane(imagewidth, imageheight, num_of_labels);

		// unary
		crf_plane.setUnaryEnergy(dataCost);
		
		// pairwise
		crf_plane.addPairwiseGaussian( 3, 3, 3 );
		crf_plane.addPairwiseBilateral( 50, 50, 15, 15, 15, im_orig, 5);
		
		int ho_on = 0, ho_det = 0, ho_cooc = 0; 

# define debug

# if 0
# ifdef debug
		//debug //subarna
		unsigned char *res1;
		char name[255];

		// Do map inference
 		crf_plane.startInference();
 		for( int it=0; it<5; it++ ) {
 			crf_plane.stepInference();
 			
			crf_plane.currentMap(map);
			sprintf(name, "debug_%d", it);
			labeltorgb(map, name); 
 		}
# endif
# endif


# ifdef ALE_UNARY_CAMVID
		//// set PN potts ho_order
		ho_on = 1; 
		crf_plane.set_ho(ho_on);
		if(ho_on)
		{
			set_ho_layers(); 
			crf_plane.ho_mem_init(imagewidth, imageheight, layers_dir, num_of_layers, ho_stats_pot, ho_seg_ext, ho_sts_ext, 0.0006, 1.0); 
	# ifndef MATLAB_DEMO
			crf_plane.readSegments(TrainFileNames[files]); 
	# else
			crf_plane.readSegments(FileName); 
	# endif  // MATLAB_DEMO
		}
# else	// ALE_UNARY_CAMVID	
		////set ho_det
		ho_det = 1;
		crf_plane.set_hodet(ho_det);
		if(ho_det)
		{
			set_det_layers(); 
			crf_plane.det_ho_mem_init(imagewidth, imageheight, det_seg_dir, det_bb_dir, det_seg_ext, det_bb_ext, 0.00005, 1.0);
			crf_plane.det_readSegmentIndex(TrainFileNames[files]); 
		}
		
		//// cooccurrence
		ho_cooc = 1;
		crf_plane.set_hocooc(ho_cooc);
		if(ho_cooc)
		{
			crf_plane.setco_occurrence(cooc_unary, cooc_pairwise, 10.0);
		}
# endif  // ALE_UNARY_CAMVID
		
		// start inference 
		clock_t start=clock(); 
		crf_plane.map(5, map);
		clock_t end=clock(); 
		printf("time taken %f\n", (end-start)/(float)CLOCKS_PER_SEC);
		
		crf_plane.del_mem_higherorder();  


		// save the output
# ifndef MATLAB_DEMO
		labeltorgb(map, TrainFileNames[files]); 
# else
		labeltorgb(map, FileName); 
# endif
		del_meminit(); 
		delete[] map; 
		delete[] im_orig; 

	}
# endif   //subarna: # else MIMAGE

//# ifndef MATLAB_DEMO
# ifndef ALE_UNARY_CAMVID
	printf("\n\nfinished with processing \n\nResults are stored at Pascal/Result/Crf/\n");
# else
	printf("\n\nfinished with processing \n\nResults are stored at CamVid/Result/Crf/\n");
# endif

	//getch(); 
//# endif
		
}

