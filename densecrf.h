
#pragma once

#include <vector>
#include <cstdlib>
//#include "potential.h"


class PairwisePotential{
public:
	virtual ~PairwisePotential();
	virtual void apply( float * out_values, const float * in_values, float * tmp, int value_size ) const = 0;
};
class SemiMetricFunction{
public:
	virtual ~SemiMetricFunction();
	// For two probabilities apply the semi metric transform: v_i = sum_j mu_ij u_j
	virtual void apply( float * out_values, const float * in_values, int value_size ) const = 0;
};


class DenseCRF
{
protected:
	friend class BipartiteDenseCRF;
	
	// Number of variables and labels
	int N_, M_;
	//subarna
	int n_;

	float *unary_, *additional_unary_, *current_, *next_, *tmp_;

	// start cooc

	float *unary_cooc_, *pair_cooc_, *current_cooc_, *next_cooc_, *tmp_cooc_; 
	float cooc_factor; 

	// end cooc


	// Store all pairwise potentials
	std::vector<PairwisePotential*> pairwise_;
	
	// Run inference and return the pointer to the result
	float* runInference( int n_iterations, /*float *un_normalized_val,*/ float relax);

	//subarna
	float* runInferenceM( int num_images, int n_iterations, /*float *un_normalized_val,*/ float relax);
	
	// Auxillary functions
	void expAndNormalize( float* out, const float* in, float scale = 1.0, float relax = 1.0 );
	void expAndNormalize_cooc( float *out_cooc_, float *in_cooc_, float scale = 1.0, float relax = 1.0 );
	
	// Don't copy this object, bad stuff will happen
	DenseCRF( DenseCRF & o ){}
public:
	// Create a dense CRF model of size N with M labels
	DenseCRF( int N, int M );
	virtual ~DenseCRF();
	// Add  a pairwise potential defined over some feature space
	// The potential will have the form:    w*exp(-0.5*|f_i - f_j|^2)
	// The kernel shape should be captured by transforming the
	// features before passing them into this function
	void addPairwiseEnergy( const float * features, int D, float w=1.0f, const SemiMetricFunction * function=NULL );
	
	// Add your own favorite pairwise potential (ownwership will be transfered to this class)
	void addPairwiseEnergy( PairwisePotential* potential );
	
	// Set the unary potential for all variables and labels (memory order is [x0l0 x0l1 x0l2 .. x1l0 x1l1 ...])
	// void setUnaryEnergy( const float * unary );
	void setUnaryEnergy( const float * unary/*, float *cooc_unary, float *cooc_pairwise*/);
	
	// Set the unary potential for a specific variable
	void setUnaryEnergy( int n, const float * unary );

	// subarna
	void setUnaryEnergyM( const float * unary );
	void setUnaryEnergyM( int n, const float * unary );
	
	// Run inference and return the probabilities
	void inference( int n_iterations, float* result, float relax=1.0 );
	
	// Run MAP inference and return the map for each pixel
	//void map( int n_iterations, short int* result, float relax=1.0 );
	
	void map( int n_iterations, short int* result, /*float *pix_prob, float *un_normalized_val,*/ float relax=1.0 );

	//subarna
	void mapM( int num_images, int n_iterations, short int* result, /*float *pix_prob, float *un_normalized_val,*/ float relax=1.0 );
		
	// Step by step inference
	void startInference();
	void stepInference( /*float *un_normalized_val,*/ float relax = 1.0 );
	void currentMap( short * result );

	//subarna
	void stepInferenceM( int num_images, /*float *un_normalized_val,*/ float relax = 1.0 );

	// start cooccurrence 
	void set_hocooc(int);
	void setco_occurrence(float *cooc_unary, float *cooc_pairwise, float coocFactor);
	char addCooc; 
	// end cooccurrence

	// start adding higher order reltated stuffs
	void set_ho(int);
	char addHO; 
	int ho_width, ho_height; 
	float *higher_order;
	char *ho_segment_ext; 
	char *ho_stats_ext; 
	int *segmentIndex, *segmentCount, **baseSegmentCounts, ***baseSegmentIndexes, *segmentationIndex, *segmentimagedata, *temp_segmentimagedata; 
	double *stats_potential, *h_norm; 
	float ho_param1, ho_param2; 
	void readSegments(char *filename); 
	void calculate_higherorder_pot(int layer); 
	void readSegmentIndex(char *, int); 
	void readStatsPotential(char *file_name, int num_of_layers); 
	void add_higher_order(); 
	void mem_init_higherorder(); 
	void set_mem_init_higherorder(); 
	void del_mem_higherorder(); 
	void ho_mem_init(int, int, char **, int, char *, char *, char *, float, float); 
	char **layers_dir_name; 
	char *stats_pot; 
	int num_layers; 
	// end higher order related stuffs

	//subarna
	void ho_mem_initM(int, int, int, char **, int, char *, char *, char *, float, float); 
	void readSegmentsM(int num_images, char **filename); 
	void readStatsPotentialM(int num_images, char **file_name, int num_of_layers, char *stats_pot, char *ho_stats_ext );
	void readSegmentIndexM(int, int, int, char *, int); 
	void calculate_higherorder_potM(int layerM); 
	void del_mem_higherorderM(int num_images); 
	
	// det start
	// start adding higher order reltated stuffs
	void set_hodet(int);
	char addDet; 
	char *det_seg_dir, *det_bb_dir; 
	char *det_segment_ext, *det_bbox_ext; 
	int det_ho_width, det_ho_height; 
	float *det_higher_order, det_param1, det_param2; 
	double *det_h_norm, *det_resp;
	int *det_segmentIndex, det_segmentCount, *det_baseSegmentCounts, **det_baseSegmentIndexes, *det_stats_potential; 
	
	void det_calculate_higherorder_pot(); 
	void det_readSegmentIndex(char *); 
	void det_set_mem_init_higherorder(); 
	void det_del_mem_higherorder(); 
	void det_ho_mem_init(int, int, char *, char *, char *, char *, float, float); 

	//void det_mem_init_higherorder(); 
	// det end
	
public: /* Debugging functions */
	// Compute the unary energy of an assignment
	void unaryEnergy( const short * ass, float * result );
	
	// Compute the pairwise energy of an assignment (half of each pairwise potential is added to each of it's endpoints)
	void pairwiseEnergy( const short * ass, float * result, int term=-1 );

	//subarna
	void unaryEnergyM( int n, const short * ass, float * result );
	void pairwiseEnergyM( int n, const short * ass, float * result, int term=-1 );
};

class DenseCRF2D:public DenseCRF
{
protected:
	// Width, height of the 2d grid
	int W_, H_;

	//subarna
	int n_;
public:
	// Create a 2d dense CRF model of size W x H with M labels
	DenseCRF2D( int W, int H, int M );
	virtual ~DenseCRF2D();
	// Add a Gaussian pairwise potential with standard deviation sx and sy
	void addPairwiseGaussian( float sx, float sy, float w, const SemiMetricFunction * function=NULL );
	
	// Add a Bilateral pairwise potential with spacial standard deviations sx, sy and color standard deviations sr,sg,sb
	void addPairwiseBilateral( float sx, float sy, float sr, float sg, float sb, const unsigned char * im, float w, const SemiMetricFunction * function=NULL );
	
	// Set the unary potential for a specific variable
	void setUnaryEnergy( int x, int y, const float * unary );
	using DenseCRF::setUnaryEnergy;

	// ---------
	// functions for MIMAGE case
	//subarna
	// Create a 2d dense CRF model of size W x H with M labels
	DenseCRF2D( int n, int W, int H, int M );
	// Add a Gaussian pairwise potential with standard deviation sx and sy
	void addPairwiseGaussianM( int n, int inter, float sx, float sy, float w, const SemiMetricFunction * function=NULL );
	
	// Add a Bilateral pairwise potential with spacial standard deviations sx, sy and color standard deviations sr,sg,sb
	void addPairwiseBilateralM( int n, int inter, float sx, float sy, float sr, float sg, float sb, const unsigned char * im, float w, const SemiMetricFunction * function=NULL );
	
	// Set the unary potential for a specific variable
	void setUnaryEnergyM(int n, int x, int y, const float * unary );
	using DenseCRF::setUnaryEnergyM;
};


// This function defines a simplified interface to the permutohedral lattice
// We assume a filter standard deviation of 1
class Permutohedral;
class Filter{
protected:
    int n1_, o1_, n2_, o2_;
    Permutohedral * permutohedral_;
    // Don't copy
    Filter( const Filter& filter ){}
public:
    // Use different source and target features
    Filter( const float * source_features, int N_source, const float * target_features, int N_target, int feature_dim );
    // Use the same source and target features
    Filter( const float * features, int N, int feature_dim );
    //
    ~Filter();
    // Filter a bunch of values
    void filter( const float * source, float * target, int value_size );
};
