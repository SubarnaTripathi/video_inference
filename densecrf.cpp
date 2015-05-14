
#pragma warning(disable : 4305)

#include "densecrf.h"
#include "fastmath.h"
#include "permutohedral.h"
#include "util.h"
#include <cmath>
#include <cstring>

PairwisePotential::~PairwisePotential() {
}
SemiMetricFunction::~SemiMetricFunction() {
}
class PottsPotential: public PairwisePotential
{
protected:
	Permutohedral lattice_;
	PottsPotential( const PottsPotential&o ){}
	int N_;
	float w_;
	float *norm_;
public:
	~PottsPotential()
	{
		deallocate( norm_ );
	}
	PottsPotential(const float* features, int D, int N, float w, bool per_pixel_normalization=true) :N_(N), w_(w) 
	{
		lattice_.init( features, D, N );
		norm_ = allocate( N );
		for ( int i=0; i<N; i++ )
			norm_[i] = 1;
		// Compute the normalization factor
		lattice_.compute( norm_, norm_, 1 );
		if ( per_pixel_normalization ) {
			// use a per pixel normalization
			for ( int i=0; i<N; i++ )
				norm_[i] = 1.0 / (norm_[i]+1e-20);
		}
		else {
			float mean_norm = 0;
			for ( int i=0; i<N; i++ )
				mean_norm += norm_[i];
			mean_norm = N / mean_norm;
			// use a per pixel normalization
			for ( int i=0; i<N; i++ )
				norm_[i] = mean_norm;
		}
	}
	void apply(float* out_values, const float* in_values, float* tmp, int value_size) const {
		lattice_.compute( tmp, in_values, value_size );
		
		for ( int i=0,k=0; i<N_; i++ )
		{
			for ( int j=0; j<value_size; j++, k++ )
			{
				out_values[k] += w_*norm_[i]*tmp[k];
			}
		}
	}
};
class SemiMetricPotential: public PottsPotential{
protected:
	const SemiMetricFunction * function_;
public:
	void apply(float* out_values, const float* in_values, float* tmp, int value_size) const {
		lattice_.compute( tmp, in_values, value_size );

		// To the metric transform
		float * tmp2 = new float[value_size];
		for ( int i=0; i<N_; i++ ) {
			float * out = out_values + i*value_size;
			float * t1  = tmp  + i*value_size;
			function_->apply( tmp2, t1, value_size );
			for ( int j=0; j<value_size; j++ )
				out[j] -= w_*norm_[i]*tmp2[j];
		}
		delete[] tmp2;
	}
	SemiMetricPotential(const float* features, int D, int N, float w, const SemiMetricFunction* function, bool per_pixel_normalization=true) :PottsPotential( features, D, N, w, per_pixel_normalization ),function_(function) {
	}
};



/////////////////////////////
/////  Alloc / Dealloc  /////
/////////////////////////////
DenseCRF::DenseCRF(int N, int M) : N_(N), M_(M) {

	// initialize higher order terms
	addHO = 0;
	addDet = 0; 
	addCooc = 0; 

	unary_ = allocate( N_*M_ );
	additional_unary_ = allocate( N_*M_ );
	current_ = allocate( N_*M_ ); 
	next_ = allocate( N_*M_ ); 
	tmp_ = allocate( 2*N_*M_ ); 
	memset( additional_unary_, 0, sizeof(float)*N_*M_ );
}
DenseCRF::~DenseCRF() 
{
	deallocate( unary_ ); 
	deallocate( additional_unary_ ); 
	deallocate( current_ ); 
	deallocate( next_ ); 
	deallocate( tmp_ );
	for( unsigned int i=0; i<pairwise_.size(); i++ )
		delete pairwise_[i];
}

DenseCRF2D::DenseCRF2D(int W, int H, int M) : DenseCRF(W*H,M), W_(W), H_(H)
{
}

//subarna
DenseCRF2D::DenseCRF2D(int n, int W, int H, int M) : DenseCRF(n*W*H,M), n_(n), W_(W), H_(H)
{
}

DenseCRF2D::~DenseCRF2D()
{
}

/////////////////////////////////
/////  Pairwise Potentials  /////
/////////////////////////////////
void DenseCRF::addPairwiseEnergy (const float* features, int D, float w, const SemiMetricFunction * function) 
{
	if (function)
		addPairwiseEnergy( new SemiMetricPotential( features, D, N_, w, function ) );
	else
		addPairwiseEnergy( new PottsPotential( features, D, N_, w ) );
}

void DenseCRF::addPairwiseEnergy ( PairwisePotential* potential )
{
	pairwise_.push_back( potential );
}

void DenseCRF2D::addPairwiseGaussian ( float sx, float sy, float w, const SemiMetricFunction * function ) 
{
	float * feature = new float [N_*2];
	for( int j=0; j<H_; j++ )
		for( int i=0; i<W_; i++ ){
			feature[(j*W_+i)*2+0] = i / sx;
			feature[(j*W_+i)*2+1] = j / sy;
		}
	addPairwiseEnergy( feature, 2, w, function );
	delete [] feature;
}

//subarna
void DenseCRF2D::addPairwiseGaussianM ( int num_images, int inter, float sx, float sy, float w, const SemiMetricFunction * function ) 
{
	int no_connection = 0 ; 
	float *feature = NULL;
	float dev = 5; //5
	
	if (0 == inter)
		feature = new float [N_*2];
	else
		feature = new float [N_*3];

	//if (0 == inter)
	//	no_connection = 1;

	int offset = N_/num_images;
	int idx;

	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		for( int j=0; j<H_; j++ )
		{
			for( int i=0; i<W_; i++ )
			{
				idx = img_idx*offset + j*W_ + i;
			    if (0 == inter)
				{
					feature[idx*2+0] = i / sx;
					feature[idx*2+1] = j / sy;
				}
				else
				{
					feature[idx*3+0] = i / sx;
					feature[idx*3+1] = j / sy;
					feature[idx*3+2] = (img_idx)*dev;

					if (1 == no_connection)
						feature[idx*3+2] = (img_idx)*99999;
				}
			}
		}
	}

	if (0 == inter)
		addPairwiseEnergy( feature, 2, w, function );
	else
		addPairwiseEnergy( feature, 3, w, function );

	delete [] feature;
}

void DenseCRF2D::addPairwiseBilateral ( float sx, float sy, float sr, float sg, float sb, const unsigned char* im, float w, const SemiMetricFunction * function ) {
	float * feature = new float [N_*5];
	for( int j=0; j<H_; j++ )
		for( int i=0; i<W_; i++ ){
			feature[(j*W_+i)*5+0] = i / sx;
			feature[(j*W_+i)*5+1] = j / sy;
			feature[(j*W_+i)*5+2] = im[(i+j*W_)*3+0] / sr;
			feature[(j*W_+i)*5+3] = im[(i+j*W_)*3+1] / sg;
			feature[(j*W_+i)*5+4] = im[(i+j*W_)*3+2] / sb;
		}
	addPairwiseEnergy( feature, 5, w, function );
	delete [] feature;
}

//subarna
void DenseCRF2D::addPairwiseBilateralM ( int num_images, int inter, float sx, float sy, float sr, float sg, float sb, const unsigned char* im, float w, const SemiMetricFunction * function ) {

	int no_connection = 0 ; //1; 
	float dev = 10; //10
	float *feature = NULL;

	if (0 == inter)
		feature = new float [N_*5];
	else
		feature = new float [N_*6];

	//if (0 == inter)
	//	no_connection = 1;
	

	int offset = N_/num_images;
	int idx;

	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		for( int j=0; j<H_; j++ )
		{
			for( int i=0; i<W_; i++ )
			{	
				idx = img_idx*offset + j*W_ + i;
				if (0 == inter)
				{
					feature[idx*5+0] = i / sx;
					feature[idx*5+1] = j / sy;
					feature[idx*5+2] = im[3*offset + (i+j*W_)*3+0] / sr;
					feature[idx*5+3] = im[3*offset + (i+j*W_)*3+1] / sg;
					feature[idx*5+4] = im[3*offset + (i+j*W_)*3+2] / sb;
				}
				else
				{
					feature[idx*6+0] = i / sx;
					feature[idx*6+1] = j / sy;
					feature[idx*6+2] = (img_idx+1)*dev;

					if (1 == no_connection)
						feature[idx*6+2] = (img_idx+1)*99999;

					feature[idx*6+3] = im[3*offset + (i+j*W_)*3+1] / sr;
					feature[idx*6+4] = im[3*offset + (i+j*W_)*3+2] / sg;
					feature[idx*6+5] = im[3*offset + (i+j*W_)*3+2] / sb;
				}
			}
		}
	}

    if (0 == inter)
		addPairwiseEnergy( feature, 5, w, function );
	else
		addPairwiseEnergy( feature, 6, w, function );

	delete [] feature;
}
//////////////////////////////
/////  Unary Potentials  /////
//////////////////////////////
void DenseCRF::setUnaryEnergy ( const float* unary/*,  float *cooc_unary, float *cooc_pairwise*/) 
{
	memcpy( unary_, unary, N_*M_*sizeof(float) );
}

void DenseCRF::setUnaryEnergy ( int n, const float* unary ) 
{
	memcpy( unary_+n*M_, unary, M_*sizeof(float) );
}
void DenseCRF2D::setUnaryEnergy ( int x, int y, const float* unary ) {
	memcpy( unary_+(x+y*W_)*M_, unary, M_*sizeof(float) );
}

//subarna
void DenseCRF::setUnaryEnergyM ( const float* unary/*,  float *cooc_unary, float *cooc_pairwise*/) 
{
	memcpy( unary_, unary, N_*M_*sizeof(float) );
}
void DenseCRF2D::setUnaryEnergyM ( int n, int x, int y, const float* unary ) {
	memcpy( unary_+(x+y*W_)*M_*n_, unary, n_*M_*sizeof(float) );
}
///////////////////////
/////  Inference  /////
///////////////////////

void DenseCRF::map ( int n_iterations, short* result, float relax) 
{
	// Run inference
	float * prob = runInference( n_iterations, relax );
	
	// Find the map
	for( int i=0; i<N_; i++ ){
		const float * p = prob + i*M_;
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++)
			if( mx < p[j] )
			{
				mx = p[j];
				imx = j;
			}
		result[i] = imx;
	}
}

//subarna
void DenseCRF::mapM (int num_images, int n_iterations, short* result, float relax) 
{
	// Run inference
	float * prob = runInferenceM( num_images, n_iterations, relax );
	
	// Find the map
	for( int i=0; i<N_; i++ ){
		const float * p = prob + i*M_;
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++)
			if( mx < p[j] )
			{
				mx = p[j];
				imx = j;
			}
		result[i] = imx;
	}
}

float* DenseCRF::runInference( int n_iterations, float relax ) 
{
	startInference();

	for( int it=0; it<n_iterations; it++)
		stepInference(relax);

	return current_;
}

//subarna
float* DenseCRF::runInferenceM( int num_images, int n_iterations, float relax ) 
{
	startInference();

	for( int it=0; it<n_iterations; it++)
		stepInferenceM(num_images, relax);

	return current_;
}
void DenseCRF::expAndNormalize ( float* out, const float* in, float scale, float relax ) 
{
	float *V = new float[ N_+10 ];
	for( int i=0; i<N_; i++ ){
		const float * b = in + i*M_;
		// Find the max and subtract it so that the exp doesn't explode
		float mx = scale*b[0];
		for( int j=1; j<M_; j++ )
			if( mx < scale*b[j] )
				mx = scale*b[j];
		float tt = 0;
		for( int j=0; j<M_; j++ )
		{
			V[j] = fast_exp( scale*b[j]-mx );
			tt += V[j];
		}
		// Make it a probability
		for( int j=0; j<M_; j++ )
			V[j] /= tt;
		
		float * a = out + i*M_;
		for( int j=0; j<M_; j++ )
			if (relax == 1)
				a[j] = V[j];
			else
				a[j] = (1-relax)*a[j] + relax*V[j];
	}
	delete[] V;
}



void DenseCRF::expAndNormalize_cooc ( float *cooc_out, float *cooc_in, float scale, float relax ) 
{
	float *V_cooc = new float[ M_+10 ];
	for( int i=0; i<M_; i++ )
	{
		const float * b_cooc = cooc_in + i*2;
		// Find the max and subtract it so that the exp doesn't explode
		float mx_cooc = scale*b_cooc[0];
		for( int j=1; j < 2; j++ )
			if( mx_cooc < scale * b_cooc[j] )
				mx_cooc = scale*b_cooc[j];

		float tt = 0;
		for( int j=0; j<2; j++ )
		{
			V_cooc[j] = fast_exp( scale*b_cooc[j]-mx_cooc );
			tt += V_cooc[j];
		}
		// Make it a probability
		for( int j=0; j<2; j++ )
			V_cooc[j] /= tt;
		
		float * a_cooc = cooc_out + i*2;
		for( int j=0; j<2; j++ )
			a_cooc[j] = V_cooc[j];		
	}

 	delete[] V_cooc;
}



///////////////////
/////  Debug  /////
///////////////////

void DenseCRF::unaryEnergy(const short* ass, float* result) {
	for( int i=0; i<N_; i++ )
		if ( 0 <= ass[i] && ass[i] < M_ )
			result[i] = unary_[ M_*i + ass[i] ];
		else
			result[i] = 0;
}
void DenseCRF::pairwiseEnergy(const short* ass, float* result, int term) 
{
	float * current = allocate( N_*M_ );
	// Build the current belief [binary assignment]
	for( int i=0,k=0; i<N_; i++ )
		for( int j=0; j<M_; j++, k++ )
			current[k] = (ass[i] == j);
	
	for( int i=0; i<N_*M_; i++ )
		next_[i] = 0;
	if (term == -1)
		for( unsigned int i=0; i<pairwise_.size(); i++ )
			pairwise_[i]->apply( next_, current, tmp_, M_ );
	else
		pairwise_[ term ]->apply( next_, current, tmp_, M_ );
	for( int i=0; i<N_; i++ )
		if ( 0 <= ass[i] && ass[i] < M_ )
			result[i] =-next_[ i*M_ + ass[i] ];
		else
			result[i] = 0;
	deallocate( current );
}


void DenseCRF::startInference()
{
	if(addCooc)
	{
		int *total_num_labels = new int[M_];
		memset(total_num_labels, 0, M_);
		for(int i = 0; i < N_; i++)
		{
			int class_label = 0; float temp_unary_cost = unary_[i*M_]; 
			for(int j = 1; j < M_; j++)
			{
				if(temp_unary_cost < unary_[i*M_+j])
				{
					temp_unary_cost = unary_[i*M_+j]; 
					class_label = j; 
				}
			}
			total_num_labels[class_label]++; 
		}

		float pairwise_cooc = 0.0; // float p1, p2, p, p12; 
		for(int i = 0; i < M_; i++)
		{
			if(total_num_labels[i] > 0)
			{
				next_cooc_[2*i+1] = total_num_labels[i];
				next_cooc_[2*i] = 1;
			}
			else
			{
				next_cooc_[2*i+1] = 1;
				next_cooc_[2*i] = 100;
			}
		}

		delete []total_num_labels; 
	}
	// Initialize using the unary energies
	expAndNormalize( current_, unary_, -1 );
	
	if(addCooc)
	{
		expAndNormalize_cooc ( current_cooc_, next_cooc_) ; 
	}
}

void DenseCRF::stepInference( float relax )
{
#ifdef SSE_DENSE_CRF
	__m128 * sse_next_ = (__m128*)next_;
	__m128 * sse_unary_ = (__m128*)unary_;
	__m128 * sse_additional_unary_ = (__m128*)additional_unary_;
#endif
	// Set the unary potential
#ifdef SSE_DENSE_CRF

	static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	
	__m128 val1, val2; 

	for( int i=0; i<(N_*M_-1)/4+1; i++ )
	{
		val1 = _mm_xor_ps(sse_unary_[i], SIGNMASK); 
		val2 = _mm_xor_ps(sse_additional_unary_[i], SIGNMASK); 
		sse_next_[i] = _mm_add_ps(val1, val2); 
	}
#else
	for( int i=0; i<N_*M_; i++ )
		next_[i] = -unary_[i] - additional_unary_[i];
#endif

	//start PN Potts
	if(addHO)
	{
		set_mem_init_higherorder(); 
		int num_of_layers = num_layers; 
		for(int i = 0; i < num_of_layers; i++)
			calculate_higherorder_pot(i); 

		//HO->calculate_higherorder_pot(i, current, next); 

		for(int i = 0; i < N_*M_; i++)
		{
			next_[i] = next_[i] - higher_order[i]; 
		}
	}

	//end PN Potts


	// start add co-occurrence terms
	if(addCooc)
	{
		int *higher_labels = new int[M_];

		for(int i = 0; i < M_; i++)
		{
			if(current_cooc_[2*i] < current_cooc_[2*i+1])
				higher_labels[i] = 1; 
			else 
				higher_labels[i] = 0; 
		}

		float *temp_prob_mult = new float[2*M_];
		float mult_prob = 0.0, mult_prob1 = 0.0; 

		for(int i = 0; i < M_; i++)
		{
			mult_prob = 0.0; mult_prob1 = 0.0; 
			for(int j = 0; j < N_; j++)
			{
				mult_prob = mult_prob + (1.0-current_[j*M_+i]);
				mult_prob1 = mult_prob1 + current_[j*M_+i]; 
			}

			temp_prob_mult[2*i] = mult_prob; 
			temp_prob_mult[2*i+1] = mult_prob1;
		
			if(temp_prob_mult[2*i] < 1e-4) temp_prob_mult[2*i] = 1e-4; 
			if(temp_prob_mult[2*i+1] < 1e-4) temp_prob_mult[2*i+1] = 1e-4; 
		}
	
		float pairwise_cooc = 0.0; float p1, p2, p, p12; 
		for(int i = 0; i < M_; i++)
		{
			pairwise_cooc = 0.0; 
			p1 = unary_cooc_[i]; 
			for(int j = 0; j < M_; j++)
			{
				p2 = unary_cooc_[j]; 
				p12 = pair_cooc_[i*M_+j]; 
				p = 1 - (1 - p12 / p2) * (1 - p12 / p1);
				if(p > 1) p = 1;
				if(p < 1e-6) p = 1e-6;
	
				if(i != j)
					pairwise_cooc = pairwise_cooc - (0.005*N_) * log(p) * current_cooc_[j*2+1];
			}
			next_cooc_[2*i+1] = -pairwise_cooc; 		
		}

		for(int i = 0; i < M_; i++)
		{
			next_cooc_[2*i] = -1.0*cooc_factor*(temp_prob_mult[2*i+1]); 
		}


		float temp_cooc_factor = 1.0; 

		for(int i = 0; i < N_; i++)
		{
			for(int j = 0; j < M_; j++)
			{
				mult_prob = 1.0; 
				next_[i*M_+j] = next_[i*M_+j] - cooc_factor*current_cooc_[2*j];
			}
		}

		delete []temp_prob_mult; 
		delete []higher_labels; 
	}
	// end co-occurrence terms

	// start det
	if(addDet)
	{
		//det_set_mem_init_higherorder(); 
		det_calculate_higherorder_pot(); 

		for(int i = 0; i < N_*M_; i++)
		{
			next_[i] = next_[i] - det_higher_order[i]; 
		}
	}
	// end det


	// pairwise potentials
	for( unsigned int i=0; i<pairwise_.size(); i++)
		pairwise_[i]->apply( next_, current_, tmp_, M_);
	
	// end pairwise

	// Exponentiate and normalize
	expAndNormalize( current_, next_, 1.0, relax );
	if(addCooc)
	{
		expAndNormalize_cooc ( current_cooc_, next_cooc_) ; 
	}
}

//subarna
void DenseCRF::stepInferenceM( int num_images, float relax )
{
#ifdef SSE_DENSE_CRF
	__m128 * sse_next_ = (__m128*)next_;
	__m128 * sse_unary_ = (__m128*)unary_;
	__m128 * sse_additional_unary_ = (__m128*)additional_unary_;
#endif
	// Set the unary potential
#ifdef SSE_DENSE_CRF

	static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	
	__m128 val1, val2; 

	for( int i=0; i<(N_*M_-1)/4+1; i++ )
	{
		val1 = _mm_xor_ps(sse_unary_[i], SIGNMASK); 
		val2 = _mm_xor_ps(sse_additional_unary_[i], SIGNMASK); 
		sse_next_[i] = _mm_add_ps(val1, val2); 
	}
#else
	for( int i=0; i<N_*M_; i++ )
		next_[i] = -unary_[i] - additional_unary_[i];
#endif

	//start PN Potts
	if(addHO)
	{
		set_mem_init_higherorder(); 
		int num_of_layers = num_layers; 

		for (int n=0; n < num_images; n++)
		{
			for(int i = 0; i < num_of_layers; i++)
			{
				//calculate_higherorder_pot(i); 
				calculate_higherorder_potM(n*num_of_layers + i);
			}
		}

		//HO->calculate_higherorder_pot(i, current, next); 

		for(int i = 0; i < N_*M_; i++)
		{
			next_[i] = next_[i] - higher_order[i]; 
		}
	}

	//end PN Potts


	// start add co-occurrence terms
	if(addCooc)
	{
		int *higher_labels = new int[M_];

		for(int i = 0; i < M_; i++)
		{
			if(current_cooc_[2*i] < current_cooc_[2*i+1])
				higher_labels[i] = 1; 
			else 
				higher_labels[i] = 0; 
		}

		float *temp_prob_mult = new float[2*M_];
		float mult_prob = 0.0, mult_prob1 = 0.0; 

		for(int i = 0; i < M_; i++)
		{
			mult_prob = 0.0; mult_prob1 = 0.0; 
			for(int j = 0; j < N_; j++)
			{
				mult_prob = mult_prob + (1.0-current_[j*M_+i]);
				mult_prob1 = mult_prob1 + current_[j*M_+i]; 
			}

			temp_prob_mult[2*i] = mult_prob; 
			temp_prob_mult[2*i+1] = mult_prob1;
		
			if(temp_prob_mult[2*i] < 1e-4) temp_prob_mult[2*i] = 1e-4; 
			if(temp_prob_mult[2*i+1] < 1e-4) temp_prob_mult[2*i+1] = 1e-4; 
		}
	
		float pairwise_cooc = 0.0; float p1, p2, p, p12; 
		for(int i = 0; i < M_; i++)
		{
			pairwise_cooc = 0.0; 
			p1 = unary_cooc_[i]; 
			for(int j = 0; j < M_; j++)
			{
				p2 = unary_cooc_[j]; 
				p12 = pair_cooc_[i*M_+j]; 
				p = 1 - (1 - p12 / p2) * (1 - p12 / p1);
				if(p > 1) p = 1;
				if(p < 1e-6) p = 1e-6;
	
				if(i != j)
					pairwise_cooc = pairwise_cooc - (0.005*N_) * log(p) * current_cooc_[j*2+1];
			}
			next_cooc_[2*i+1] = -pairwise_cooc; 		
		}

		for(int i = 0; i < M_; i++)
		{
			next_cooc_[2*i] = -1.0*cooc_factor*(temp_prob_mult[2*i+1]); 
		}


		float temp_cooc_factor = 1.0; 

		for(int i = 0; i < N_; i++)
		{
			for(int j = 0; j < M_; j++)
			{
				mult_prob = 1.0; 
				next_[i*M_+j] = next_[i*M_+j] - cooc_factor*current_cooc_[2*j];
			}
		}

		delete []temp_prob_mult; 
		delete []higher_labels; 
	}
	// end co-occurrence terms

	// start det
	if(addDet)
	{
		//det_set_mem_init_higherorder(); 
		det_calculate_higherorder_pot(); 

		for(int i = 0; i < N_*M_; i++)
		{
			next_[i] = next_[i] - det_higher_order[i]; 
		}
	}
	// end det


	// pairwise potentials
	for( unsigned int i=0; i<pairwise_.size(); i++)
		pairwise_[i]->apply( next_, current_, tmp_, M_);
	
	// end pairwise

	// Exponentiate and normalize
	expAndNormalize( current_, next_, 1.0, relax );
	if(addCooc)
	{
		expAndNormalize_cooc ( current_cooc_, next_cooc_) ; 
	}
}

void DenseCRF::currentMap( short * result )
{
	// Find the map
	for( int i=0; i<N_; i++ ){
		const float * p = current_ + i*M_;
		// Find the max and subtract it so that the exp doesn't explode
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++ )
			if( mx < p[j] ){
				mx = p[j];
				imx = j;
			}
		result[i] = imx;
	}
}


// start co-occurrence

void DenseCRF::set_hocooc(int add_Cooc)
{
	addCooc = 0;
	if(add_Cooc) 
	{
		addCooc = 1;	
		
		unary_cooc_ = new float[M_];
		memset(unary_cooc_, 0, sizeof(float)*M_);
		pair_cooc_ = new float[M_*M_];
		memset(pair_cooc_, 0, sizeof(float)*M_*M_);
		current_cooc_ = allocate(2*M_); 
		next_cooc_ = allocate(2*M_);	
	}
}

void DenseCRF::setco_occurrence(float *cooc_unary, float *cooc_pairwise,  float coocFactor) 
{
	if(addCooc)
	{
		memcpy( unary_cooc_, cooc_unary, M_*sizeof(float) );
		for(int i = 0; i < M_; i++)
			unary_cooc_[i] = unary_cooc_[i]/600.0;

		memcpy(pair_cooc_, cooc_pairwise, M_*M_*sizeof(float));

		for(int i = 0; i < M_*M_; i++)
			pair_cooc_[i] = pair_cooc_[i]/600.0;

		cooc_factor = coocFactor; 
	}
}

// end co-occurrence



// start adding higher order potential related stuffs

void DenseCRF::set_ho(int add_ho)
{
	addHO = 0;
	if(add_ho) 
		addHO = 1; 		
}

void DenseCRF::ho_mem_init(int imagewidth, int imageheight, char **layers_name, int num_of_layers, char *ho_stats, char *ho_seg_ext, char *ho_sts_ext, float hoParam1, float hoParam2)
{
	if(addHO)
	{
		num_layers = num_of_layers; 
		layers_dir_name = new char*[num_layers]; 
		for(int i = 0; i < num_layers; i++)
			layers_dir_name[i] = layers_name[i]; 

		stats_pot = ho_stats; 

		ho_segment_ext = ho_seg_ext; 
		ho_stats_ext = ho_sts_ext; 

		ho_width = imagewidth; ho_height = imageheight; 
		higher_order = new float[ho_width*ho_height*M_];

		for(int i = 0; i < ho_width*ho_height*M_; i++)
		{
			higher_order[i] = 0.0; 
		}

		temp_segmentimagedata = new int[N_];
		segmentimagedata = new int[N_];
		baseSegmentIndexes = new int**[num_layers]; 
		segmentCount = new int[num_layers];
		baseSegmentCounts = new int*[num_layers];
		ho_param1 = hoParam1; 
		ho_param2 = hoParam2; 
	}
}

//subarna
void DenseCRF::ho_mem_initM(int num_images, int imagewidth, int imageheight, char **layers_name, int num_of_layers, char *ho_stats, char *ho_seg_ext, char *ho_sts_ext, float hoParam1, float hoParam2)
{
	if(addHO)
	{
		num_layers = num_of_layers; 
		layers_dir_name = new char*[num_layers]; 
		for(int i = 0; i < num_layers; i++)
			layers_dir_name[i] = layers_name[i]; 

		stats_pot = ho_stats; 

		ho_segment_ext = ho_seg_ext; 
		ho_stats_ext = ho_sts_ext; 

		ho_width = imagewidth; ho_height = imageheight; 
		higher_order = new float[num_images*ho_width*ho_height*M_];

		for(int i = 0; i < num_images*ho_width*ho_height*M_; i++)
		{
			higher_order[i] = 0.0; 
		}

		temp_segmentimagedata = new int[N_];
		segmentimagedata = new int[N_];
		baseSegmentIndexes = new int**[num_images*num_layers]; 
		segmentCount = new int[num_images*num_layers];
		baseSegmentCounts = new int*[num_images*num_layers];
		ho_param1 = hoParam1; 
		ho_param2 = hoParam2; 
	}
}

void DenseCRF::readSegments(char *filename)
{
	if (addHO)
	{
		char *seg_file_name = new char[200];

		// read segments
		for(int i = 0; i < num_layers; i++)
		{
			sprintf(seg_file_name, "%s%s.%s", layers_dir_name[i], filename, ho_segment_ext);
			readSegmentIndex(seg_file_name, i);
		}

		// read stats potential
		sprintf(seg_file_name, "%s%s.%s", stats_pot, filename, ho_stats_ext);
		readStatsPotential(seg_file_name, num_layers); 
			
		delete []seg_file_name; 
	}
}

//subarna
void DenseCRF::readSegmentsM(int num_images, char **filename)
{
	if (addHO)
	{
		char *seg_file_name = new char[1024];

		for (int n=0; n <num_images; n++ )
		{
			// read segments
			for(int i = 0; i < num_layers; i++)
			{
				sprintf(seg_file_name, "%s%s.%s", layers_dir_name[i], filename[n], ho_segment_ext);
				readSegmentIndexM(n, n*num_layers, num_images, seg_file_name, i);
			}

			// read stats potential
			////////////////////////////////////////////////////////////////////////////
			//sprintf(seg_file_name, "%s%s.%s", stats_pot, filename[n], ho_stats_ext);

			////////////////////////////////////////////////////////////////////////////
			// subarna: start debugging from here
			//readStatsPotentialM(num_images, seg_file_name, num_layers); 
			
		}

		readStatsPotentialM(num_images, filename, num_layers, stats_pot, ho_stats_ext ); 
			
		delete []seg_file_name; 
	}
}

void DenseCRF::readSegmentIndex(char *file_name, int layer)
{
	if(addHO)
	{
		FILE *fin = fopen(file_name, "rb");
		int width, height, bands; 
		fread(&width, 4, 1, fin); fread(&height, 4, 1, fin); fread(&bands, 4, 1, fin);
		fread(temp_segmentimagedata, sizeof(int), width * height * bands, fin);
		fclose(fin);

		for(int i = 0; i < height; i++)
			for(int j = 0; j < width; j++)
				segmentimagedata[i*width+j] = temp_segmentimagedata[(height-(i+1))*width+j];
		
		segmentCount[layer] = 0;
		int points = width*height; 
		for(int j = 0; j < points; j++) 
			if(segmentimagedata[j]+1 > segmentCount[layer]) 
				segmentCount[layer] = segmentimagedata[j] + 1;
	
		baseSegmentCounts[layer] = new int[segmentCount[layer]];

		memset(baseSegmentCounts[layer], 0, segmentCount[layer] * sizeof(int));

		for(int j = 0; j < points; j++) 
			baseSegmentCounts[layer][segmentimagedata[j]]++;

		int temp_count = 0; 

		for(int i = 0; i < segmentCount[layer]; i++)
			temp_count += baseSegmentCounts[layer][i]; 
	
		baseSegmentIndexes[layer] = new int *[segmentCount[layer]];
		for(int j = 0; j < segmentCount[layer]; j++) 
			baseSegmentIndexes[layer][j] = new int[baseSegmentCounts[layer][j]];

		segmentationIndex = new int[segmentCount[layer]];
		memset(segmentationIndex, 0, segmentCount[layer] * sizeof(int));

		for(int j = 0; j < points; j++)
		{
			baseSegmentIndexes[layer][segmentimagedata[j]][segmentationIndex[segmentimagedata[j]]] = j;
			segmentationIndex[segmentimagedata[j]]++;
		}
		delete []segmentationIndex;
	}
}

//subarna
void DenseCRF::readSegmentIndexM(int curr, int cum_layer, int num_images, char *file_name, int layer)
{
	//subarna
	long offset;
	//long idx;

	if(addHO)
	{
		FILE *fin = fopen(file_name, "rb");
		int width, height, bands; 
		fread(&width, 4, 1, fin); fread(&height, 4, 1, fin); fread(&bands, 4, 1, fin);
		fread(temp_segmentimagedata, sizeof(int), width * height * bands, fin);
		fclose(fin);

		//subarna
		offset = curr*height*width;

		for(int i = 0; i < height; i++)
		{
			for(int j = 0; j < width; j++)
				segmentimagedata[offset + i*width+j] = temp_segmentimagedata[(height-(i+1))*width+j];
		}
		
		//subarna
		segmentCount[cum_layer+layer] = 0;
		int points = width*height; 
				
		long offset = curr*points;
		for(int j = 0; j < points; j++) 
		{
			if(segmentimagedata[offset+j]+1 > segmentCount[cum_layer+layer]) 
				segmentCount[cum_layer+layer] = segmentimagedata[offset+j] + 1;
		}
	
		baseSegmentCounts[cum_layer+layer] = new int[segmentCount[cum_layer+layer]];

		memset(baseSegmentCounts[cum_layer+layer], 0, segmentCount[cum_layer+layer] * sizeof(int));

		
		for(int j = 0; j < points; j++) 
			baseSegmentCounts[cum_layer+layer][segmentimagedata[offset+j]]++;

		int temp_count = 0; 

		for(int i = 0; i < segmentCount[cum_layer+layer]; i++)
			temp_count += baseSegmentCounts[cum_layer+layer][i]; 
	
		baseSegmentIndexes[cum_layer+layer] = new int *[segmentCount[cum_layer+layer]];
		for(int j = 0; j < segmentCount[cum_layer+layer]; j++) 
			baseSegmentIndexes[cum_layer+layer][j] = new int[baseSegmentCounts[cum_layer+layer][j]];

		segmentationIndex = new int[segmentCount[cum_layer+layer]];
		memset(segmentationIndex, 0, segmentCount[cum_layer+layer] * sizeof(int));

		for(int j = 0; j < points; j++)
		{
			baseSegmentIndexes[cum_layer+layer][segmentimagedata[offset+j]][segmentationIndex[segmentimagedata[offset+j]]] = offset + j;
			segmentationIndex[segmentimagedata[offset+j]]++;
		}
		delete []segmentationIndex;
	}
}

void DenseCRF::readStatsPotential(char *file_name, int num_of_layers)
{
	if(addHO)
	{
		int totalsegmentCount = 0; 
	
		for(int i = 0; i < num_of_layers; i++)
		{
			totalsegmentCount += segmentCount[i]; 
		}

		stats_potential = new double[totalsegmentCount*M_];
		for(int i = 0; i < totalsegmentCount*M_; i++)
			stats_potential[i] = 0.0; 

		FILE *fin1 = fopen(file_name, "rb");
		fread(stats_potential, sizeof(double), totalsegmentCount*M_, fin1);
		fclose(fin1); 

		float sum = 0; int start_loc = 0; 
		for(int i1 = 0; i1 < num_of_layers; i1++)
		{
			for(int i = 0; i < segmentCount[i1]; i++)
			{
				sum = 0.0; 
				for(int k = 0; k < M_; k++) 
				{
					sum += exp(stats_potential[start_loc+k]);
				}
				for(int k = 0; k < M_; k++) 
					stats_potential[start_loc+k] = exp(stats_potential[start_loc+k]) / sum; 
				start_loc = start_loc + M_; 
			}
		}
	}
}

//subarna
void DenseCRF::readStatsPotentialM(int num_images, char **file_name, int num_of_layers, char *stats_pot, char *ho_stats_ext )
{
	if(addHO)
	{
		char *seg_file_name = new char[1024];

		int totalsegmentCount = 0; 

		int *totsegCount = (int*)malloc(num_images*sizeof(int));
	
		/////////////
		for (int n = 0; n < num_images; n++ )
		{
			totsegCount[n] = 0;
			for(int i = 0; i < num_of_layers; i++)
			{
				totsegCount[n] += segmentCount[n*num_of_layers+i]; 
				
			}
			totalsegmentCount += totsegCount[n];
		}

		stats_potential = new double[totalsegmentCount*M_];
		for(int i = 0; i < totalsegmentCount*M_; i++)
			stats_potential[i] = 0.0; 

		int so_far = 0;
		for (int n=0; n < num_images; n++)
		{
			sprintf(seg_file_name, "%s%s.%s", stats_pot, file_name[n], ho_stats_ext);

			FILE *fin1 = fopen(seg_file_name, "rb");
			fread(stats_potential + so_far, sizeof(double), totsegCount[n]*M_, fin1);
			fclose(fin1); 

			so_far += totsegCount[n]*M_; //subarna
		}


		//subarna
		float sum = 0; int start_loc = 0; 
		for (int n=0; n < num_images; n++ )
		{
			for(int i1 = 0; i1 < num_of_layers; i1++)
			{
				for(int i = 0; i < segmentCount[n*num_of_layers + i1]; i++)
				{
					sum = 0.0; 
					for(int k = 0; k < M_; k++) 
					{
						sum += exp(stats_potential[start_loc+k]);
					}
					for(int k = 0; k < M_; k++) 
						stats_potential[start_loc+k] = exp(stats_potential[start_loc+k]) / sum; 
					
					start_loc = start_loc + M_; 
				}
			}
		}

		delete []seg_file_name; 
		delete []totsegCount;
	}
}


void DenseCRF::calculate_higherorder_pot(int layer)
{
	if(addHO)
	{
		h_norm = new double[segmentCount[layer]*M_];

		float norm_val = 0.0; 
		int basesegmentcounts = 0; 
		int curr_pix_label = 0, curr_pix_index; // int x, y; 
	//	int neigh_pix_index, neigh_pix_label; 
	
		double higher_order_prob; 
	
		for(int i = 0; i < segmentCount[layer]; i++)
			for(int j = 0; j < M_; j++)
				h_norm[i*M_+j] = 1.0; 
	
		for(int i = 0; i < segmentCount[layer]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layer][i];
			higher_order_prob = 1.0;
			for(int j = 0; j < M_; j++)
			{
				higher_order_prob = 1.0; 
				for(int k = 0; k < basesegmentcounts; k++)
				{
					curr_pix_index = baseSegmentIndexes[layer][i][k];
					higher_order_prob = higher_order_prob * current_[curr_pix_index*M_+j]; 
				}
				h_norm[i*M_+j] = higher_order_prob; 
			}
		}


		double alpha = 0.5, maxcost, weight, costdata = 0.0; int start_loc = 0; 

		for(int i = 0; i < layer; i++)
			start_loc = start_loc + segmentCount[i]; 

		start_loc = start_loc * M_; 
	
		for(int i = 0; i < segmentCount[layer]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layer][i];
		
			weight = 0.3 * basesegmentcounts; 
			maxcost = -weight * log(alpha); 
		
			for(int j = 0; j < basesegmentcounts; j++)
			{
				curr_pix_index = baseSegmentIndexes[layer][i][j]; 
				for(int k = 0; k < M_; k++)
				{
					higher_order_prob = h_norm[i*M_+k]/(current_[curr_pix_index*M_+k]+0.0001);
					costdata = - weight * log(stats_potential[start_loc+k]); 
					higher_order[curr_pix_index*M_+k] += (ho_param1*costdata - ho_param2*higher_order_prob);
				}		
			}	
			start_loc = start_loc+M_; 
		}
		delete []h_norm; 	
	}
}

//subarna
void DenseCRF::calculate_higherorder_potM(int layerM)
{
	if(addHO)
	{
		int start_loc = 0; 
		
		h_norm = new double[segmentCount[layerM]*M_];  //subarna

		float norm_val = 0.0; 
		int basesegmentcounts = 0; 
		int curr_pix_label = 0, curr_pix_index; // int x, y; 
	//	int neigh_pix_index, neigh_pix_label; 
	
		double higher_order_prob; 
	
		for(int i = 0; i < segmentCount[layerM]; i++)
		{
			for(int j = 0; j < M_; j++)
				h_norm[i*M_+j] = 1.0; 
		}
	
		for(int i = 0; i < segmentCount[layerM]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layerM][i];
			higher_order_prob = 1.0;
			for(int j = 0; j < M_; j++)
			{
				higher_order_prob = 1.0; 
				for(int k = 0; k < basesegmentcounts; k++)
				{
					curr_pix_index = baseSegmentIndexes[layerM][i][k];
					higher_order_prob = higher_order_prob * current_[curr_pix_index*M_+j]; 
				}
				h_norm[i*M_+j] = higher_order_prob; 
			}
		}

		double alpha = 0.5, maxcost, weight, costdata = 0.0; 
		start_loc = 0; //int start_loc = 0; 

		for(int i = 0; i < layerM; i++)
			start_loc = start_loc + segmentCount[i]; 

		start_loc = start_loc * M_; 
	
		for(int i = 0; i < segmentCount[layerM]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layerM][i];
		
			weight = 0.3 * basesegmentcounts; 
			maxcost = -weight * log(alpha); 
		
			for(int j = 0; j < basesegmentcounts; j++)
			{
				curr_pix_index = baseSegmentIndexes[layerM][i][j]; 
				for(int k = 0; k < M_; k++)
				{
					higher_order_prob = h_norm[i*M_+k]/(current_[curr_pix_index*M_+k]+0.0001);
					costdata = - weight * log(stats_potential[start_loc+k]); 
					higher_order[curr_pix_index*M_+k] += (ho_param1*costdata - ho_param2*higher_order_prob);
				}		
			}	
			start_loc = start_loc+M_; 
		}
		
		
		delete []h_norm; 	
	}
}

void DenseCRF::set_mem_init_higherorder()
{
	if(addHO)
	{
		for(int i = 0; i < N_*M_; i++)
			higher_order[i] = 0.0; 
	}
}


void DenseCRF::del_mem_higherorder()
{
	if(addHO)
	{
		delete []higher_order;
		delete []segmentimagedata;
		delete []temp_segmentimagedata;
		for(int i = 0; i < num_layers; i++)
		{
			for(int j = 0; j < segmentCount[i]; j++)
			{
				delete []baseSegmentIndexes[i][j];
			}
			delete []baseSegmentIndexes[i]; 
		}

		delete []baseSegmentIndexes; 

		for(int i = 0; i < num_layers; i++)
		{
			delete []baseSegmentCounts[i]; 
		}
		delete []baseSegmentCounts;
		delete[]segmentCount; 
		delete []stats_potential; 		
	}
	if(addCooc)
	{
		delete []unary_cooc_; 
		delete []pair_cooc_; 	
	}
	if(addDet)
	{
		delete []det_higher_order; 
		delete []det_baseSegmentCounts;
		delete []det_baseSegmentIndexes;
		delete []det_stats_potential; 
		delete []det_h_norm; 
	}
}

//subarna
void DenseCRF::del_mem_higherorderM(int num_images)
{
	if(addHO)
	{
		delete []higher_order;
		delete []segmentimagedata;
		delete []temp_segmentimagedata;

		//subarna
		for(int i = 0; i < num_layers*num_images; i++)
		{
			for(int j = 0; j < segmentCount[i]; j++)
			{
				delete []baseSegmentIndexes[i][j];
			}
			delete []baseSegmentIndexes[i]; 
		}

		delete []baseSegmentIndexes; 


		//subarna
		for(int i = 0; i < num_images*num_layers; i++)
		{
			delete []baseSegmentCounts[i]; 
		}
		delete []baseSegmentCounts;
		delete[]segmentCount; 
		delete []stats_potential; 		
	}
	if(addCooc)
	{
		delete []unary_cooc_; 
		delete []pair_cooc_; 	
	}
	if(addDet)
	{
		delete []det_higher_order; 
		delete []det_baseSegmentCounts;
		delete []det_baseSegmentIndexes;
		delete []det_stats_potential; 
		delete []det_h_norm; 
	}
}

void DenseCRF::add_higher_order()
{

}

// end higher order potential related stuffs


// start det 

void DenseCRF::set_hodet(int add_det)
{
	addDet = 0; 
	if(add_det) addDet = 1; 
}

void DenseCRF::det_ho_mem_init(int imagewidth, int imageheight, char *det_seg, char *det_bb, char *det_seg_ext, char *det_bb_ext, float detParam1, float detParam2)
{
	if(addDet)
	{
		det_seg_dir = det_seg; 
		det_bb_dir = det_bb; 
		det_segment_ext = det_seg_ext; 
		det_bbox_ext = det_bb_ext; 
		det_ho_width = imagewidth; 
		det_ho_height = imageheight; 
		det_higher_order = new float[det_ho_width*det_ho_height*M_];
		for(int i = 0; i < det_ho_width*det_ho_height*M_; i++)
		{
			det_higher_order[i] = 0.0; 
		}

		det_param1 = detParam1; 
		det_param2 = detParam2; 
	}	
}

void DenseCRF::det_readSegmentIndex(char *fileName)
{
	if(addDet)
	{
		char *file_name = new char[200]; 
		sprintf(file_name, "%s%s.%s",det_seg_dir, fileName, det_segment_ext);
		FILE *fin = fopen(file_name, "rb");
		fread(&det_segmentCount, sizeof(int), 1, fin);
		det_baseSegmentCounts = new int[det_segmentCount];
		det_baseSegmentIndexes = new int *[det_segmentCount];
		int **det_tempbaseSegmentIndexes = new int*[det_segmentCount];
		for(int i = 0; i < det_segmentCount; i++)
		{
			fread(&det_baseSegmentCounts[i], sizeof(int), 1, fin);
			det_baseSegmentIndexes[i] = new int[det_baseSegmentCounts[i]];
			det_tempbaseSegmentIndexes[i] = new int[det_baseSegmentCounts[i]];
			if(det_baseSegmentCounts[i] != 0) fread(det_tempbaseSegmentIndexes[i], sizeof(int), det_baseSegmentCounts[i], fin);
	
			for(int j = 0; j < det_baseSegmentCounts[i]; j++)
			{
				int temp_id1 = det_tempbaseSegmentIndexes[i][j], x1, y1; 
				x1 = temp_id1%det_ho_width; y1 = temp_id1/det_ho_width; 
				det_baseSegmentIndexes[i][j] = (det_ho_height-(y1+1))*det_ho_width+x1; 
			}
		}
		fclose(fin);

		sprintf(file_name, "%s%s.%s", det_bb_dir, fileName, det_bbox_ext);
		unsigned char det_segmented; int det_objects; 
		FILE *fin1 = fopen(file_name, "rb");

		fread(&det_segmented, sizeof(unsigned char), 1, fin1);
		fread(&det_objects, sizeof(int), 1, fin1);

		det_stats_potential = new int[det_objects];
		
		int x1, y1, x2, y2; unsigned char type; 

		for(int i = 0; i < det_objects; i++)
		{
			fread(&type, sizeof(unsigned char), 1, fin1); 
			fread(&x1, sizeof(int), 1, fin1); fread(&y1, sizeof(int), 1, fin1); fread(&x2, sizeof(int), 1, fin1); fread(&y2, sizeof(int), 1, fin1);
			det_stats_potential[i] = (int)type; 
		}	
	
		det_resp = new double[det_objects];
		if(det_objects > 0) fread(det_resp, sizeof(double), det_objects, fin1);
		fclose(fin1);	

		det_h_norm = new double[det_segmentCount*M_];
	}
}

void DenseCRF::det_calculate_higherorder_pot()
{
	if(addDet)
	{
		float norm_val = 0.0; 
	
		int basesegmentcounts = 0; 
		int curr_pix_label = 0, curr_pix_index; //int x, y; 
		//int neigh_pix_index, neigh_pix_label; 
	
		double higher_order_prob; 
	
		for(int i = 0; i < det_segmentCount; i++)
			for(int j = 0; j < M_; j++)
				det_h_norm[i*M_+j] = 1.0; 
	
		for(int i = 0; i < det_segmentCount; i++)
		{
			basesegmentcounts = det_baseSegmentCounts[i];
			higher_order_prob = 1.0;
			for(int j = 0; j < M_; j++)
			{
				higher_order_prob = 1.0; 
				for(int k = 0; k < basesegmentcounts; k++)
				{
					curr_pix_index = det_baseSegmentIndexes[i][k];
					higher_order_prob = higher_order_prob * current_[curr_pix_index*M_+j]; 
				}
				det_h_norm[i*M_+j] = higher_order_prob; 
			}
		}


		for(int i = 0; i < det_ho_width*det_ho_height*M_; i++)
		{
			det_higher_order[i] = 0.0; 
		}

		double alpha = 0.5, maxcost, weight, costdata = 0.0;

		for(int i = 0; i < det_segmentCount; i++)
		{
			basesegmentcounts = det_baseSegmentCounts[i];
		
			weight = 0.3 * basesegmentcounts; 
			maxcost = -weight * log(alpha); 
			costdata = 6.0 * basesegmentcounts * (det_resp[i]+1.2);

			if(costdata < 0) costdata = 0; 

			for(int j = 0; j < basesegmentcounts; j++)
			{
				curr_pix_index = det_baseSegmentIndexes[i][j]; 
	
				for(int k = 0; k < M_; k++)
				{
					higher_order_prob = det_h_norm[i*M_+k]/(current_[curr_pix_index*M_+k]+0.0001);
					if(k != det_stats_potential[i])
					{
						det_higher_order[curr_pix_index*M_+k] += det_param1*costdata-det_param2*higher_order_prob; 
					}				
				}
			}
		}
	}
}

// end det