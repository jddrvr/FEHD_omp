#ifndef FEHD_H
#define FEHD_H
#include <complex>
#include "dataContainers.h"
#include "utility.h"
// Computes the gradient for multiple rotations of an AR model.
void compGradient(float *AR, std::vector<float> &gradient, std::vector<float> angleArray, std::vector<float> GCvals,
		  paramContainer params, int numComps, std::vector<float> Q, float *rotatedModels,
		  float *workArray, std::complex<float> *Tf,
		  std::complex<float> *Swhole, std::complex<float> *tmp,
		  std::complex<float> *Spartial, std::complex<float> *d_wholeSpec, 
		  float *det_whole, float *det_partial,float *eigValsWhole,
		  float *eigValsPartial,std::complex<float> *modelTranspose,
		  std::complex<float> *tmpVector,std::complex<float> *work,
		  float *rwork); 
void runFEHD(dataList, std::vector<float> &, paramContainer);

void runFEHDstep(std::vector<float> &, matrix &, dataList, paramContainer, int);

void singleQ(std::vector<float> &, std::vector<float>);

//void multMat(float *,float *,int,int);



#endif
