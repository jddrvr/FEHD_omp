#ifndef GCOPENMP_H
#define GCOPENMP_H
#include <vector>
#include <complex>
#include "utility.h"
#include "dataContainers.h"

void granger(float *AR, std::vector<float> angleArray, std::vector<float> &GCvals, paramContainer params,
	     int numComps, std::vector<float> Q, float *rotatedModels,
	     float *workArray, std::complex<float> *Tf,
	     std::complex<float> *Swhole, std::complex<float> *tmp,
	     std::complex<float> *Spartial, std::complex<float> *d_wholeSpec, 
	     float *det_whole, float *det_partial,float *eigValsWhole,
	     float *eigValsPartial,std::complex<float> *modelTranspose,
	     std::complex<float> *tmpVector, std::complex<float> *work,
	     float *rwork); 
#endif

