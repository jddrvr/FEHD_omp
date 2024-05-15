#include "kernelsCPU.h"
#include <vector>
#include <math.h>
#include <complex>
#include <cblas.h>
#include <cmath>
#include <algorithm>
void generateRotationMatrices(std::vector<float> angleArray, std::vector<float> &Q, int numComps, int numParticles)
{
  float sinVal;
  float cosVal;

  float Qcol1;
  float Qcol2;

  std::fill(Q.begin(),Q.end(),0.0f);
 
    
#pragma omp parallel for default(shared) schedule(dynamic,10) private(sinVal,cosVal,Qcol1,Qcol2)
  for(int particle=0;particle<numParticles;particle++)
    {
      // initialize to zero
      
      for(int row=0;row<numComps;row++)
	Q[particle*numComps*numComps+row*numComps+row] = 1.0f;

      // Cycle through the angles
      for(int angle=numComps-2; angle>=0; angle--)
	{
	  sinVal = sinf(angleArray[particle*(numComps-1)+angle]);
	  cosVal = cosf(angleArray[particle*(numComps-1)+angle]);

	  for(int k=0;k<numComps;k++)
	    {
	      Qcol1 = Q[particle*numComps*numComps+k*numComps+angle];
	      Qcol2 = Q[particle*numComps*numComps+k*numComps+numComps-1];
	      Q[particle*numComps*numComps+angle+k*numComps] = cosVal*Qcol1+sinVal*Qcol2;
	      Q[particle*numComps*numComps+(numComps-1)+k*numComps] = -sinVal*Qcol1+cosVal*Qcol2;
	    }
	}
    }
  return;
}

void transposeBlockMatrices(std::vector<float> &outMatrix, std::vector<float> inMatrix, int numComps, int numParticles, int numLags)
{
#pragma omp parallel for default(shared) schedule(dynamic,10)
  for(int particle=0;particle<numParticles;particle++)
    for(int lag=0;lag<numLags;lag++)
	for(int col=0;col<numComps;col++)
	  for(int row=0;row<numComps;row++)
	    outMatrix[particle*numLags*numComps*numComps+lag*numComps+col*numLags*numComps+row] =
	      inMatrix[particle*numLags*numComps*numComps+lag*numComps+row*numLags*numComps+col];


  return;
}

void compTransferFunc(std::vector<float> models, std::vector<std::complex<float>> &Tf, int numComps, paramContainer params)
{

  
  float deltaT = 1.0f/(float)(params.sampRate);

  // The calculations could be simplified greatly here, however I prioritized how I store the
  // matrices over this. This is probably just another way of saying I am being lazy.

  // I think that the issue is that there are too many private variables here.
  // First idea is to vectorize them and then just lett it share.
  // Evidence, from above, is that this is not the problem.
  float argmtBASE = -2.0*M_PI*deltaT;
  std::vector<float> freqVal(params.numFreqs,0.0);
  for(int freq=0;freq < params.numFreqs;freq++)
    {
      freqVal[freq] = (params.freqHi-params.freqLo)/((float)(params.numFreqs-1))*(float)(freq)+params.freqLo;
      //printf("%f ",freqVal[freq]);
    }
  std::vector<int> lagList(params.lagList);
  std::sort(lagList.begin(),lagList.end());

  std::vector<std::complex<float>> modelTranspose(models.size(),std::complex<float>(0.0,0.0));

#pragma omp parallel for default(shared) schedule(dynamic,10)
  for(int col = 0;col<params.numParticles*numComps;col++)
    for(int row=0;row<params.numLags*numComps;row++)
      modelTranspose[row*params.numParticles*numComps+col]=std::complex<float>(models[col*numComps*params.numLags+row],0.0);

      
  // The size of this should be numFreqs*numParticles*numComps*numComps
  std::vector<std::complex<float>> tmpVector(params.numParticles*numComps*numComps*params.numFreqs,std::complex<float>(0.0,0.0));
  std::complex<float> argmt;
#pragma omp parallel for default(shared) schedule(dynamic,1) private(argmt)
  for(int freq=0;freq<params.numFreqs;freq++)
    {
      std::fill(tmpVector.data()+freq*params.numParticles*numComps*numComps,
		tmpVector.data()+(freq+1)*params.numParticles*numComps*numComps,0.0);
      //printf("Round and round \n");
      for(int lag=0;lag<params.numLags;lag++)
	{
	  argmt = std::exp(std::complex<float>(0,argmtBASE*((float)(lagList[lag]))*freqVal[freq]));
	  
	  cblas_caxpy(params.numParticles*numComps*numComps,&argmt,
	  	      modelTranspose.data()+lag*params.numParticles*numComps*numComps,1,
	  	      tmpVector.data()+freq*params.numParticles*numComps*numComps,1);
	}
    }

  
    // Transpose back, and assemble into final form.
#pragma omp parallel for default(shared) schedule(dynamic,10)
  for(int particle=0;particle<params.numParticles;particle++)
    for(int freq=0;freq<params.numFreqs;freq++)
      for(int col=0;col<numComps;col++)
	for(int row=0;row<numComps;row++)
	  {
	    Tf[(particle*params.numFreqs+freq)*numComps*numComps+col*numComps+row]=
	      tmpVector[freq*params.numParticles*numComps*numComps+particle*numComps+row*params.numParticles*numComps+col];
	    if(col==row)
	      Tf[(particle*params.numFreqs+freq)*numComps*numComps+col*numComps+row]+=std::complex<float>(1.0,0.0);
	  }
  return;
  
}

void scale_columns(std::vector<std::complex<float>> &TF, int numComps, int numParticles, int numFreqs)
{
  int particle;
  int freq;
  int row;
  
#pragma omp parallel for default(shared) schedule(dynamic,10) private(particle,freq,row)
  for(particle=0;particle<numParticles;particle++)
    {
      for(freq=0;freq<numFreqs;freq++)
	{
	  for(row=0;row<numComps-1;row++)
	    {
	      TF[particle*numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]=
		TF[particle*numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]/
		TF[particle*numFreqs*numComps*numComps+freq*numComps*numComps+(numComps+1)*(numComps-1)];
	    }
	}
    }

  return;
}

void shrinkArrays(std::vector<std::complex<float>> S1, std::vector<std::complex<float>> &S2, int numComps, int numParticles, int numFreqs)
{
  int particle;
  int freq;
  int col;
  int row;

#pragma omp parallel for default(shared) schedule(dynamic,10) private(particle,freq,col,row)
  for(particle=0;particle<numParticles;particle++)    
    for(freq=0;freq<numFreqs;freq++)
      for(col=0;col<numComps-1;col++)
	for(row=0;row<numComps-1;row++)
	  S2[particle*numFreqs*(numComps-1)*(numComps-1)+freq*(numComps-1)*(numComps-1)+col*(numComps-1)+row]=
	    S1[particle*numFreqs*numComps*numComps+freq*numComps*numComps+col*numComps+row];

  return;
}

void prodEigs(std::vector<float> W, std::vector<float> &det, int numComps, int numParticles,int numFreqs)
{
  int particle;
  int freq;
  int comp;

#pragma omp parallel for default(shared) schedule(dynamic,10) private(particle,freq,comp)
  for(particle=0;particle<numParticles;particle++)
    for(freq=0;freq<numFreqs;freq++)
      for(comp=0;comp<numComps;comp++)
	det[particle*numFreqs+freq]=det[particle*numFreqs+freq]*W[particle*numFreqs*numComps+freq*numComps+comp];
}

void det2GC(std::vector<float> detPartial,std::vector<float> detWhole, std::vector<float> &GC, int numParticles, int numFreqs)
{
  int particle;
  int freq;
  float ratio;
#pragma omp parallel for default(shared) schedule(dynamic,10) private(particle,freq,ratio)
  for(particle=0;particle<numParticles;particle++)
    {
      for(freq=0;freq<numFreqs;freq++)
	{
	  ratio = detPartial[particle*numFreqs+freq]/detWhole[particle*numFreqs+freq];
	  if(ratio > 1.0f)
	    GC[particle] = GC[particle]+logf(ratio);
	  else if(ratio < 0.9f)
	    {
	      printf("determinants are real close to zero, and the ratio is off. \n");
	      printf("idx: %i, upper determinant: %e, lower determinant: %e \n",particle*numFreqs+freq,
		     detPartial[particle*numFreqs+freq], detWhole[particle*numFreqs+freq]);
	    }
	}
      
    }
      
  return;
}
