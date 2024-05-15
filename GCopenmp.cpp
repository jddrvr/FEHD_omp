#include "GCopenmp.h"
#include "dataContainers.h"
#include "kernelsCPU.h"
#include "utility.h"
#include <cblas.h>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_float_real(z) z.real()
#define lapack_complex_float_imag(z) z.imag()
#include <lapacke.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>


void granger(float *AR, std::vector<float> angleArray, std::vector<float> &GCvals, paramContainer params,
	     int numComps, std::vector<float> Q, float *rotatedModels,
	     float *workArray, std::complex<float> *Tf,
	     std::complex<float> *Swhole, std::complex<float> *tmp,
	     std::complex<float> *Spartial, std::complex<float> *d_wholeSpec, 
	     float *det_whole, float *det_partial,float *eigValsWhole,
	     float *eigValsPartial,std::complex<float> *modelTranspose,
	     std::complex<float> *tmpVector, std::complex<float> *work,
	     float *rwork) 
{

  // Create the Q matrices, given the angleArray passed.
  float alpha=1.0f;
  float beta = 0.0f;

  std::fill(GCvals.begin(),GCvals.end(),0.0);
  
  //printf("Entered \n");
  //std::vector<float> Q(params.numParticles*numComps*numComps,0);

  float deltaT = 1.0f/(float)(params.sampRate);
  float sinVal;
  float cosVal;

  float Qcol1;
  float Qcol2;
  float argmtBASE = -2.0*M_PI*deltaT;
  std::vector<float> freqVal(params.numFreqs,0.0);
  for(int freq=0;freq < params.numFreqs;freq++)
    {
      freqVal[freq] = (params.freqHi-params.freqLo)/((float)(params.numFreqs-1))*(float)(freq)+params.freqLo;
      //printf("%f ",freqVal[freq]);
    }
  //std::vector<int> lagList(params.lagList);
  std::sort(params.lagList.begin(),params.lagList.end());
  std::complex<float> alphaC(1.0,0.0);
  std::complex<float> betaC(0.0,0.0);
  std::complex<float> alphaC2(-1.0,0.0);
  std::complex<float> betaC2(1.0,0.0);
  std::complex<float> argmt;



  float ratio;

  std::fill(Q.begin(),Q.end(),0.0f);

  int thread_num;

  int lwork = (numComps-1)*(numComps-1);
  int rsize = 3*(numComps-1)-2;

#pragma omp parallel default(shared) private(sinVal,cosVal,Qcol1,Qcol2,argmt,ratio,thread_num) //num_threads(24)
  {    
#pragma omp for schedule(dynamic,10) 
    for(int particle=0; particle<params.numParticles; particle++)
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

	// Give a sample of Q

	
	cblas_sgemm(CblasColMajor, CblasTrans,CblasNoTrans,
		    numComps*params.numLags,numComps,numComps,
		    alpha,AR,numComps,Q.data()+particle*numComps*numComps,numComps,
		    beta,rotatedModels+particle*params.numLags*numComps*numComps,params.numLags*numComps);


	// Transpose the individual matrices
	for(int lag=0;lag<params.numLags;lag++)
	  for(int col=0;col<numComps;col++)
	    for(int row=0;row<numComps;row++)
	      workArray[particle*params.numLags*numComps*numComps+lag*numComps+col*params.numLags*numComps+row] =
		rotatedModels[particle*params.numLags*numComps*numComps+lag*numComps+row*params.numLags*numComps+col];

	// Multiply by the second Q
	cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		    numComps*params.numLags,numComps,numComps,
		    alpha,
		    workArray+particle*params.numLags*numComps*numComps,numComps*params.numLags,
		    Q.data()+particle*numComps*numComps,numComps,
		    beta,
		    rotatedModels+particle*params.numLags*numComps*numComps,numComps*params.numLags);

    // For each particle, numLags matrices need to be transposed.
	for(int col = 0;col<numComps;col++)
	  for(int row=0;row<params.numLags*numComps;row++)
	    modelTranspose[particle*numComps+row*params.numParticles*numComps+col]=
			   std::complex<float>(rotatedModels[particle*params.numLags*numComps*numComps+
							     col*params.numLags*numComps+row],0.0);
      

	for(int freq=0;freq<params.numFreqs;freq++)
	  {
	    for(int col=0;col<numComps;col++)
	      for(int row=0;row<numComps;row++)
		tmpVector[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+col*numComps+row]=std::complex<float>(0,0);
	    
	    for(int lag=0;lag<params.numLags;lag++)
	      {
		argmt = std::exp(std::complex<float>(0,argmtBASE*((float)(params.lagList[lag]))*freqVal[freq]));
		
		for(int col=0;col<numComps;col++)
		  cblas_caxpy(numComps,&argmt,modelTranspose+particle*numComps+lag*params.numParticles*numComps*numComps+
			col*params.numParticles*numComps,1,
			      tmpVector+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+col*numComps,1);
	    
	      }
	  }



	for(int freq=0;freq<params.numFreqs;freq++)
	  for(int col=0;col<numComps;col++)
	    for(int row=0;row<numComps;row++)
	      {
		Tf[(particle*params.numFreqs+freq)*numComps*numComps+row*numComps+col]=
		  -tmpVector[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+col*numComps+row];
		if(col==row)
		  Tf[(particle*params.numFreqs+freq)*numComps*numComps+col*numComps+row]+=std::complex<float>(1.0,0.0);
	      }
	// Tf is formed here. Print it here, for the zeroth particle and frequnecy
	//if(particle==0)
	//  for(int row=0;row<numComps;row++)
	//    {
	//      for(int col=0;col<numComps;col++)
	//	printf("%f+i*%f ",Tf[col*numComps+row].real(),Tf[col*numComps+row].imag());
	//      printf("\n");
	//    }
      
	for(int freq=0;freq<params.numFreqs;freq++)
	  cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,
		      numComps,numComps,numComps,
		      &alphaC,
		      Tf+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps,
		      Tf+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps,
		      &betaC,
		      Swhole+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps);
	
	for(int freq=0;freq<params.numFreqs;freq++)
	    for(int row=0;row<numComps-1;row++)
	      {
		Swhole[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]=
		  Swhole[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]/
		  Swhole[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+(numComps+1)*(numComps-1)];
	      }

      
	std::copy(Swhole+(particle*params.numFreqs)*numComps*numComps,
		  Swhole+(particle*params.numFreqs+params.numFreqs)*numComps*numComps,
		  tmp+(particle*params.numFreqs)*numComps*numComps);

	for(int freq=0;freq<params.numFreqs;freq++)
	  cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		      numComps-1,numComps-1,1,
		      &alphaC2,
		      Swhole+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+(numComps-1)*numComps,numComps,
		      Swhole+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+numComps-1,numComps,
		      &betaC2,
		      tmp+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps);

	std::copy(tmp+(particle*params.numFreqs)*numComps*numComps,
		  tmp+(particle*params.numFreqs+params.numFreqs)*numComps*numComps,
		  Swhole+(particle*params.numFreqs)*numComps*numComps);
	
	for(int freq=0;freq<params.numFreqs;freq++)
	  for(int row=0;row<numComps-1;row++)
	    {
	      Tf[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]=
		Tf[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+numComps*(numComps-1)+row]/
		Tf[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+(numComps+1)*(numComps-1)];
	    }


	std::copy(Tf+(particle*params.numFreqs)*numComps*numComps,
		  Tf+(particle*params.numFreqs+params.numFreqs)*numComps*numComps,
		  tmp+(particle*params.numFreqs)*numComps*numComps);
      
	for(int freq=0;freq<params.numFreqs;freq++)
	  cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		      numComps-1,numComps-1,1,
		      &alphaC2,
		      Tf+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+(numComps-1)*numComps,numComps,
		      Tf+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+(numComps-1),numComps,
		      &betaC2,
		      tmp+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps);
	
	for(int freq=0;freq<params.numFreqs;freq++)
	  cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,
		      numComps-1,numComps-1,numComps-1,
		      &alphaC,
		      tmp+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps,
		      tmp+particle*params.numFreqs*numComps*numComps+freq*numComps*numComps,numComps,
		      &betaC,
		      Spartial+particle*params.numFreqs*(numComps-1)*(numComps-1)+freq*(numComps-1)*(numComps-1),numComps-1);
	// There might be a better way to do this shrinkage.
	for(int freq=0;freq<params.numFreqs;freq++)
	  for(int col=0;col<numComps-1;col++)
	    for(int row=0;row<numComps-1;row++)
	      d_wholeSpec[particle*params.numFreqs*(numComps-1)*(numComps-1)+freq*(numComps-1)*(numComps-1)+col*(numComps-1)+row]=
		Swhole[particle*params.numFreqs*numComps*numComps+freq*numComps*numComps+col*numComps+row];

	
	thread_num=omp_get_thread_num();
		
	for(int freq=0;freq<params.numFreqs;freq++)      
	  {
	    auto info = LAPACKE_cheev_work(LAPACK_COL_MAJOR,'N','U',numComps-1,d_wholeSpec+
					   particle*params.numFreqs*(numComps-1)*
					   (numComps-1)+freq*(numComps-1)*(numComps-1),
					   numComps-1,
					   eigValsWhole+particle*params.numFreqs*(numComps-1)+freq*(numComps-1),
					   work+thread_num*lwork,lwork,rwork+thread_num*rsize);
	    
	    LAPACKE_cheev_work(LAPACK_COL_MAJOR,'N','U',numComps-1,Spartial+
	    		  particle*params.numFreqs*(numComps-1)*
	    		  (numComps-1)+freq*(numComps-1)*(numComps-1),
	    		  numComps-1,
	    		  eigValsPartial+particle*params.numFreqs*(numComps-1)+freq*(numComps-1),
			  work+thread_num*lwork,lwork,rwork+thread_num*rsize);
	  }
		
	for(int freq=0;freq<params.numFreqs;freq++)
	  {
	    det_partial[particle*params.numFreqs+freq]=1.0;
	    det_whole[particle*params.numFreqs+freq]=1.0;
	    
	    for(int comp=0;comp<numComps-1;comp++)
	      {
		det_partial[particle*params.numFreqs+freq]=det_partial[particle*params.numFreqs+freq]*
		  eigValsPartial[particle*params.numFreqs*(numComps-1)+
				 freq*(numComps-1)+comp];
		det_whole[particle*params.numFreqs+freq]=det_whole[particle*params.numFreqs+freq]*
		  eigValsWhole[particle*params.numFreqs*(numComps-1)+
			       freq*(numComps-1)+comp];
	      }
	  }

	GCvals[particle] = 0.0;
	for(int freq=0;freq<params.numFreqs;freq++)
	  {
	    ratio = det_partial[particle*params.numFreqs+freq]/det_whole[particle*params.numFreqs+freq];
	    //if(ratio > 1.0f)
	    GCvals[particle] = GCvals[particle]+logf(ratio);
	    //else if(ratio < 0.9f)
	    //{
	    //  printf("determinants are real close to zero, and the ratio is off. \n");
	    //  printf("idx: %i, upper determinant: %e, lower determinant: %e \n",particle*params.numFreqs+freq,
	    //	     det_partial[particle*params.numFreqs+freq], det_whole[particle*params.numFreqs+freq]);
	    //}
	  }      
      }
  }

  return;
}
