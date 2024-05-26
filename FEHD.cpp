#include <iostream>
#include <complex>
#include "FEHD.h"
#include "dataContainers.h"
#include "mkAR.h"
#include <algorithm>
#include <math.h>
#include <fstream>
#include "GCopenmp.h"
#include "kernelsCPU.h"
#include "utility.h"
#include <chrono>
#include <cblas.h>
#include "timeSeriesOPs.h"
#include <omp.h>

// Compute the gradient of the GC objective function at a point in angle-space.
void compGradient(float *AR, std::vector<float> &gradient, std::vector<float> angleArray, std::vector<float> GCvals,
		  paramContainer params, int numComps, std::vector<float> Q, float *rotatedModels,
		  float *workArray, std::complex<float> *Tf,
		  std::complex<float> *Swhole, std::complex<float> *tmp,
		  std::complex<float> *Spartial, std::complex<float> *d_wholeSpec, 
		  float *det_whole, float *det_partial,float *eigValsWhole,
		  float *eigValsPartial,std::complex<float> *modelTranspose,
		  std::complex<float> *tmpVector,std::complex<float> *work,
		  float *rwork)
{
  
  const int numVars = numComps-1;
  const float  h_val = 0.001f; // This is for the gradient spacing.
  
  std::vector<float> angle(angleArray); // Copies the value. At this time I am not sure why.
  std::vector<float> GCvalsUTIL(params.numParticles,0); // Output at each angle perturbation.
  // The gradient is computed by a difference in each angle of rotation.

  for(int angleArrayIndex=0;angleArrayIndex<numVars;angleArrayIndex++)
    {
      for(int particle=0;particle<params.numParticles;particle++)
	{
	  angle[particle*numVars+angleArrayIndex] += h_val;
	}
      
      granger(AR,angle,GCvalsUTIL,params,numComps,Q,rotatedModels,workArray,Tf,Swhole,tmp,Spartial,
	      d_wholeSpec,
	      det_whole,det_partial,eigValsWhole,eigValsPartial,modelTranspose,
	      tmpVector,work,rwork);
      
      for(int particle=0;particle<params.numParticles;particle++)
	{
	  gradient[particle*numVars+angleArrayIndex] = (GCvalsUTIL[particle]-// This is saxpy
							GCvals[particle])/h_val;
	  angle[particle*numVars+angleArrayIndex] -= h_val;
	}      
    }  
}


void runFEHD(dataList dataArray, std::vector<float> &Lmat, paramContainer params)
{
  // Set the parameters for sgemm. 
  float alpha=1.0f;
  float beta=0.0f;

  
  std::vector<float> bestAngle;
  matrix Rdecor; // Another example - it's the pca function.
  std::vector<float> Q; // The rotation matrix, resized at each step.
  std::vector<float> T(params.numPCs*params.numPCs,0); // The "work" transformation
  std::vector<float> oneArrayData; // Holds the data without epoch boundaries.
  std::vector<float> transformedData; // Holds the new data without epoch boundaries.  
  std::vector<float> newTrans(params.numPCs*params.numChannels); // Another worker.

  for(int numComps = params.numPCs;numComps>1;numComps--)
    {
      
      Rdecor.elements.clear();
      bestAngle.resize(numComps-1);
     
      runFEHDstep(bestAngle, Rdecor, dataArray, params, numComps);
      //printf("got here \n");
      
      Q.resize(numComps*numComps);

      singleQ(Q,bestAngle);

      // Compute the total transformation for this step
      std::fill(T.begin(),T.end(), 0.0f);

      for(int i=numComps;i<params.numPCs;i++)
	T[i*params.numPCs+i] = 1.0f;

      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		  numComps,numComps,numComps,
		  alpha,Q.data(),numComps,
		  Rdecor.elements.data(),numComps,
		  beta,T.data(),params.numPCs);

      oneArrayData.resize(numComps*params.numEpochs*params.epochPts);
      convertDataListToRawArray(dataArray, oneArrayData.data());

      // Transform the data

      
      transformedData.resize(numComps*params.numEpochs*params.epochPts);

      cblas_sgemm(CblasColMajor, CblasNoTrans,CblasNoTrans,
		  numComps,params.numEpochs*params.epochPts,numComps,
		  alpha,T.data(),params.numPCs,
		  oneArrayData.data(), numComps,
		  beta,transformedData.data(), numComps);

      // Convert it back into dataArray

      convertRawArrayToDataList(transformedData.data(), dataArray, numComps, params.epochPts, params.numEpochs);

      // Remove the last component
      removeComponent(dataArray, numComps-1);

      // Update the transformation

      cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		  params.numPCs,params.numChannels,params.numPCs,
		  alpha, T.data(),params.numPCs,
		  Lmat.data(),params.numPCs,
		  beta, newTrans.data(), params.numPCs);

      Lmat = newTrans;
    }
  
  return;
}

// Runs a step of the algorithm. Randomly initializes the rotation angles,
// rotates the AR model and computes the Granger Causality of the last component
// on those above. The component found to have the minimum influence
// is subsequently removed from the time series. 
void runFEHDstep(std::vector<float> &bestAngle, matrix &L, dataList dataArray ,paramContainer params,int numComps)
{

  // Create random initial conditions
  std::vector<float> angleArray;

  srand((unsigned)time(0));

  for(int indx=0;indx<params.numParticles*(numComps-1);indx++)
    angleArray.push_back((float)(rand()%314-157)/100.0f);

  for(int indx=0;indx<numComps-1;indx++)
    angleArray[indx]=0;
  
  ARmodel A;
  dataList residuals;

  std::sort(params.lagList.begin(),params.lagList.end());

  
  // Assume you are dealing with a sorted lagList
  mkAR(dataArray, params.lagList, A, residuals);

  dataList ortho_residuals;

  orthonormalizeR(residuals, ortho_residuals, L);
  
  rotate_model(A, L);

  std::vector<float> AR(params.numLags*numComps*numComps);
  for(int lag=0;lag<params.numLags;lag++)
    for(int row=0;row<numComps;row++)
      for(int col=0;col<numComps;col++)
	AR[lag*numComps*numComps+col*numComps+row] = A.lagMatrices[lag].elements[col*numComps+row];
  
  // Number of iterations
  int numIts = 50;
  std::vector<float> h = {0.01f, 1.0f, 10.0f};
  // Declare a couple of arrays to hold the GCs

  std::vector<float> GCvals(params.numParticles,0);
  std::vector<float> GCvals1(params.numParticles,0);
  std::vector<float> GCvals2(params.numParticles,0);
  std::vector<float> GCvals3(params.numParticles,0);
  
  // And some more angle arrays

  std::vector<float> angleArray1(params.numParticles*(numComps-1),0);
  std::vector<float> angleArray2(params.numParticles*(numComps-1),0);
  std::vector<float> angleArray3(params.numParticles*(numComps-1),0);

  // And a gradient

  std::vector<float> gradient(params.numParticles*(numComps-1),0);

  // I want to declare all of the arrays here, before entering for the first time.
  // This might make using the recycler difficult, as it will require a dynamically
  // sized version of all of this. NOT TRUE. Just set the particles parameter, it will ignore extra array.

  // Create a structure that contains all of the following vectors
  //GCarrays GC;
  //sizeGCArrays(GC);

  // Determine system memory here.
    
  std::vector<float> Q(numComps*numComps*params.numParticles);
  std::vector<std::complex<float>> Swhole(numComps*numComps*params.numFreqs*params.numParticles);
  std::vector<std::complex<float>> tmp(numComps*numComps*params.numFreqs*params.numParticles);
  std::vector<std::complex<float>> modelTranspose(params.numParticles*params.numLags*numComps*numComps);
  std::vector<std::complex<float>> tmpVector(params.numParticles*numComps*numComps*params.numFreqs);
  std::vector<float> rotatedModels(params.numParticles*params.numLags*numComps*numComps);
  std::vector<float> workArray(params.numParticles*params.numLags*numComps*numComps);
  std::vector<std::complex<float>> Tf(params.numParticles*params.numFreqs*numComps*numComps);
  std::vector<std::complex<float>> Spartial(params.numParticles*params.numFreqs*(numComps-1)*(numComps-1));
  std::vector<float> eigValsWhole(params.numParticles*params.numFreqs*(numComps-1));
  std::vector<float> eigValsPartial(params.numParticles*params.numFreqs*(numComps-1));
  std::vector<std::complex<float>> d_wholeSpec(params.numParticles*params.numFreqs*(numComps-1)*(numComps-1));
  std::vector<float> det_partial(params.numParticles*params.numFreqs,1);
  std::vector<float> det_whole(params.numParticles*params.numFreqs,1);
  int num_threads=omp_get_max_threads();
  int lwork=8*(numComps-1)*(numComps-1);
  std::vector<std::complex<float>> work(num_threads*lwork);
  std::vector<float> rwork(num_threads*(3*(numComps-1)-2));
 
  granger(AR.data(),angleArray,GCvals,params,numComps,Q,rotatedModels.data(),workArray.data(),Tf.data(),
	  Swhole.data(),tmp.data(),Spartial.data(),d_wholeSpec.data(),det_whole.data(),det_partial.data(),
	  eigValsWhole.data(),eigValsPartial.data(),modelTranspose.data(),tmpVector.data(),work.data(),
	  rwork.data());
 
  //printf("%i \n", num_threads);
  std::vector<float> candidates(4,0);
  int minIndx;
  for(int iter=0;iter<numIts;iter++)
    {
      compGradient(AR.data(), gradient, angleArray, GCvals, params,
		   numComps, Q, rotatedModels.data(), workArray.data(), Tf.data(),
		   Swhole.data(), tmp.data(), Spartial.data(), d_wholeSpec.data(),
		   det_whole.data(), det_partial.data(), eigValsWhole.data(),
		   eigValsPartial.data(), modelTranspose.data(), tmpVector.data(),
		   work.data(),rwork.data());
      
      angleArray1 = angleArray;
      angleArray2 = angleArray;
      angleArray3 = angleArray;
      cblas_saxpy(params.numParticles*(numComps-1),-h[0],gradient.data(),1,angleArray1.data(),1);
      cblas_saxpy(params.numParticles*(numComps-1),-h[1],gradient.data(),1,angleArray2.data(),1);
      cblas_saxpy(params.numParticles*(numComps-1),-h[2],gradient.data(),1,angleArray3.data(),1);
      
      // Evaluate GC for each of the shifted angle arrays
      granger(AR.data(),angleArray1,GCvals1,params,numComps,Q,rotatedModels.data(),
	      workArray.data(),Tf.data(),Swhole.data(),tmp.data(),Spartial.data(),
	      d_wholeSpec.data(),det_whole.data(),det_partial.data(),eigValsWhole.data(),
	      eigValsPartial.data(),modelTranspose.data(),tmpVector.data(),work.data(),rwork.data());
      granger(AR.data(),angleArray2,GCvals2,params,numComps,Q,rotatedModels.data(),
	      workArray.data(),Tf.data(),Swhole.data(),tmp.data(),Spartial.data(),
	      d_wholeSpec.data(),det_whole.data(),det_partial.data(),eigValsWhole.data(),
	      eigValsPartial.data(),modelTranspose.data(),tmpVector.data(),work.data(),rwork.data());
      granger(AR.data(),angleArray3,GCvals3,params,numComps,Q,rotatedModels.data(),
	      workArray.data(),Tf.data(),Swhole.data(),tmp.data(),Spartial.data(),
	      d_wholeSpec.data(),det_whole.data(),det_partial.data(),eigValsWhole.data(),
	      eigValsPartial.data(),modelTranspose.data(),tmpVector.data(),work.data(),rwork.data());
      
      // This loop determines the minimum of four values for each particle.
      // This loop can be made much more elegant. I don't know that it will work any better, but this is just ugly.
      for(int particle=0;particle<params.numParticles;particle++)
	{
	  candidates[0] = GCvals[particle];
	  candidates[1] = GCvals1[particle];
	  candidates[2] = GCvals2[particle];
	  candidates[3] = GCvals3[particle];

	  minIndx = std::distance(candidates.begin(),min_element(candidates.begin(),candidates.end()));

	  if(minIndx == 1)
	    {
	      GCvals[particle] = GCvals1[particle];
	      std::copy(angleArray1.data()+particle*(numComps-1),angleArray1.data()+
			particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	  if(minIndx == 2)
	    {
	      GCvals[particle] = GCvals2[particle];
	      std::copy(angleArray2.data()+particle*(numComps-1),angleArray2.data()+
			particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	  if(minIndx == 3)
	    {
	      GCvals[particle] = GCvals3[particle];
	      std::copy(angleArray3.data()+particle*(numComps-1),angleArray3.data()+
			particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	}
      
      //if(verbose)
      	printf("iteration = %i, particle = %li, value = %e \n",
             iter,std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin(),
             GCvals[std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin()]);
  
    }

  // Return the best angle.

  long unsigned int indexVal = std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin();

  //printf("%li \n",indexVal);

  std::copy(angleArray.data()+indexVal*(numComps-1),angleArray.data()+indexVal*(numComps-1)+
	    numComps-1,bestAngle.begin());

  return;
 
}



void singleQ(std::vector<float> &Q, std::vector<float> angle)
{

  int numVars = angle.size();
  
  float Qcol1;
  float Qcol2;

  float sinVal;
  float cosVal;

  // Create an identity matrix.
  for(int indx=0;indx<(numVars+1)*(numVars+1);indx++)
    {
      Q[indx] = 0;
    }
  for(int row=0;row<(numVars+1);row++)
    {
      Q[row*(numVars+1)+row] = 1.0f;
    }

  // Multiply the individual rotations together. 
  for(int varIndx=0; varIndx<numVars; varIndx++) // Cycle through the angles. 
    {
      sinVal = sinf(angle[varIndx]);// Assign the cos and sin to variables.
      cosVal = cosf(angle[varIndx]);

      for(int k=0;k<numVars+1;k++) // Do the matrix multiplication 
	{
	  Qcol1 = Q[k*(numVars+1)+varIndx]; //Q(particle,row i, column k
	  Qcol2 = Q[k*(numVars+1)+numVars];//(particle,row M-1,column k
	  Q[k*(numVars+1)+varIndx] = cosVal*Qcol1-sinVal*Qcol2;
	  Q[k*(numVars+1)+numVars] = sinVal*Qcol1+cosVal*Qcol2;
	}

    }
  
  return;
  
}

  
